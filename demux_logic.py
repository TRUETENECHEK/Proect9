from __future__ import annotations

import csv
import os
from itertools import product
from typing import Any, Iterable, Mapping, Sequence

import regex

from data_loader import load_fasta_dict
from matcher_core import (
    AMBIGUOUS_LAYOUT_STATUS,
    INVALID_LAYOUT_STATUS,
    MISSING_ADAPTERS_STATUS,
    PASS_STATUS,
    LayoutCandidate,
    LayoutDecision,
    filter_valid_pairs,
)
from reference_builder import ReferenceBundle, ReferenceVariant, build_reference_bundle
from sequence_matcher import reverse_complement


MISSING_BARCODE_STATUS = "MissingBarcode"
AMBIGUOUS_BARCODE_STATUS = "AmbiguousBarcode"
PRIMER_WARNING_STATUS = "PrimerWarning"
DISCARDED_STATUS = "Discarded"

RESULT_FIELDNAMES = [
    "read_id",
    "barcode_id",
    "sample_id",
    "score",
    "status",
    "trimmed_seq",
    "left_barcode_id",
    "right_barcode_id",
    "adapter_penalty",
    "barcode_penalty",
    "primer_penalty",
    "layout_penalty",
    "total_penalty",
    "trim_start",
    "trim_end",
]

_REFERENCE_BUNDLE_CACHE: dict[tuple[int | None, int | None, int | None], ReferenceBundle] = {}


def process_read(
    read_id,
    read_seq,
    adapters=None,
    barcodes=None,
    config=None,
    mapping_dict=None,
    max_errors=2,
    quality=None,
    reference_bundle=None,
):
    """
    Process one read using an adapter-first layout pipeline.

    The old call shape is still accepted. For the new architecture, pass a
    ReferenceBundle via ``reference_bundle`` or as ``adapters`` to avoid
    rebuilding references for every read.
    """
    config = config or {}
    bundle = _resolve_reference_bundle(reference_bundle or adapters, adapters, barcodes, config)
    result = _base_result(read_id, read_seq, quality)

    if bundle is None:
        return _finish_result(
            result,
            status=MISSING_ADAPTERS_STATUS,
            reason="ReferenceBundle is unavailable",
        )

    adapter_matches = _find_reference_matches(
        read_seq=read_seq,
        variants=_variants_for_role(bundle, "adapter"),
        max_errors=max_errors,
    )
    adapter_decision = _select_adapter_layout(adapter_matches, config, len(read_seq))
    if adapter_decision.status != PASS_STATUS:
        return _finish_from_layout_decision(result, adapter_decision)

    adapter_candidate = adapter_decision.selected_candidate
    full_layout_status, full_layout = _full_layout_for_adapter_candidate(adapter_candidate, config)
    if full_layout_status != PASS_STATUS:
        return _finish_result(result, status=full_layout_status)

    barcode_matches_by_site = _find_expected_barcodes(
        read_seq=read_seq,
        bundle=bundle,
        config=config,
        full_layout=full_layout,
        adapter_candidate=adapter_candidate,
        max_errors=_barcode_max_errors(config, max_errors),
    )
    missing_sites = [site for site, matches in barcode_matches_by_site.items() if not matches]
    if missing_sites:
        return _finish_result(
            result,
            status=MISSING_BARCODE_STATUS,
            adapter_candidate=adapter_candidate,
            reason=f"Missing barcode matches for: {', '.join(missing_sites)}",
        )

    pair_decisions = _score_barcode_pairs(
        adapter_candidate=adapter_candidate,
        barcode_matches_by_site=barcode_matches_by_site,
        config=config,
        read_length=len(read_seq),
    )
    if not pair_decisions:
        return _finish_result(
            result,
            status=INVALID_LAYOUT_STATUS,
            adapter_candidate=adapter_candidate,
            reason="No barcode pair produced a valid layout",
        )

    pair_decisions.sort(key=lambda item: item["sort_key"])
    best_pair = pair_decisions[0]
    if _barcode_pair_is_ambiguous(best_pair, pair_decisions[1:], config):
        return _finish_barcode_pair_result(
            result=result,
            status=AMBIGUOUS_BARCODE_STATUS,
            pair_info=best_pair,
            read_seq=read_seq,
            quality=quality,
            include_trim=False,
            reason="Best barcode pair is ambiguous",
        )

    trim_start, trim_end = _trim_bounds_from_candidate(best_pair["candidate"])
    if trim_start is None or trim_end is None or trim_start >= trim_end:
        return _finish_barcode_pair_result(
            result=result,
            status=INVALID_LAYOUT_STATUS,
            pair_info=best_pair,
            read_seq=read_seq,
            quality=quality,
            include_trim=False,
            reason="Barcode pair does not define a positive amplicon interval",
        )

    primer_qc = _run_primer_qc(
        read_seq=read_seq,
        trim_start=trim_start,
        trim_end=trim_end,
        config=config,
        default_max_errors=max_errors,
    )
    status = PRIMER_WARNING_STATUS if primer_qc["warning"] else PASS_STATUS

    result = _finish_barcode_pair_result(
        result=result,
        status=status,
        pair_info=best_pair,
        read_seq=read_seq,
        quality=quality,
        include_trim=True,
        primer_penalty=primer_qc["penalty"],
        trim_start=trim_start,
        trim_end=trim_end,
        reason=primer_qc["reason"],
    )

    sample_id, should_discard = _sample_id_for_barcode_pair(
        mapping_dict=mapping_dict,
        left_barcode_id=result["left_barcode_id"],
        right_barcode_id=result["right_barcode_id"],
    )
    result["sample_id"] = sample_id
    if should_discard:
        result["status"] = DISCARDED_STATUS
        result["sample_id"] = None
        result["trimmed_seq"] = ""
        result["trimmed_quality"] = None

    return result


def _resolve_reference_bundle(
    reference_bundle: Any,
    adapters: Any,
    barcodes: Any,
    config: Mapping[str, Any],
) -> ReferenceBundle | None:
    if isinstance(reference_bundle, ReferenceBundle):
        return reference_bundle
    if isinstance(adapters, ReferenceBundle):
        return adapters

    cache_key = (
        id(config) if config is not None else None,
        id(adapters) if adapters is not None else None,
        id(barcodes) if barcodes is not None else None,
    )
    if cache_key in _REFERENCE_BUNDLE_CACHE:
        return _REFERENCE_BUNDLE_CACHE[cache_key]

    if not config:
        return None

    try:
        bundle = build_reference_bundle(config, base_dir=".")
    except Exception:
        return None

    _REFERENCE_BUNDLE_CACHE[cache_key] = bundle
    return bundle


def _variants_for_role(bundle: ReferenceBundle, role: str) -> list[ReferenceVariant]:
    reference_type = "barcode_core" if role == "barcode" else role
    return [
        variant
        for variant in bundle.variants
        if variant.reference_type == reference_type and _usable_sequence(variant.sequence)
    ]


def _variants_for_barcode_site(bundle: ReferenceBundle, barcode_site: str) -> list[ReferenceVariant]:
    side = _side_from_site(barcode_site)
    return [
        variant
        for variant in _variants_for_role(bundle, "barcode")
        if variant.side == side
    ]


def _find_reference_matches(
    read_seq: str,
    variants: Sequence[ReferenceVariant],
    max_errors: int,
    window: tuple[int, int] | None = None,
) -> list[dict[str, Any]]:
    start_offset, search_seq = _windowed_sequence(read_seq, window)
    records: dict[tuple[Any, ...], dict[str, Any]] = {}

    for variant in variants:
        pattern = regex.compile(f"(?e)(?:{regex.escape(variant.sequence)}){{e<={max_errors}}}")
        for match in pattern.finditer(search_seq, overlapped=True):
            substitutions, insertions, deletions = match.fuzzy_counts
            start = start_offset + match.start()
            end = start_offset + match.end()
            record = _match_record_from_variant(
                variant=variant,
                start=start,
                end=end,
                errors=(substitutions, insertions, deletions),
            )
            records.setdefault(_record_identity(record), record)

    return sorted(records.values(), key=_record_sort_key)


def _match_record_from_variant(
    variant: ReferenceVariant,
    start: int,
    end: int,
    errors: tuple[int, int, int],
) -> dict[str, Any]:
    role = "barcode" if variant.reference_type == "barcode_core" else variant.reference_type
    site = f"{variant.side}_{role}"
    penalty = sum(errors)
    return {
        "site": site,
        "role": role,
        "side": variant.side,
        "start": start,
        "end": end,
        "errors": errors,
        "penalty": penalty,
        "ref_id": variant.source_id,
        "match_id": variant.source_id,
        "element_type": variant.reference_type,
        "reference_type": variant.reference_type,
        "sequence": variant.sequence,
        "transform": variant.transform,
        "strand": variant.sign,
        "sign": variant.sign,
        "direction": variant.direction,
    }


def _select_adapter_layout(
    adapter_matches: Sequence[Mapping[str, Any]],
    config: Mapping[str, Any],
    read_length: int,
) -> LayoutDecision:
    adapter_config = _adapter_only_config(config)
    return filter_valid_pairs(
        adapter_matches,
        [],
        adapter_config,
        read_length=read_length,
    )


def _adapter_only_config(config: Mapping[str, Any]) -> dict[str, Any]:
    copied_config = dict(config)
    layout_rules = dict(config.get("layout_rules", {}) or {})
    adapter_layouts: list[list[str]] = []

    for layout in layout_rules.get("valid_layouts", []) or []:
        adapter_layout = [
            str(site)
            for site in layout
            if _site_role(str(site), layout_rules) == "adapter"
        ]
        if len(adapter_layout) >= 2 and adapter_layout not in adapter_layouts:
            adapter_layouts.append(adapter_layout)

    layout_rules["valid_layouts"] = adapter_layouts
    layout_rules["allowed_distances"] = {}
    copied_config["layout_rules"] = layout_rules
    return copied_config


def _full_layout_for_adapter_candidate(
    adapter_candidate: LayoutCandidate | None,
    config: Mapping[str, Any],
) -> tuple[str, tuple[str, ...] | None]:
    if adapter_candidate is None:
        return INVALID_LAYOUT_STATUS, None

    adapter_layout = tuple(adapter_candidate.layout)
    layout_rules = config.get("layout_rules", {}) or {}
    matching_layouts = [
        tuple(str(site) for site in layout)
        for layout in layout_rules.get("valid_layouts", []) or []
        if tuple(
            str(site)
            for site in layout
            if _site_role(str(site), layout_rules) == "adapter"
        )
        == adapter_layout
    ]

    if len(matching_layouts) == 1:
        return PASS_STATUS, matching_layouts[0]
    if len(matching_layouts) > 1:
        return AMBIGUOUS_LAYOUT_STATUS, None
    return INVALID_LAYOUT_STATUS, None


def _find_expected_barcodes(
    read_seq: str,
    bundle: ReferenceBundle,
    config: Mapping[str, Any],
    full_layout: tuple[str, ...],
    adapter_candidate: LayoutCandidate,
    max_errors: int,
) -> dict[str, list[dict[str, Any]]]:
    adapter_by_site = {
        match.site: match
        for match in adapter_candidate.selected_matches
        if match.role == "adapter"
    }
    barcode_matches_by_site: dict[str, list[dict[str, Any]]] = {}

    for barcode_site in _barcode_sites(full_layout, config):
        anchor = _barcode_anchor(full_layout, barcode_site, adapter_by_site, config)
        if anchor is None:
            barcode_matches_by_site[barcode_site] = []
            continue

        variants = _variants_for_barcode_site(bundle, barcode_site)
        window = _barcode_search_window(
            read_length=len(read_seq),
            variants=variants,
            anchor=anchor,
            config=config,
            max_errors=max_errors,
        )
        matches = _find_reference_matches(
            read_seq=read_seq,
            variants=variants,
            max_errors=max_errors,
            window=window,
        )
        barcode_matches_by_site[barcode_site] = [
            match
            for match in matches
            if _match_satisfies_anchor(match, anchor)
        ]

    return barcode_matches_by_site


def _barcode_sites(full_layout: Sequence[str], config: Mapping[str, Any]) -> list[str]:
    layout_rules = config.get("layout_rules", {}) or {}
    return [
        str(site)
        for site in full_layout
        if _site_role(str(site), layout_rules) == "barcode"
    ]


def _barcode_anchor(
    full_layout: tuple[str, ...],
    barcode_site: str,
    adapter_by_site: Mapping[str, Any],
    config: Mapping[str, Any],
) -> dict[str, Any] | None:
    layout_rules = config.get("layout_rules", {}) or {}
    barcode_index = full_layout.index(barcode_site)
    barcode_side = _side_from_site(barcode_site)

    for adapter_site, adapter_match in adapter_by_site.items():
        if _side_from_site(adapter_site) != barcode_side:
            continue
        adapter_index = full_layout.index(adapter_site)
        direction = "downstream" if adapter_index < barcode_index else "upstream"
        min_distance, max_distance = _distance_limits(
            first_site=adapter_site if direction == "downstream" else barcode_site,
            second_site=barcode_site if direction == "downstream" else adapter_site,
            layout_rules=layout_rules,
            default_max=_default_barcode_window(config),
        )
        return {
            "adapter": adapter_match,
            "direction": direction,
            "min_distance": min_distance,
            "max_distance": max_distance,
        }

    return None


def _barcode_search_window(
    read_length: int,
    variants: Sequence[ReferenceVariant],
    anchor: Mapping[str, Any],
    config: Mapping[str, Any],
    max_errors: int,
) -> tuple[int, int]:
    adapter = anchor["adapter"]
    min_distance = int(anchor["min_distance"])
    max_distance = int(anchor["max_distance"])
    max_barcode_length = max((len(variant.sequence) for variant in variants), default=0)
    slack = max_errors + int(config.get("barcode_window_slack", 2))

    if anchor["direction"] == "downstream":
        start = adapter.end + min_distance
        end = adapter.end + max_distance + max_barcode_length + slack
    else:
        start = adapter.start - max_distance - max_barcode_length - slack
        end = adapter.start - min_distance

    return max(0, start), min(read_length, end)


def _match_satisfies_anchor(match: Mapping[str, Any], anchor: Mapping[str, Any]) -> bool:
    adapter = anchor["adapter"]
    min_distance = int(anchor["min_distance"])
    max_distance = int(anchor["max_distance"])

    if anchor["direction"] == "downstream":
        distance = int(match["start"]) - adapter.end
    else:
        distance = adapter.start - int(match["end"])

    return min_distance <= distance <= max_distance


def _distance_limits(
    first_site: str,
    second_site: str,
    layout_rules: Mapping[str, Any],
    default_max: int,
) -> tuple[int, int]:
    for rule in _distance_rules(layout_rules):
        if str(rule.get("from")) == first_site and str(rule.get("to")) == second_site:
            return int(rule.get("min", 0)), int(rule.get("max", default_max))

    first_side = _side_from_site(first_site)
    second_side = _side_from_site(second_site)
    for rule in _distance_rules(layout_rules):
        from_site = str(rule.get("from"))
        to_site = str(rule.get("to"))
        if {
            _side_from_site(from_site),
            _side_from_site(to_site),
        } == {first_side, second_side} and first_side == second_side:
            return int(rule.get("min", 0)), int(rule.get("max", default_max))

    return 0, default_max


def _distance_rules(layout_rules: Mapping[str, Any]) -> Iterable[dict[str, Any]]:
    raw_rules = layout_rules.get("allowed_distances", {}) or {}
    if isinstance(raw_rules, Mapping):
        for name, raw_rule in raw_rules.items():
            if not isinstance(raw_rule, Mapping):
                continue
            rule = dict(raw_rule)
            if ("from" not in rule or "to" not in rule) and "_to_" in str(name):
                from_site, to_site = str(name).split("_to_", 1)
                rule.setdefault("from", from_site)
                rule.setdefault("to", to_site)
            yield rule
        return

    if isinstance(raw_rules, Sequence) and not isinstance(raw_rules, (str, bytes)):
        for raw_rule in raw_rules:
            if isinstance(raw_rule, Mapping):
                yield dict(raw_rule)


def _score_barcode_pairs(
    adapter_candidate: LayoutCandidate,
    barcode_matches_by_site: Mapping[str, Sequence[Mapping[str, Any]]],
    config: Mapping[str, Any],
    read_length: int,
) -> list[dict[str, Any]]:
    sites = list(barcode_matches_by_site)
    decisions: list[dict[str, Any]] = []

    for barcode_pair in product(*(barcode_matches_by_site[site] for site in sites)):
        layout_decision = filter_valid_pairs(
            adapter_candidate.selected_matches,
            list(barcode_pair),
            config,
            read_length=read_length,
        )
        if layout_decision.selected_candidate is None:
            continue

        candidate = layout_decision.selected_candidate
        if layout_decision.status not in (PASS_STATUS, AMBIGUOUS_LAYOUT_STATUS):
            continue

        decisions.append(
            {
                "decision": layout_decision,
                "candidate": candidate,
                "barcode_pair": barcode_pair,
                "sort_key": (
                    candidate.total_penalty,
                    _barcode_pair_identity(candidate),
                    _barcode_pair_positions(candidate),
                ),
            }
        )

    return decisions


def _barcode_pair_is_ambiguous(
    best_pair: Mapping[str, Any],
    other_pairs: Sequence[Mapping[str, Any]],
    config: Mapping[str, Any],
) -> bool:
    if not other_pairs:
        return False

    layout_rules = config.get("layout_rules", {}) or {}
    ambiguity_margin = max(0.0, float(layout_rules.get("ambiguity_margin", 0)))
    best_penalty = best_pair["candidate"].total_penalty
    next_penalty = other_pairs[0]["candidate"].total_penalty
    delta = next_penalty - best_penalty
    if delta != 0 and delta >= ambiguity_margin:
        return False

    return _barcode_pair_identity(best_pair["candidate"]) != _barcode_pair_identity(
        other_pairs[0]["candidate"]
    ) or _barcode_pair_positions(best_pair["candidate"]) != _barcode_pair_positions(
        other_pairs[0]["candidate"]
    )


def _finish_from_layout_decision(
    result: dict[str, Any],
    decision: LayoutDecision,
) -> dict[str, Any]:
    return _finish_result(
        result,
        status=decision.status,
        adapter_candidate=decision.selected_candidate,
        reason=decision.reason,
    )


def _finish_barcode_pair_result(
    result: dict[str, Any],
    status: str,
    pair_info: Mapping[str, Any],
    read_seq: str,
    quality: Any,
    include_trim: bool,
    primer_penalty: float = 0.0,
    trim_start: int | None = None,
    trim_end: int | None = None,
    reason: str = "",
) -> dict[str, Any]:
    candidate = pair_info["candidate"]
    left_barcode = _match_by_site(candidate, "left_barcode")
    right_barcode = _match_by_site(candidate, "right_barcode")
    penalties = _penalty_summary(candidate, primer_penalty)

    result.update(
        {
            "status": status,
            "barcode_id": _barcode_id(left_barcode, right_barcode),
            "left_barcode_id": _match_ref_id(left_barcode),
            "right_barcode_id": _match_ref_id(right_barcode),
            "adapter_penalty": penalties["adapter_penalty"],
            "barcode_penalty": penalties["barcode_penalty"],
            "primer_penalty": penalties["primer_penalty"],
            "layout_penalty": penalties["layout_penalty"],
            "total_penalty": penalties["total_penalty"],
            "score": -penalties["total_penalty"],
            "trim_start": trim_start,
            "trim_end": trim_end,
            "diagnostic_message": reason,
        }
    )

    if include_trim and trim_start is not None and trim_end is not None:
        result["trimmed_seq"] = read_seq[trim_start:trim_end]
        result["trimmed_quality"] = quality[trim_start:trim_end] if quality is not None else None

    return result


def _finish_result(
    result: dict[str, Any],
    status: str,
    adapter_candidate: LayoutCandidate | None = None,
    reason: str = "",
) -> dict[str, Any]:
    if adapter_candidate is not None:
        penalties = _penalty_summary(adapter_candidate, primer_penalty=0.0)
        result.update(
            {
                "adapter_penalty": penalties["adapter_penalty"],
                "layout_penalty": penalties["layout_penalty"],
                "total_penalty": penalties["total_penalty"],
                "score": -penalties["total_penalty"],
            }
        )

    result["status"] = status
    result["diagnostic_message"] = reason
    return result


def _base_result(read_id: str, read_seq: str, quality: Any) -> dict[str, Any]:
    return {
        "read_id": read_id,
        "barcode_id": "Unknown",
        "sample_id": None,
        "score": 0.0,
        "status": DISCARDED_STATUS,
        "trimmed_seq": "",
        "trimmed_quality": None,
        "original_seq": read_seq,
        "original_quality": quality,
        "left_barcode_id": None,
        "right_barcode_id": None,
        "adapter_penalty": 0.0,
        "barcode_penalty": 0.0,
        "primer_penalty": 0.0,
        "layout_penalty": 0.0,
        "total_penalty": 0.0,
        "trim_start": None,
        "trim_end": None,
        "diagnostic_message": "",
    }


def _penalty_summary(candidate: LayoutCandidate, primer_penalty: float) -> dict[str, float]:
    adapter_penalty = sum(
        match.edit_penalty for match in candidate.selected_matches if match.role == "adapter"
    )
    barcode_penalty = sum(
        match.edit_penalty for match in candidate.selected_matches if match.role == "barcode"
    )
    layout_penalty = sum(
        value
        for key, value in candidate.penalty_breakdown.items()
        if key != "edit"
    )
    total_penalty = adapter_penalty + barcode_penalty + layout_penalty + primer_penalty

    return {
        "adapter_penalty": float(adapter_penalty),
        "barcode_penalty": float(barcode_penalty),
        "primer_penalty": float(primer_penalty),
        "layout_penalty": float(layout_penalty),
        "total_penalty": float(total_penalty),
    }


def _trim_bounds_from_candidate(candidate: LayoutCandidate) -> tuple[int | None, int | None]:
    barcode_matches = sorted(
        [match for match in candidate.selected_matches if match.role == "barcode"],
        key=lambda match: (match.start, match.end),
    )
    if len(barcode_matches) < 2:
        return None, None
    return barcode_matches[0].end, barcode_matches[-1].start


def _run_primer_qc(
    read_seq: str,
    trim_start: int,
    trim_end: int,
    config: Mapping[str, Any],
    default_max_errors: int,
) -> dict[str, Any]:
    primer_references, primer_config = _primer_references(config)
    if not primer_references:
        return {"warning": False, "penalty": 0.0, "reason": ""}

    max_errors = int(primer_config.get("max_errors", default_max_errors))
    matches = _find_simple_reference_matches(
        read_seq=read_seq,
        references=primer_references,
        max_errors=max_errors,
        window=(trim_start, trim_end),
    )
    if not matches:
        missing_penalty = float(primer_config.get("missing_penalty", 0))
        return {
            "warning": True,
            "penalty": missing_penalty,
            "reason": "Configured primer QC anchors were not found inside the amplicon",
        }

    best_by_ref: dict[str, Mapping[str, Any]] = {}
    for match in matches:
        ref_id = str(match["ref_id"])
        existing = best_by_ref.get(ref_id)
        if existing is None or _record_sort_key(match) < _record_sort_key(existing):
            best_by_ref[ref_id] = match

    required = str(primer_config.get("required", "any")).lower()
    missing_refs = set(primer_references) - set(best_by_ref)
    warning = required == "all" and bool(missing_refs)
    penalty = sum(
        _edit_penalty_from_errors(match["errors"], config)
        for match in best_by_ref.values()
    )
    return {
        "warning": warning,
        "penalty": penalty,
        "reason": (
            f"Missing primer QC anchors: {', '.join(sorted(missing_refs))}"
            if warning
            else ""
        ),
    }


def _primer_references(config: Mapping[str, Any]) -> tuple[dict[str, str], Mapping[str, Any]]:
    primer_config = (
        config.get("primer_qc")
        or config.get("primers")
        or config.get("primer")
        or {}
    )
    if not primer_config:
        return {}, {}

    if isinstance(primer_config, Mapping):
        if primer_config.get("enabled") is False:
            return {}, primer_config

        raw_sequences = primer_config.get("sequences")
        if isinstance(raw_sequences, Mapping):
            return {str(key): str(value) for key, value in raw_sequences.items()}, primer_config

        fasta_path = primer_config.get("fasta_path") or primer_config.get("path")
        if fasta_path and os.path.exists(str(fasta_path)):
            return load_fasta_dict(str(fasta_path)), primer_config

    if isinstance(primer_config, Mapping):
        sequence_items = {
            str(key): str(value)
            for key, value in primer_config.items()
            if isinstance(value, str) and set(value.upper()) <= set("ACGTN")
        }
        if sequence_items:
            return sequence_items, primer_config

    return {}, primer_config if isinstance(primer_config, Mapping) else {}


def _find_simple_reference_matches(
    read_seq: str,
    references: Mapping[str, str],
    max_errors: int,
    window: tuple[int, int] | None = None,
) -> list[dict[str, Any]]:
    start_offset, search_seq = _windowed_sequence(read_seq, window)
    records: list[dict[str, Any]] = []

    for ref_id, sequence in references.items():
        for transform, search_sequence in (
            ("forward", sequence),
            ("reverse_complement", reverse_complement(sequence)),
        ):
            if not _usable_sequence(search_sequence):
                continue
            pattern = regex.compile(f"(?e)(?:{regex.escape(search_sequence)}){{e<={max_errors}}}")
            for match in pattern.finditer(search_seq, overlapped=True):
                errors = match.fuzzy_counts
                records.append(
                    {
                        "site": "primer",
                        "role": "primer",
                        "side": None,
                        "start": start_offset + match.start(),
                        "end": start_offset + match.end(),
                        "errors": errors,
                        "penalty": sum(errors),
                        "ref_id": str(ref_id),
                        "match_id": str(ref_id),
                        "element_type": "primer",
                        "transform": transform,
                        "strand": "-" if transform == "reverse_complement" else "+",
                        "direction": "reverse" if transform == "reverse_complement" else "forward",
                    }
                )

    return sorted(records, key=_record_sort_key)


def _sample_id_for_barcode_pair(
    mapping_dict: Mapping[Any, Any] | None,
    left_barcode_id: str | None,
    right_barcode_id: str | None,
) -> tuple[str | None, bool]:
    if not left_barcode_id or not right_barcode_id:
        return None, False

    barcode_id = f"{left_barcode_id}_{right_barcode_id}"
    if not mapping_dict:
        return barcode_id, False

    if barcode_id in mapping_dict:
        return str(mapping_dict[barcode_id]), False

    tuple_key = (left_barcode_id, right_barcode_id)
    if tuple_key in mapping_dict:
        return str(mapping_dict[tuple_key]), False

    return None, True


def _edit_penalty_from_errors(errors: Sequence[int], config: Mapping[str, Any]) -> float:
    layout_scoring = config.get("layout_scoring", {}) or {}
    edit_penalties = layout_scoring.get("edit_penalties", {}) or {}
    scoring_weights = config.get("scoring_weights", {}) or {}
    mismatch = abs(float(edit_penalties.get("mismatch", scoring_weights.get("mismatch", 1))))
    insertion = abs(float(edit_penalties.get("insertion", scoring_weights.get("insertion", 1))))
    deletion = abs(float(edit_penalties.get("deletion", scoring_weights.get("deletion", 1))))
    substitutions, insertions, deletions = errors
    return substitutions * mismatch + insertions * insertion + deletions * deletion


def _site_role(site: str, layout_rules: Mapping[str, Any]) -> str | None:
    metadata = (layout_rules.get("match_sites", {}) or {}).get(site, {}) or {}
    role = metadata.get("role") or metadata.get("kind")
    if role:
        return str(role)
    if "adapter" in site:
        return "adapter"
    if "barcode" in site:
        return "barcode"
    if "primer" in site:
        return "primer"
    return None


def _side_from_site(site: str | None) -> str | None:
    if site is None:
        return None
    site = str(site)
    if site.startswith("left_") or site == "left":
        return "left"
    if site.startswith("right_") or site == "right":
        return "right"
    return None


def _barcode_max_errors(config: Mapping[str, Any], default_max_errors: int) -> int:
    return int(
        config.get(
            "barcode_max_errors",
            config.get("max_barcode_errors", min(default_max_errors, 1)),
        )
    )


def _default_barcode_window(config: Mapping[str, Any]) -> int:
    layout_rules = config.get("layout_rules", {}) or {}
    return int(config.get("barcode_window", layout_rules.get("barcode_window", 30)))


def _windowed_sequence(read_seq: str, window: tuple[int, int] | None) -> tuple[int, str]:
    if window is None:
        return 0, read_seq
    start, end = window
    start = max(0, int(start))
    end = min(len(read_seq), int(end))
    if start >= end:
        return start, ""
    return start, read_seq[start:end]


def _usable_sequence(sequence: str | None) -> bool:
    if not sequence:
        return False
    return bool(set(str(sequence).upper()) - {"N"})


def _record_identity(record: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        record.get("site"),
        record.get("ref_id"),
        record.get("start"),
        record.get("end"),
        record.get("transform"),
        record.get("errors"),
    )


def _record_sort_key(record: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        record.get("penalty", 0),
        record.get("start", 0),
        record.get("end", 0),
        record.get("site") or "",
        record.get("ref_id") or "",
        record.get("transform") or "",
    )


def _match_by_site(candidate: LayoutCandidate, site: str):
    for match in candidate.selected_matches:
        if match.site == site:
            return match
    return None


def _match_ref_id(match: Any) -> str | None:
    if match is None:
        return None
    return match.ref_id or match.site


def _barcode_id(left_barcode: Any, right_barcode: Any) -> str:
    left_id = _match_ref_id(left_barcode)
    right_id = _match_ref_id(right_barcode)
    if not left_id or not right_id:
        return "Unknown"
    return f"{left_id}_{right_id}"


def _barcode_pair_identity(candidate: LayoutCandidate) -> tuple[str | None, str | None]:
    return (
        _match_ref_id(_match_by_site(candidate, "left_barcode")),
        _match_ref_id(_match_by_site(candidate, "right_barcode")),
    )


def _barcode_pair_positions(candidate: LayoutCandidate) -> tuple[tuple[str, int, int], ...]:
    return tuple(
        sorted(
            (match.site, match.start, match.end)
            for match in candidate.selected_matches
            if match.role == "barcode"
        )
    )


def save_process_results(results, output_file="demultiplex_results.csv"):
    with open(output_file, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=RESULT_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for res in results:
            writer.writerow(res)
