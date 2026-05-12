from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations, product
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


PASS_STATUS = "Pass"
MISSING_ADAPTERS_STATUS = "MissingAdapters"
INVALID_LAYOUT_STATUS = "InvalidLayout"
AMBIGUOUS_LAYOUT_STATUS = "AmbiguousLayout"


@dataclass(frozen=True)
class CandidateMatch:
    site: str
    role: str
    side: str | None
    start: int
    end: int
    strand: str | None
    direction: str | None
    record: Any
    ref_id: str | None = None
    element_type: str | None = None
    transform: str | None = None
    edit_penalty: float = 0.0
    source_key: str | None = None


@dataclass(frozen=True)
class LayoutCandidate:
    adapter_pair: tuple[CandidateMatch, CandidateMatch]
    selected_matches: tuple[CandidateMatch, ...]
    ordered_sites: tuple[str, ...]
    layout: tuple[str, ...]
    total_penalty: float
    penalty_breakdown: Mapping[str, float] = field(default_factory=dict)
    reason: str = ""

    @property
    def matches_by_site(self) -> dict[str, Any]:
        return {match.site: match.record for match in self.selected_matches}

    def as_legacy_elements(self) -> dict[str, Any]:
        return {match.site: _legacy_record(match) for match in self.selected_matches}


@dataclass(frozen=True)
class LayoutDecision:
    status: str
    selected_candidate: LayoutCandidate | None
    all_candidates_count: int
    reason: str

    def __iter__(self):
        """
        Backward compatibility for old callers:
        filter_status, filtered_elements = filter_valid_pairs(...)
        """
        if self.status == PASS_STATUS and self.selected_candidate is not None:
            yield PASS_STATUS
            yield self.selected_candidate.as_legacy_elements()
            return

        yield "Discarded"
        yield {}


def filter_valid_pairs(
    adapter_matches: Any,
    barcode_primer_matches: Any = None,
    config: Mapping[str, Any] | None = None,
    read_length: int | None = None,
) -> LayoutDecision:
    """
    Build and score layout candidates from adapter matches plus barcode/primer matches.

    New API:
        filter_valid_pairs(adapter_match_records, barcode_or_primer_records, config)

    The old dict API is still accepted for existing demux_logic callers. The returned
    LayoutDecision can be unpacked as the old ("Pass"/"Discarded", elements) tuple.
    """
    if config is None and _looks_like_config(barcode_primer_matches):
        config = barcode_primer_matches
        barcode_primer_matches = None

    config = _load_default_config() if config is None else config
    layout_rules = _layout_rules(config)
    scoring = _layout_scoring(config)
    read_length = _configured_read_length(read_length, config, layout_rules)

    adapters, barcode_primers = _normalize_inputs(
        adapter_matches=adapter_matches,
        barcode_primer_matches=barcode_primer_matches,
        layout_rules=layout_rules,
        scoring=scoring,
    )

    if len(adapters) < 2:
        return LayoutDecision(
            status=MISSING_ADAPTERS_STATUS,
            selected_candidate=None,
            all_candidates_count=0,
            reason=f"Expected at least 2 adapter matches, got {len(adapters)}",
        )

    adapter_pairs = _build_allowed_adapter_pairs(adapters, layout_rules)
    if not adapter_pairs:
        return LayoutDecision(
            status=INVALID_LAYOUT_STATUS,
            selected_candidate=None,
            all_candidates_count=0,
            reason="No adapter pairs matched allowed_adapter_pairs",
        )

    candidates = _build_layout_candidates(
        adapter_pairs=adapter_pairs,
        barcode_primer_matches=barcode_primers,
        layout_rules=layout_rules,
        scoring=scoring,
        read_length=read_length,
    )

    if not candidates:
        return LayoutDecision(
            status=INVALID_LAYOUT_STATUS,
            selected_candidate=None,
            all_candidates_count=0,
            reason="No adapter pair could be completed into a valid layout",
        )

    candidates = sorted(candidates, key=_candidate_sort_key)
    selected = candidates[0]
    ambiguity_margin = max(0.0, float(layout_rules.get("ambiguity_margin", 0)))

    if len(candidates) > 1:
        delta = candidates[1].total_penalty - selected.total_penalty
        if delta == 0 or delta < ambiguity_margin:
            return LayoutDecision(
                status=AMBIGUOUS_LAYOUT_STATUS,
                selected_candidate=selected,
                all_candidates_count=len(candidates),
                reason=(
                    "Best layout is ambiguous: "
                    f"next candidate is {delta:g} penalty units away"
                ),
            )

    return LayoutDecision(
        status=PASS_STATUS,
        selected_candidate=selected,
        all_candidates_count=len(candidates),
        reason="Selected minimum-penalty layout",
    )


def _normalize_inputs(
    adapter_matches: Any,
    barcode_primer_matches: Any,
    layout_rules: Mapping[str, Any],
    scoring: Mapping[str, Any],
) -> tuple[list[CandidateMatch], list[CandidateMatch]]:
    if isinstance(adapter_matches, Mapping) and barcode_primer_matches is None:
        normalized = [
            match
            for match in _iter_mapping_records(adapter_matches, layout_rules, scoring)
            if match is not None
        ]
        return (
            [match for match in normalized if match.role == "adapter"],
            [match for match in normalized if match.role != "adapter"],
        )

    adapters = [
        match
        for match in _iter_records(
            adapter_matches,
            layout_rules=layout_rules,
            scoring=scoring,
            default_role="adapter",
        )
        if match is not None and match.role == "adapter"
    ]
    barcode_primers = [
        match
        for match in _iter_records(
            barcode_primer_matches,
            layout_rules=layout_rules,
            scoring=scoring,
        )
        if match is not None and match.role != "adapter"
    ]
    return adapters, barcode_primers


def _iter_mapping_records(
    records_by_site: Mapping[str, Any],
    layout_rules: Mapping[str, Any],
    scoring: Mapping[str, Any],
) -> Iterable[CandidateMatch | None]:
    for site, value in records_by_site.items():
        if value is None:
            continue
        if _is_record_collection(value):
            for record in value:
                yield _to_candidate_match(
                    record,
                    layout_rules=layout_rules,
                    scoring=scoring,
                    forced_site=str(site),
                )
            continue
        yield _to_candidate_match(
            value,
            layout_rules=layout_rules,
            scoring=scoring,
            forced_site=str(site),
        )


def _iter_records(
    records: Any,
    layout_rules: Mapping[str, Any],
    scoring: Mapping[str, Any],
    default_role: str | None = None,
) -> Iterable[CandidateMatch | None]:
    if records is None:
        return
    if isinstance(records, Mapping):
        yield from _iter_mapping_records(records, layout_rules, scoring)
        return
    if not _is_record_collection(records):
        yield _to_candidate_match(
            records,
            layout_rules=layout_rules,
            scoring=scoring,
            default_role=default_role,
        )
        return

    for record in records:
        yield _to_candidate_match(
            record,
            layout_rules=layout_rules,
            scoring=scoring,
            default_role=default_role,
        )


def _to_candidate_match(
    record: Any,
    layout_rules: Mapping[str, Any],
    scoring: Mapping[str, Any],
    forced_site: str | None = None,
    default_role: str | None = None,
) -> CandidateMatch | None:
    if isinstance(record, CandidateMatch):
        return record

    start = _field(record, "start")
    end = _field(record, "end")
    if start is None or end is None:
        return None

    metadata = _match_metadata(record, layout_rules, forced_site, default_role)
    site = metadata.get("site")
    role = metadata.get("role")
    if not site or not role:
        return None

    return CandidateMatch(
        site=str(site),
        role=str(role),
        side=_optional_str(metadata.get("side")),
        start=int(start),
        end=int(end),
        strand=_optional_str(_field(record, "strand") or _field(record, "sign")),
        direction=_optional_str(_field(record, "direction")),
        record=record,
        ref_id=_optional_str(
            _field(record, "ref_id")
            or _field(record, "match_id")
            or _field(record, "source_id")
            or _field(record, "id")
            or _field(record, "name")
        ),
        element_type=_optional_str(
            _field(record, "element_type")
            or _field(record, "reference_type")
            or _field(record, "type")
        ),
        transform=_optional_str(_field(record, "transform")),
        edit_penalty=_edit_penalty(record, scoring),
        source_key=forced_site,
    )


def _match_metadata(
    record: Any,
    layout_rules: Mapping[str, Any],
    forced_site: str | None,
    default_role: str | None,
) -> dict[str, Any]:
    configured_sites = layout_rules.get("match_sites", {}) or {}
    identity_values = [
        forced_site,
        _field(record, "site"),
        _field(record, "layout_site"),
        _field(record, "element_type"),
        _field(record, "reference_type"),
        _field(record, "ref_id"),
        _field(record, "match_id"),
        _field(record, "source_id"),
        _field(record, "id"),
        _field(record, "name"),
    ]

    matched_site, matched_meta = _lookup_configured_site(identity_values, configured_sites)
    site = forced_site or matched_site or _first_present(identity_values)
    role = (
        matched_meta.get("role")
        or matched_meta.get("kind")
        or _field(record, "role")
        or _field(record, "reference_type")
        or default_role
        or _infer_role(site)
    )
    side = matched_meta.get("side") or _field(record, "side") or _metadata_field(record, "side")

    if side is None:
        side = _infer_side_from_site(site, configured_sites)

    return {"site": site, "role": role, "side": side}


def _lookup_configured_site(
    values: Sequence[Any],
    configured_sites: Mapping[str, Any],
) -> tuple[str | None, Mapping[str, Any]]:
    normalized_values = {_normalize(value) for value in values if value is not None}
    for site, raw_metadata in configured_sites.items():
        metadata = raw_metadata or {}
        aliases = {site, *(metadata.get("aliases", []) or [])}
        normalized_aliases = {_normalize(alias) for alias in aliases}
        if normalized_values & normalized_aliases:
            return str(site), metadata
    return None, {}


def _build_allowed_adapter_pairs(
    adapters: Sequence[CandidateMatch],
    layout_rules: Mapping[str, Any],
) -> list[tuple[CandidateMatch, CandidateMatch]]:
    adapter_pairs: list[tuple[CandidateMatch, CandidateMatch]] = []

    for left, right in combinations(adapters, 2):
        first, second = sorted((left, right), key=_match_sort_key)
        if first.side is not None and first.side == second.side:
            continue
        if _adapter_pair_allowed(first, second, layout_rules):
            adapter_pairs.append((first, second))

    return adapter_pairs


def _adapter_pair_allowed(
    first: CandidateMatch,
    second: CandidateMatch,
    layout_rules: Mapping[str, Any],
) -> bool:
    allowed_pairs = layout_rules.get("allowed_adapter_pairs", []) or []
    if not allowed_pairs:
        return True

    return any(_adapter_pair_rule_matches(rule, first, second) for rule in allowed_pairs)


def _adapter_pair_rule_matches(
    rule: Any,
    first: CandidateMatch,
    second: CandidateMatch,
) -> bool:
    if isinstance(rule, Sequence) and not isinstance(rule, (str, bytes, Mapping)):
        if len(rule) != 2:
            return False
        first_rule = {"side": rule[0]}
        second_rule = {"side": rule[1]}
        return _match_allowed_values(first, first_rule) and _match_allowed_values(second, second_rule)

    if not isinstance(rule, Mapping):
        return False

    first_rule = rule.get("first") or rule.get("upstream") or rule.get("a")
    second_rule = rule.get("second") or rule.get("downstream") or rule.get("b")
    if first_rule is None or second_rule is None:
        return False

    if _match_allowed_values(first, first_rule) and _match_allowed_values(second, second_rule):
        return True

    if rule.get("unordered", False):
        return _match_allowed_values(first, second_rule) and _match_allowed_values(second, first_rule)

    return False


def _match_allowed_values(match: CandidateMatch, rule: Any) -> bool:
    if isinstance(rule, str):
        rule = {"side": rule}
    if not isinstance(rule, Mapping):
        return False

    for key, expected in rule.items():
        if key == "unordered":
            continue
        actual = _candidate_value(match, key)
        if not _value_allowed(actual, expected):
            return False
    return True


def _build_layout_candidates(
    adapter_pairs: Sequence[tuple[CandidateMatch, CandidateMatch]],
    barcode_primer_matches: Sequence[CandidateMatch],
    layout_rules: Mapping[str, Any],
    scoring: Mapping[str, Any],
    read_length: int | None,
) -> list[LayoutCandidate]:
    valid_layouts = [tuple(layout) for layout in layout_rules.get("valid_layouts", []) or []]
    all_matches = _unique_matches(
        [match for pair in adapter_pairs for match in pair] + list(barcode_primer_matches)
    )
    site_counts = _site_counts(all_matches)
    candidates: list[LayoutCandidate] = []

    if not valid_layouts:
        for pair in adapter_pairs:
            selected = tuple(sorted(pair, key=_match_sort_key))
            penalty, breakdown = _score_candidate(
                selected=selected,
                layout=tuple(match.site for match in selected),
                all_site_counts=site_counts,
                layout_rules=layout_rules,
                scoring=scoring,
                read_length=read_length,
            )
            candidates.append(
                LayoutCandidate(
                    adapter_pair=pair,
                    selected_matches=selected,
                    ordered_sites=tuple(match.site for match in selected),
                    layout=tuple(match.site for match in selected),
                    total_penalty=penalty,
                    penalty_breakdown=breakdown,
                )
            )
        return candidates

    barcode_primers_by_site = _group_by_site(barcode_primer_matches)
    candidate_by_identity: dict[tuple[Any, ...], LayoutCandidate] = {}
    for pair in adapter_pairs:
        adapter_sites = {match.site for match in pair}
        for layout in valid_layouts:
            if not adapter_sites.issubset(set(layout)):
                continue

            layout_adapter_sites = {
                site
                for site in layout
                if _configured_role(site, layout_rules) == "adapter" or site in adapter_sites
            }
            if layout_adapter_sites != adapter_sites:
                continue

            choices: list[Sequence[CandidateMatch]] = []
            missing_required_site = False
            for site in layout:
                if site in adapter_sites:
                    continue
                site_matches = barcode_primers_by_site.get(site, [])
                if not site_matches:
                    missing_required_site = True
                    break
                choices.append(site_matches)
            if missing_required_site:
                continue

            for selected_non_adapters in product(*choices) if choices else [()]:
                selected = tuple(sorted((*pair, *selected_non_adapters), key=_match_sort_key))
                penalty, breakdown = _score_candidate(
                    selected=selected,
                    layout=layout,
                    all_site_counts=site_counts,
                    layout_rules=layout_rules,
                    scoring=scoring,
                    read_length=read_length,
                )
                candidate = LayoutCandidate(
                    adapter_pair=pair,
                    selected_matches=selected,
                    ordered_sites=tuple(match.site for match in selected),
                    layout=layout,
                    total_penalty=penalty,
                    penalty_breakdown=breakdown,
                )
                key = _selected_identity(selected)
                existing = candidate_by_identity.get(key)
                if existing is None or _candidate_sort_key(candidate) < _candidate_sort_key(existing):
                    candidate_by_identity[key] = candidate

    return list(candidate_by_identity.values())


def _score_candidate(
    selected: Sequence[CandidateMatch],
    layout: tuple[str, ...],
    all_site_counts: Mapping[str, int],
    layout_rules: Mapping[str, Any],
    scoring: Mapping[str, Any],
    read_length: int | None,
) -> tuple[float, dict[str, float]]:
    breakdown = {
        "edit": sum(match.edit_penalty for match in selected),
        "wrong_order": 0.0,
        "distance": 0.0,
        "extra_flank": 0.0,
        "duplicate_site": 0.0,
        "repeat_primer": 0.0,
    }

    ordered_sites = tuple(match.site for match in sorted(selected, key=_match_sort_key))
    if ordered_sites != layout:
        breakdown["wrong_order"] += _penalty_value(scoring, "wrong_order_penalty")

    selected_by_site = {match.site: match for match in selected}
    breakdown["distance"] += _distance_penalty(selected_by_site, layout_rules, scoring)

    if read_length is not None and selected:
        ordered = sorted(selected, key=_match_sort_key)
        flank_bases = max(0, ordered[0].start) + max(0, read_length - ordered[-1].end)
        breakdown["extra_flank"] += flank_bases * _penalty_value(scoring, "extra_flank_penalty")

    duplicate_sites = sum(max(0, all_site_counts.get(match.site, 0) - 1) for match in selected)
    breakdown["duplicate_site"] += duplicate_sites * _penalty_value(
        scoring,
        "duplicate_site_penalty",
    )

    repeated_primers = _repeat_primer_count(selected)
    breakdown["repeat_primer"] += repeated_primers * _penalty_value(
        scoring,
        "repeat_primer_penalty",
    )

    return sum(breakdown.values()), breakdown


def _distance_penalty(
    selected_by_site: Mapping[str, CandidateMatch],
    layout_rules: Mapping[str, Any],
    scoring: Mapping[str, Any],
) -> float:
    total = 0.0
    per_base = _penalty_value(scoring, "distance_penalty_per_base")
    if per_base == 0:
        return 0.0

    for rule in _distance_rules(layout_rules):
        from_site = rule.get("from")
        to_site = rule.get("to")
        if not from_site or not to_site:
            continue
        left = selected_by_site.get(str(from_site))
        right = selected_by_site.get(str(to_site))
        if left is None or right is None:
            continue

        distance = right.start - left.end
        min_distance = rule.get("min")
        max_distance = rule.get("max")
        if min_distance is not None and distance < int(min_distance):
            total += (int(min_distance) - distance) * per_base
        if max_distance is not None and distance > int(max_distance):
            total += (distance - int(max_distance)) * per_base

    return total


def _distance_rules(layout_rules: Mapping[str, Any]) -> Iterable[dict[str, Any]]:
    raw_rules = layout_rules.get("allowed_distances", {}) or {}
    if isinstance(raw_rules, Sequence) and not isinstance(raw_rules, (str, bytes, Mapping)):
        for rule in raw_rules:
            if isinstance(rule, Mapping):
                yield dict(rule)
        return

    if not isinstance(raw_rules, Mapping):
        return

    for name, raw_rule in raw_rules.items():
        if not isinstance(raw_rule, Mapping):
            continue
        rule = dict(raw_rule)
        if ("from" not in rule or "to" not in rule) and "_to_" in str(name):
            from_site, to_site = str(name).split("_to_", 1)
            rule.setdefault("from", from_site)
            rule.setdefault("to", to_site)
        yield rule


def _repeat_primer_count(selected: Sequence[CandidateMatch]) -> int:
    seen: dict[str, int] = {}
    repeats = 0
    for match in selected:
        if match.role != "primer":
            continue
        key = match.ref_id or match.site
        seen[key] = seen.get(key, 0) + 1
        if seen[key] > 1:
            repeats += 1
    return repeats


def _edit_penalty(record: Any, scoring: Mapping[str, Any]) -> float:
    errors = _field(record, "errors")
    edit_penalties = scoring.get("edit_penalties", {}) or {}

    if errors is not None:
        substitutions, insertions, deletions = errors
        return (
            int(substitutions) * _penalty_value(edit_penalties, "mismatch")
            + int(insertions) * _penalty_value(edit_penalties, "insertion")
            + int(deletions) * _penalty_value(edit_penalties, "deletion")
        )

    penalty = _field(record, "penalty")
    if penalty is not None:
        return abs(float(penalty))

    score = _field(record, "score")
    if score is not None:
        return max(0.0, -float(score))

    return 0.0


def _penalty_value(scoring: Mapping[str, Any], name: str, default: float = 0.0) -> float:
    return abs(float(scoring.get(name, default)))


def _layout_rules(config: Mapping[str, Any]) -> Mapping[str, Any]:
    if "layout_rules" in config:
        return config.get("layout_rules", {}) or {}
    return config


def _layout_scoring(config: Mapping[str, Any]) -> Mapping[str, Any]:
    layout_scoring = dict(config.get("layout_scoring", {}) or {})
    scoring_weights = config.get("scoring_weights", {}) or {}

    if "edit_penalties" not in layout_scoring:
        layout_scoring["edit_penalties"] = {
            "mismatch": abs(float(scoring_weights.get("mismatch", 1))),
            "insertion": abs(float(scoring_weights.get("insertion", 1))),
            "deletion": abs(float(scoring_weights.get("deletion", 1))),
        }

    for key in (
        "wrong_order_penalty",
        "distance_penalty_per_base",
        "extra_flank_penalty",
        "duplicate_site_penalty",
        "repeat_primer_penalty",
    ):
        if key not in layout_scoring and key in scoring_weights:
            layout_scoring[key] = abs(float(scoring_weights[key]))

    return layout_scoring


def _configured_read_length(
    read_length: int | None,
    config: Mapping[str, Any],
    layout_rules: Mapping[str, Any],
) -> int | None:
    configured = read_length or layout_rules.get("read_length") or config.get("read_length")
    if configured is None:
        return None
    return int(configured)


def _load_default_config() -> Mapping[str, Any]:
    config_path = Path(__file__).with_name("config.yaml")
    if not config_path.exists():
        return {}
    try:
        import yaml

        with config_path.open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {}
    except Exception:
        return {}


def _looks_like_config(value: Any) -> bool:
    return isinstance(value, Mapping) and any(
        key in value
        for key in (
            "layout_rules",
            "layout_scoring",
            "scoring_weights",
            "sequences",
            "reference_builder",
        )
    )


def _configured_role(site: str, layout_rules: Mapping[str, Any]) -> str | None:
    metadata = (layout_rules.get("match_sites", {}) or {}).get(site, {}) or {}
    role = metadata.get("role") or metadata.get("kind") or _infer_role(site)
    return str(role) if role else None


def _group_by_site(matches: Iterable[CandidateMatch]) -> dict[str, list[CandidateMatch]]:
    grouped: dict[str, list[CandidateMatch]] = {}
    for match in matches:
        grouped.setdefault(match.site, []).append(match)
    return grouped


def _site_counts(matches: Iterable[CandidateMatch]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for match in matches:
        counts[match.site] = counts.get(match.site, 0) + 1
    return counts


def _unique_matches(matches: Iterable[CandidateMatch]) -> list[CandidateMatch]:
    unique: dict[tuple[Any, ...], CandidateMatch] = {}
    for match in matches:
        unique.setdefault(
            (
                match.site,
                match.start,
                match.end,
                match.ref_id,
                match.element_type,
                id(match.record),
            ),
            match,
        )
    return list(unique.values())


def _candidate_sort_key(candidate: LayoutCandidate) -> tuple[Any, ...]:
    return (
        candidate.total_penalty,
        tuple(_match_sort_key(match) for match in candidate.selected_matches),
    )


def _selected_identity(selected: Sequence[CandidateMatch]) -> tuple[Any, ...]:
    return tuple(
        (
            match.site,
            match.start,
            match.end,
            match.ref_id,
            match.element_type,
            id(match.record),
        )
        for match in selected
    )


def _match_sort_key(match: CandidateMatch) -> tuple[Any, ...]:
    return (match.start, match.end, match.site, match.ref_id or "")


def _candidate_value(match: CandidateMatch, key: str) -> Any:
    if key == "site":
        return match.site
    if key == "role":
        return match.role
    if key == "side":
        return match.side
    if key == "strand":
        return match.strand
    if key == "direction":
        return match.direction
    if key == "ref_id":
        return match.ref_id
    if key == "element_type":
        return match.element_type
    if key == "transform":
        return match.transform
    return _field(match.record, key)


def _value_allowed(actual: Any, expected: Any) -> bool:
    if expected in (None, "*", "any"):
        return True
    if isinstance(expected, Sequence) and not isinstance(expected, (str, bytes, Mapping)):
        return any(_value_allowed(actual, item) for item in expected)
    return _normalize(actual) == _normalize(expected)


def _field(record: Any, name: str) -> Any:
    if record is None:
        return None
    if isinstance(record, Mapping):
        if name in record:
            return record[name]
        metadata = record.get("metadata")
        if isinstance(metadata, Mapping) and name in metadata:
            return metadata[name]
        return None
    if hasattr(record, name):
        return getattr(record, name)
    metadata = getattr(record, "metadata", None)
    if isinstance(metadata, Mapping):
        return metadata.get(name)
    return None


def _metadata_field(record: Any, name: str) -> Any:
    metadata = _field(record, "metadata")
    if isinstance(metadata, Mapping):
        return metadata.get(name)
    return None


def _first_present(values: Sequence[Any]) -> str | None:
    for value in values:
        if value is not None:
            return str(value)
    return None


def _infer_role(value: Any) -> str | None:
    normalized = _normalize(value)
    if "adapter" in normalized:
        return "adapter"
    if "barcode" in normalized:
        return "barcode"
    if "primer" in normalized:
        return "primer"
    return None


def _infer_side_from_site(site: Any, configured_sites: Mapping[str, Any]) -> str | None:
    normalized_site = _normalize(site)
    for configured_site, metadata in configured_sites.items():
        aliases = {configured_site, *(metadata.get("aliases", []) or [])}
        if normalized_site in {_normalize(alias) for alias in aliases}:
            return _optional_str(metadata.get("side"))
    return None


def _normalize(value: Any) -> str:
    return "" if value is None else str(value).strip().lower()


def _optional_str(value: Any) -> str | None:
    if value is None:
        return None
    return str(value)


def _is_record_collection(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray, Mapping))


def _legacy_record(match: CandidateMatch) -> Any:
    if isinstance(match.record, Mapping):
        return match.record
    return {
        "start": match.start,
        "end": match.end,
        "penalty": match.edit_penalty,
        "match_id": match.ref_id or match.site,
        "is_rc": match.strand in ("-", "reverse"),
    }
