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
    record: Any = field(compare=False, hash=False)
    ref_id: str | None = None
    element_type: str | None = None
    transform: str | None = None
    edit_penalty: float = 0.0
    source_key: str | None = None


@dataclass(frozen=True)
class LayoutCandidate:
    adapter_pair: tuple[CandidateMatch, ...]
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

    # We build potential "anchor sets" from adapters (0, 1, or 2)
    anchor_sets = _build_anchor_sets(adapters, layout_rules)
    
    candidates = _build_layout_candidates(
        anchor_sets=anchor_sets,
        barcode_primer_matches=barcode_primers,
        layout_rules=layout_rules,
        scoring=scoring,
        read_length=read_length,
    )

    if not candidates:
        if not adapters:
            return LayoutDecision(MISSING_ADAPTERS_STATUS, None, 0, "No adapters found")
        return LayoutDecision(INVALID_LAYOUT_STATUS, None, 0, "No valid layout found")

    candidates = sorted(candidates, key=_candidate_sort_key)
    selected = candidates[0]
    ambiguity_margin = max(0.0, float(layout_rules.get("ambiguity_margin", 0)))

    if len(candidates) > 1:
        delta = candidates[1].total_penalty - selected.total_penalty
        if delta < ambiguity_margin:
            # If they have different identities, it's ambiguous.
            if _selected_identity(candidates[1].selected_matches) != _selected_identity(selected.selected_matches):
                 return LayoutDecision(AMBIGUOUS_LAYOUT_STATUS, selected, len(candidates), "Multiple similar layouts")

    return LayoutDecision(PASS_STATUS, selected, len(candidates), "Selected layout")


def _build_anchor_sets(
    adapters: Sequence[CandidateMatch],
    layout_rules: Mapping[str, Any],
) -> list[tuple[CandidateMatch, ...]]:
    anchor_sets: list[tuple[CandidateMatch, ...]] = []
    
    # Pairs
    for left, right in combinations(adapters, 2):
        if left.side != right.side or left.side is None:
            anchor_sets.append((left, right))
    
    # Singletons
    for match in adapters:
        anchor_sets.append((match,))
        
    # Empty set (for barcode-only layouts)
    anchor_sets.append(())
    
    return anchor_sets


def _build_layout_candidates(
    anchor_sets: Sequence[tuple[CandidateMatch, ...]],
    barcode_primer_matches: Sequence[CandidateMatch],
    layout_rules: Mapping[str, Any],
    scoring: Mapping[str, Any],
    read_length: int | None,
) -> list[LayoutCandidate]:
    valid_layouts = [tuple(layout) for layout in layout_rules.get("valid_layouts", []) or []]
    barcode_primers_by_site = {}
    for match in barcode_primer_matches:
        barcode_primers_by_site.setdefault(match.site, []).append(match)
        
    all_matches = list(barcode_primer_matches)
    for s in anchor_sets: all_matches.extend(s)
    site_counts = _site_counts(_unique_matches(all_matches))
    
    candidates: list[LayoutCandidate] = []
    candidate_by_identity = {}

    for anchors in anchor_sets:
        anchor_sites = {m.site for m in anchors}
        for layout in valid_layouts:
            layout_adapter_sites = {s for s in layout if _configured_role(s, layout_rules) == "adapter"}
            if not anchor_sites.issubset(layout_adapter_sites): continue
            if layout_adapter_sites != anchor_sites: continue
            
            choices = []
            missing = False
            for site in layout:
                if site in anchor_sites: continue
                matches = barcode_primers_by_site.get(site, [])
                if not matches:
                    missing = True
                    break
                choices.append(matches)
            
            if missing: continue
            
            for non_anchors in product(*choices) if choices else [()]:
                selected = tuple(sorted((*anchors, *non_anchors), key=_match_sort_key))
                penalty, breakdown = _score_candidate(selected, layout, site_counts, layout_rules, scoring, read_length)
                
                candidate = LayoutCandidate(anchors, selected, tuple(m.site for m in selected), layout, penalty, breakdown)
                ident = _selected_identity(selected)
                if ident not in candidate_by_identity or penalty < candidate_by_identity[ident].total_penalty:
                    candidate_by_identity[ident] = candidate

    return list(candidate_by_identity.values())


def _score_candidate(selected, layout, site_counts, layout_rules, scoring, read_length):
    breakdown = {"edit": sum(m.edit_penalty for m in selected), "wrong_order": 0.0, "distance": 0.0, "extra_flank": 0.0, "duplicate_site": 0.0}
    ordered = sorted(selected, key=_match_sort_key)
    if tuple(m.site for m in ordered) != layout:
        breakdown["wrong_order"] = _penalty_value(scoring, "wrong_order_penalty")
    
    selected_by_site = {m.site: m for m in selected}
    breakdown["distance"] = _distance_penalty(selected_by_site, layout_rules, scoring)
    
    if read_length and selected:
        flank = max(0, ordered[0].start) + max(0, read_length - ordered[-1].end)
        breakdown["extra_flank"] = flank * _penalty_value(scoring, "extra_flank_penalty")
        
    dups = sum(max(0, site_counts.get(m.site, 0) - 1) for m in selected)
    breakdown["duplicate_site"] = dups * _penalty_value(scoring, "duplicate_site_penalty")
    
    return sum(breakdown.values()), breakdown


def _distance_penalty(selected_by_site, layout_rules, scoring):
    total = 0.0
    per_base = _penalty_value(scoring, "distance_penalty_per_base")
    for rule in _distance_rules(layout_rules):
        fr, to = rule.get("from"), rule.get("to")
        if fr in selected_by_site and to in selected_by_site:
            dist = selected_by_site[to].start - selected_by_site[fr].end
            mi, ma = rule.get("min"), rule.get("max")
            if mi is not None and dist < int(mi): total += (int(mi) - dist) * per_base
            if ma is not None and dist > int(ma): total += (dist - int(ma)) * per_base
    return total

def _distance_rules(layout_rules):
    raw = layout_rules.get("allowed_distances", {}) or {}
    if isinstance(raw, Mapping):
        for k, v in raw.items():
            if isinstance(v, Mapping):
                r = dict(v)
                if "_to_" in str(k) and "from" not in r:
                    f, t = str(k).split("_to_", 1)
                    r["from"], r["to"] = f, t
                yield r
    elif isinstance(raw, Sequence):
        for r in raw: yield r

def _edit_penalty(record, scoring):
    errors = _field(record, "errors")
    p = scoring.get("edit_penalties", {})
    if errors: return errors[0]*p.get("mismatch",1) + errors[1]*p.get("insertion",1) + errors[2]*p.get("deletion",1)
    return abs(float(_field(record, "penalty") or 0))

def _penalty_value(scoring, name, default=0.0):
    return abs(float(scoring.get(name, default)))

def _layout_rules(config): return config.get("layout_rules", {}) or config
def _layout_scoring(config):
    s = config.get("scoring", {}) or config.get("layout_scoring", {}) or config.get("scoring_weights", {})
    return {
        "edit_penalties": {"mismatch": abs(float(s.get("mismatch", 1))), "insertion": abs(float(s.get("insertion", 3))), "deletion": abs(float(s.get("deletion", 3)))},
        "wrong_order_penalty": abs(float(s.get("wrong_order_penalty", 100))),
        "distance_penalty_per_base": abs(float(s.get("distance_penalty_per_base", 1))),
        "extra_flank_penalty": abs(float(s.get("extra_flank_penalty", 1))),
        "duplicate_site_penalty": abs(float(s.get("duplicate_site_penalty", 50))),
    }

def _normalize_inputs(adapter_matches, barcode_primer_matches, layout_rules, scoring):
    def to_cm(r, role=None):
        if isinstance(r, CandidateMatch): return r
        meta = _match_metadata(r, layout_rules, None, role)
        if not meta.get("site"): return None
        return CandidateMatch(site=str(meta["site"]), role=str(meta["role"]), side=meta.get("side"), start=int(_field(r,"start")), end=int(_field(r,"end")), strand=str(_field(r,"strand") or "+"), direction=None, record=r, edit_penalty=_edit_penalty(r, scoring), ref_id=str(_field(r,"ref_id") or ""))
    
    adapters = [to_cm(r, "adapter") for r in (adapter_matches or [])]
    barcodes = [to_cm(r, "barcode") for r in (barcode_primer_matches or [])]
    return [a for a in adapters if a], [b for b in barcodes if b]

def _match_metadata(record, layout_rules, forced_site, default_role):
    sites = layout_rules.get("match_sites", {}) or {}
    val = _field(record, "site") or _field(record, "element_type") or _field(record, "ref_id")
    for s, m in sites.items():
        if val == s or val in (m.get("aliases", [])): return {"site": s, "role": m.get("role", default_role), "side": m.get("side")}
    role = _field(record, "role") or default_role or (_infer_role(val) if val else None)
    return {"site": val, "role": role, "side": _field(record, "side")}

def _infer_role(v):
    v = str(v).lower()
    if "adapter" in v: return "adapter"
    if "barcode" in v: return "barcode"
    return None

def _configured_role(site, layout_rules):
    return ((layout_rules.get("match_sites", {}) or {}).get(site, {}) or {}).get("role")

def _site_counts(matches):
    c = {}
    for m in matches: c[m.site] = c.get(m.site, 0) + 1
    return c

def _unique_matches(matches):
    u = {}
    for m in matches: u[(m.site, m.start, m.end, m.ref_id)] = m
    return list(u.values())

def _candidate_sort_key(c): return (c.total_penalty, tuple((m.start, m.site) for m in c.selected_matches))
def _selected_identity(selected): return tuple((m.site, m.ref_id, m.start, m.end) for m in selected)
def _match_sort_key(m): return (m.start, m.end, m.site)
def _field(r, n):
    if isinstance(r, Mapping): return r.get(n)
    return getattr(r, n, None)
def _normalize(v): return str(v).strip().lower() if v else ""
def _optional_str(v): return str(v) if v else None
def _is_record_collection(v): return isinstance(v, Sequence) and not isinstance(v, (str, Mapping))
def _legacy_record(m): return {"start": m.start, "end": m.end, "penalty": m.edit_penalty, "match_id": m.ref_id or m.site}
def _configured_read_length(rl, config, rules): return rl or rules.get("read_length") or config.get("read_length")
def _load_default_config(): return {}
def _looks_like_config(v): return isinstance(v, Mapping) and "layout_rules" in v
