from __future__ import annotations
import csv
import os
from typing import Any, Mapping, Sequence
import regex
from data_loader import load_fasta_dict
from matcher_core import (
    PASS_STATUS,
    LayoutCandidate,
    LayoutDecision,
    filter_valid_pairs,
)
from reference_builder import ReferenceBundle, ReferenceVariant
from sequence_matcher import reverse_complement

RESULT_FIELDNAMES = [
    "read_id", "barcode_id", "sample_id", "score", "status", "trimmed_seq",
    "left_barcode_id", "right_barcode_id", "adapter_penalty", "barcode_penalty",
    "primer_penalty", "layout_penalty", "total_penalty", "trim_start", "trim_end"
]

def process_read(read_id, read_seq, adapters=None, barcodes=None, config=None, mapping_dict=None, max_errors=2, quality=None, reference_bundle=None):
    config = config or {}
    bundle = reference_bundle or adapters
    if not isinstance(bundle, ReferenceBundle): return _base_result(read_id, read_seq, quality, "Fail_NoBundle")
    
    adapter_variants = [v for v in bundle.variants if v.reference_type == "adapter"]
    barcode_variants = [v for v in bundle.variants if v.reference_type == "barcode_core"]
    
    all_adapter_matches = _find_matches(read_seq, adapter_variants, config.get("reference_builder",{}).get("max_adapter_errors", 3))
    all_barcode_matches = _find_matches(read_seq, barcode_variants, config.get("reference_builder",{}).get("max_barcode_errors", 3))
    
    decision = filter_valid_pairs(all_adapter_matches, all_barcode_matches, config, len(read_seq))
    
    if decision.status != PASS_STATUS:
        return _finish_err(read_id, read_seq, quality, decision)
    
    cand = decision.selected_candidate
    res = _base_result(read_id, read_seq, quality, PASS_STATUS)
    
    left_bc = next((m for m in cand.selected_matches if m.site == "left_barcode"), None)
    right_bc = next((m for m in cand.selected_matches if m.site == "right_barcode"), None)
    
    l_id = left_bc.ref_id if left_bc else "None"
    r_id = right_bc.ref_id if right_bc else "None"
    
    res["left_barcode_id"] = l_id
    res["right_barcode_id"] = r_id
    res["barcode_id"] = f"{l_id}_{r_id}"
    
    res["total_penalty"] = cand.total_penalty
    res["adapter_penalty"] = sum(m.edit_penalty for m in cand.selected_matches if m.role == "adapter")
    res["barcode_penalty"] = sum(m.edit_penalty for m in cand.selected_matches if m.role == "barcode")
    res["layout_penalty"] = cand.total_penalty - res["adapter_penalty"] - res["barcode_penalty"]
    res["score"] = -cand.total_penalty
    
    bc_matches = sorted([m for m in cand.selected_matches if m.role == "barcode"], key=lambda x: x.start)
    if len(bc_matches) >= 2:
        ts, te = bc_matches[0].end, bc_matches[-1].start
        if ts < te:
            res["trim_start"], res["trim_end"] = ts, te
            res["trimmed_seq"] = read_seq[ts:te]
            if quality: res["trimmed_quality"] = quality[ts:te]
    
    if mapping_dict:
        res["sample_id"] = mapping_dict.get(res["barcode_id"])
    
    return res

def _find_matches(seq, variants, max_err):
    matches = []
    for v in variants:
        pat = regex.compile(f"(?e)(?:{regex.escape(v.sequence)}){{e<={max_err}}}")
        for m in pat.finditer(seq, overlapped=True):
            matches.append({
                "site": f"{v.side}_{'barcode' if v.reference_type=='barcode_core' else v.reference_type}",
                "role": 'barcode' if v.reference_type=='barcode_core' else v.reference_type,
                "side": v.side,
                "start": m.start(),
                "end": m.end(),
                "errors": m.fuzzy_counts,
                "ref_id": v.source_id,
                "transform": v.transform
            })
    return matches

def _base_result(rid, seq, qual, status):
    return {
        "read_id": rid, "barcode_id": "Unknown", "sample_id": None, "score": 0.0,
        "status": status, "trimmed_seq": "", "left_barcode_id": None, "right_barcode_id": None,
        "adapter_penalty": 0.0, "barcode_penalty": 0.0, "primer_penalty": 0.0,
        "layout_penalty": 0.0, "total_penalty": 0.0, "trim_start": None, "trim_end": None,
        "original_seq": seq, "original_quality": qual
    }

def _finish_err(rid, seq, qual, decision):
    res = _base_result(rid, seq, qual, decision.status)
    res["diagnostic_message"] = decision.reason
    if decision.selected_candidate:
        cand = decision.selected_candidate
        res["total_penalty"] = cand.total_penalty
        res["adapter_penalty"] = sum(m.edit_penalty for m in cand.selected_matches if m.role == "adapter")
        res["barcode_penalty"] = sum(m.edit_penalty for m in cand.selected_matches if m.role == "barcode")
        res["layout_penalty"] = cand.total_penalty - res["adapter_penalty"] - res["barcode_penalty"]
    return res
