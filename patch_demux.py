import re

with open('demux_logic.py', 'r') as f:
    content = f.read()

new_func = """def process_read(read_id, read_seq, adapters, barcodes, config, mapping_dict=None, max_errors=2, quality=None):
    weights_dict = config.get("scoring_weights", {})
    layout_rules = config.get("layout_rules", {})
    
    candidates_dict = {
        "left_adapter": None,
        "left_barcode": None,
        "right_adapter": None,
        "right_barcode": None
    }
    
    left_fwd = adapters.get("left_fwd", {})
    right_rev = adapters.get("right_rev", {})
    
    la_seq = left_fwd.get("adapter")
    if la_seq and la_seq != "NNNNNNNNNN":
        match_la = find_best_match(read_seq, la_seq, max_errors)
        if match_la:
            start, end, errors, is_rc = match_la
            score = calculate_element_score(errors, weights_dict)
            candidates_dict["left_adapter"] = {
                "start": start, "end": end, "score": score, 
                "match_id": "left_adapter", "is_rc": is_rc
            }
            
    ra_seq = right_rev.get("adapter")
    if ra_seq and ra_seq != "NNNNNNNNNN":
        match_ra = find_best_match(read_seq, ra_seq, max_errors)
        if match_ra:
            start, end, errors, is_rc = match_ra
            score = calculate_element_score(errors, weights_dict)
            candidates_dict["right_adapter"] = {
                "start": start, "end": end, "score": score, 
                "match_id": "right_adapter", "is_rc": is_rc
            }

    left_bc_ids = left_fwd.get("barcodes", [])
    best_lb_id, best_lb_match, best_lb_score = None, None, -float('inf')
    for b_id, b_seq in barcodes.items():
        if b_id in left_bc_ids:
            match_b = find_best_match(read_seq, b_seq, max_errors=1)
            if match_b:
                start, end, errors, is_rc = match_b
                score = calculate_element_score(errors, weights_dict)
                if score > best_lb_score:
                    best_lb_score = score
                    best_lb_id = b_id
                    best_lb_match = (start, end, score, is_rc)
                    
    if best_lb_match:
        candidates_dict["left_barcode"] = {
            "start": best_lb_match[0], "end": best_lb_match[1], 
            "score": best_lb_match[2], "match_id": best_lb_id, "is_rc": best_lb_match[3]
        }
        
    right_bc_ids = right_rev.get("barcodes", [])
    best_rb_id, best_rb_match, best_rb_score = None, None, -float('inf')
    for b_id, b_seq in barcodes.items():
        if b_id in right_bc_ids:
            match_b = find_best_match(read_seq, b_seq, max_errors=1)
            if match_b:
                start, end, errors, is_rc = match_b
                score = calculate_element_score(errors, weights_dict)
                if score > best_rb_score:
                    best_rb_score = score
                    best_rb_id = b_id
                    best_rb_match = (start, end, score, is_rc)
                    
    if best_rb_match:
        candidates_dict["right_barcode"] = {
            "start": best_rb_match[0], "end": best_rb_match[1], 
            "score": best_rb_match[2], "match_id": best_rb_id, "is_rc": best_rb_match[3]
        }

    filter_status, filtered_elements = filter_valid_pairs(candidates_dict)
    
    if filter_status == "Discarded":
        return {
            "read_id": read_id,
            "barcode_id": "Unknown",
            "sample_id": None,
            "score": 0,
            "status": "Discarded",
            "trimmed_seq": "",
            "trimmed_quality": None,
            "original_seq": read_seq,
            "original_quality": quality
        }
        
    total_score, status = evaluate_read_layout(filtered_elements, layout_rules, weights_dict)
    
    trimmed_seq = ""
    trimmed_quality = None
    lb = filtered_elements.get("left_barcode")
    rb = filtered_elements.get("right_barcode")
    
    if lb and rb:
        if lb["start"] < rb["start"]:
            start_pos = lb["end"]
            end_pos = rb["start"]
        else:
            start_pos = rb["end"]
            end_pos = lb["start"]
            
        if start_pos < end_pos:
            trimmed_seq = read_seq[start_pos:end_pos]
            if quality:
                trimmed_quality = quality[start_pos:end_pos]

    barcode_id = f"{lb['match_id']}_{rb['match_id']}" if (lb and rb) else "Unknown"
    
    sample_id = None
    if mapping_dict is not None:
        if barcode_id in mapping_dict:
            sample_id = mapping_dict[barcode_id]
        else:
            return {
                "read_id": read_id,
                "barcode_id": barcode_id,
                "sample_id": None,
                "score": total_score,
                "status": "Discarded",
                "trimmed_seq": "",
                "trimmed_quality": None,
                "original_seq": read_seq,
                "original_quality": quality
            }
        
    return {
        "read_id": read_id,
        "barcode_id": barcode_id,
        "sample_id": sample_id,
        "score": total_score,
        "status": status,
        "trimmed_seq": trimmed_seq,
        "trimmed_quality": trimmed_quality,
        "original_seq": read_seq,
        "original_quality": quality
    }"""

old_func_pattern = re.compile(r'def process_read\(read_id, read_seq, adapters, barcodes, config, max_errors=2\):.*?return \{(.*?)\}', re.DOTALL)

new_content = old_func_pattern.sub(new_func, content)

# Also update save_process_results to include sample_id
new_content = new_content.replace(
    'fieldnames=["read_id", "barcode_id", "score", "status", "trimmed_seq"]',
    'fieldnames=["read_id", "barcode_id", "sample_id", "score", "status", "trimmed_seq"]'
)

with open('demux_logic.py', 'w') as f:
    f.write(new_content)
