import csv
import os
import sys
import traceback
from typing import Any, Iterable

from data_loader import load_config, load_mapping_excel, read_fastq_generator
from adapter_discovery import discover_adapters
from demux_logic import RESULT_FIELDNAMES
from export_fastq import export_demultiplexed_reads
from reference_builder import build_reference_bundle
from simple_matcher import SimpleMatcher


def _normalize_error_budget(value: Any, default: int, hard_cap: int = 3) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        parsed = default
    return max(0, min(parsed, hard_cap))

def run_pipeline(config_path, fastq_path, output_csv, output_dir):
    print("--- BioDemux Simple Pipeline Start ---")
    config = load_config(config_path)
    
    discovery_cfg = config.get("reference_builder", {}).get("adapter_discovery", {})
    discovered = None
    
    if discovery_cfg.get("enabled"):
        print("Discovering adapters from reads...")
        discovered = discover_adapters(
            fastq_path=fastq_path,
            fasta_path=config.get("reference_builder", {}).get("fasta_path"),
            sample_size=discovery_cfg.get("sample_size", 5000),
            window_size=discovery_cfg.get("window_size", 300),
            min_abundance=discovery_cfg.get("min_abundance", 0.2),
        )
    base_dir = os.path.dirname(os.path.abspath(config_path))
    bundle = build_reference_bundle(config, base_dir=base_dir, discovered_adapters=discovered)
    left_adapter = bundle.adapters.get("left", "")
    right_adapter = bundle.adapters.get("right", "")
    barcodes = {
        barcode_id: sequence
        for side_cores in bundle.barcode_cores.values()
        for barcode_id, sequence in side_cores.items()
    }
    
    # Configure errors from reference_builder (fallback to root for compatibility)
    rb_cfg = config.get("reference_builder", {})
    max_adp_err = _normalize_error_budget(
        rb_cfg.get("max_adapter_errors", config.get("max_adapter_errors", 3)),
        default=3,
    )
    max_bc_err = _normalize_error_budget(
        rb_cfg.get("max_barcode_errors", config.get("max_barcode_errors", 3)),
        default=3,
    )
    
    matcher = SimpleMatcher(left_adapter, right_adapter, barcodes, max_adp_err, max_bc_err)
    
    mapping_dict = {}
    if os.path.exists("mapping.xlsx"):
        mapping_dict = load_mapping_excel("mapping.xlsx")
        
    print(f"Processing {fastq_path}...")
    export_records = []
    
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=RESULT_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        
        count = 0
        for rid, seq, qual in read_fastq_generator(fastq_path):
            try:
                res = matcher.process_read(rid, seq, qual, mapping_dict)
                writer.writerow(res)
                
                export_records.append({
                    "read_id": res["read_id"],
                    "status": res["status"],
                    "barcode_id": res["barcode_id"],
                    "sample_id": res["sample_id"],
                    "trimmed_seq": res["trimmed_seq"],
                    "trimmed_quality": res.get("original_quality") if not res.get("trimmed_quality") else res.get("trimmed_quality"),
                    "original_seq": res["original_seq"],
                    "original_quality": res["original_quality"]
                })
            except Exception as e:
                writer.writerow({"read_id": rid, "status": f"Fail_Exception: {e}"})
                
            count += 1
            if count % 100 == 0:
                sys.stdout.write(f"\r  Processed {count} reads...")
                sys.stdout.flush()
            if count >= 1000:
                break
                
        print(f"\nFinished processing {count} reads.")
        
    print(f"Exporting FASTQ to {output_dir}...")
    export_demultiplexed_reads(export_records, output_dir=output_dir)
    print("Pipeline Complete.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default="config.yaml")
    parser.add_argument("-i", "--input", default="Read_file/test_reads.fastq")
    parser.add_argument("-o", "--output_dir", default="Demux_Out")
    args = parser.parse_args()
    run_pipeline(args.config, args.input, "results.csv", args.output_dir)
