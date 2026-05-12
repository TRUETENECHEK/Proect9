#!/usr/bin/env python3
"""mock_data_generator.py

Generate synthetic FASTQ reads for testing.

The script reads adapter sequences and barcode lists from ``config.yaml``
which is expected to be located in the same directory.
It creates 10 000 reads composed of:
    left_adapter + left_barcode + random_amplicon + right_barcode + right_adapter
and optionally introduces sequencing errors (substitutions and indels).

The output is written to ``mock_test_reads.fastq`` in standard FASTQ format.

Dependencies: only the Python standard library and Biopython.
"""

import random
import yaml
from pathlib import Path
from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Configuration handling
# ---------------------------------------------------------------------------
CONFIG_PATH = Path(__file__).with_name("config.yaml")

def load_config(path: Path = CONFIG_PATH) -> dict:
    """Load adapters and barcodes from ``config.yaml``.

    The expected structure (as used by the rest of the repository) is::

        sequences:
          left_fwd:
            adapter: "NNNN..."
            barcodes:
              - "left_bc_1"
              - "left_bc_2"
          right_rev:
            adapter: "NNNN..."
            barcodes:
              - "right_bc_1"
              - "right_bc_2"
    """
    with path.open("r") as fh:
        cfg = yaml.safe_load(fh)
    return cfg.get("sequences", {})

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------
BASES = ["A", "C", "G", "T"]

def random_amplicon(length: int = 100) -> str:
    """Return a random DNA sequence of *length* bases.

    In a real experiment the primer would already be part of the amplicon –
    we therefore treat the full random stretch as the region that follows the
    left barcode and precedes the right barcode.
    """
    return "".join(random.choices(BASES, k=length))

def introduce_errors(seq: str, error_rate: float = 0.05) -> str:
    """Introduce random substitutions and indels into *seq*.

    *error_rate* is the probability that any given position will be altered.
    When an alteration occurs we choose uniformly between a mismatch,
    deletion or insertion (the insertion and deletion share the same probability
    as the mismatch; the exact distribution is not critical for a mock data set).
    """
    i = 0
    new_seq = []
    while i < len(seq):
        if random.random() < error_rate:
            # decide what kind of error to introduce
            choice = random.choice(["mismatch", "deletion", "insertion"])
            if choice == "mismatch":
                original = seq[i]
                alternatives = [b for b in BASES if b != original]
                new_seq.append(random.choice(alternatives))
                i += 1
            elif choice == "deletion":
                # skip the current base – i.e. delete it
                i += 1
            else:  # insertion
                # insert a random base before the current one, keep current base
                new_seq.append(random.choice(BASES))
                # do not advance i – the original base will still be processed
        else:
            new_seq.append(seq[i])
            i += 1
    return "".join(new_seq)

def generate_ideal_read(
    left_adapter: str,
    left_barcode: str,
    right_barcode: str,
    right_adapter: str,
    amplicon_len: int = 100,
) -> str:
    """Construct a perfect read without sequencing errors.

    The order follows the *valid_layouts* definition in ``config.yaml`` –
    left adapter → left barcode → amplicon → right barcode → right adapter.
    """
    amplicon = random_amplicon(amplicon_len)
    return f"{left_adapter}{left_barcode}{amplicon}{right_barcode}{right_adapter}"

# ---------------------------------------------------------------------------
# Main generation routine
# ---------------------------------------------------------------------------
def main(
    n_reads: int = 10_000,
    error_rate: float = 0.05,
    output_path: Path = Path(__file__).with_name("mock_test_reads.fastq"),
):
    cfg = load_config()
    left_cfg = cfg.get("left_fwd", {})
    right_cfg = cfg.get("right_rev", {})

    left_adapter = left_cfg.get("adapter", "")
    right_adapter = right_cfg.get("adapter", "")
    left_barcodes: List[str] = left_cfg.get("barcodes", [])
    right_barcodes: List[str] = right_cfg.get("barcodes", [])

    if not (left_adapter and right_adapter and left_barcodes and right_barcodes):
        raise ValueError("Config file does not contain required adapters / barcodes.")

    records = []
    for i in range(1, n_reads + 1):
        left_bc = random.choice(left_barcodes)
        right_bc = random.choice(right_barcodes)
        seq = generate_ideal_read(
            left_adapter=left_adapter,
            left_barcode=left_bc,
            right_barcode=right_bc,
            right_adapter=right_adapter,
        )
        seq = introduce_errors(seq, error_rate=error_rate)
        # Phred quality – we use a constant high quality (40) which maps to 'I'
        qual_char = chr(33 + 40)  # 'I'
        qual = qual_char * len(seq)
        record = SeqRecord(Seq(seq), id=f"read_{i}", description="", letter_annotations={"phred_quality": [40] * len(seq)})
        # Biopython will translate the integer list to the appropriate ASCII chars.
        records.append(record)

    with output_path.open("w") as out_handle:
        SeqIO.write(records, out_handle, "fastq")
    print(f"Generated {n_reads} synthetic reads → {output_path}")

if __name__ == "__main__":
    main()
