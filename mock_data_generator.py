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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Mapping

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from reference_builder import ReferenceBundle, ReferenceVariant
from sequence_matcher import complement, reverse_complement

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


MOCK_LEFT_ADAPTER = "AACCGGTTAGCACCTGAACT"
MOCK_RIGHT_ADAPTER = "TTGGAACCTAGTGGACTTCA"
MOCK_LEFT_BARCODES = {
    "LB1": "CGTACGATCGTA",
    "LB2": "CGTACGATCGTG",
}
MOCK_RIGHT_BARCODES = {
    "RB1": "GATCTAGCTAGC",
    "RB2": "CTAGGATCCGTA",
}
MOCK_AMPLICON = "ATATCCGGATTAAGGCCATATGCGCGTTAACCGATAT"
MOCK_PRIMER = "CCGGATTAA"
MOCK_AMBIGUOUS_LEFT_BARCODE = "CGTACGATCGTT"


@dataclass(frozen=True)
class MockReadCase:
    """A compact process_read smoke-test case."""

    read_id: str
    sequence: str
    expected_status: str
    max_errors: int = 2
    expected_trimmed_seq: str | None = None
    expected_left_barcode_id: str | None = None
    expected_right_barcode_id: str | None = None
    min_penalties: Mapping[str, float] = field(default_factory=dict)


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


def introduce_deterministic_error(
    seq: str,
    error_type: str,
    index: int,
    inserted_base: str = "A",
) -> str:
    """Return *seq* with exactly one mismatch, insertion or deletion."""
    error_type = error_type.lower()
    if error_type not in {"mismatch", "insertion", "deletion"}:
        raise ValueError("error_type must be one of: mismatch, insertion, deletion")

    if error_type == "insertion":
        if not 0 <= index <= len(seq):
            raise IndexError("insertion index is outside the sequence")
        return f"{seq[:index]}{inserted_base}{seq[index:]}"

    if not 0 <= index < len(seq):
        raise IndexError("edit index is outside the sequence")

    if error_type == "deletion":
        return f"{seq[:index]}{seq[index + 1:]}"

    replacement = next(base for base in BASES if base != seq[index].upper())
    return f"{seq[:index]}{replacement}{seq[index + 1:]}"


def build_mock_config() -> dict[str, Any]:
    """Return a small architecture config for fast process_read smoke tests."""
    return {
        "barcode_max_errors": 1,
        "barcode_window": 2,
        "barcode_window_slack": 0,
        "layout_scoring": {
            "edit_penalties": {
                "mismatch": 2,
                "insertion": 3,
                "deletion": 3,
            },
            "wrong_order_penalty": 100,
            "distance_penalty_per_base": 1,
            "extra_flank_penalty": 1,
            "duplicate_site_penalty": 50,
            "repeat_primer_penalty": 50,
        },
        "layout_rules": {
            "ambiguity_margin": 0,
            "match_sites": {
                "left_adapter": {"role": "adapter", "side": "left"},
                "right_adapter": {"role": "adapter", "side": "right"},
                "left_barcode": {"role": "barcode", "side": "left"},
                "right_barcode": {"role": "barcode", "side": "right"},
            },
            "allowed_adapter_pairs": [
                {
                    "first": {"side": "left", "strand": "any", "direction": "any"},
                    "second": {"side": "right", "strand": "any", "direction": "any"},
                },
                {
                    "first": {"side": "right", "strand": "any", "direction": "any"},
                    "second": {"side": "left", "strand": "any", "direction": "any"},
                },
            ],
            "valid_layouts": [
                ["left_adapter", "left_barcode", "right_barcode", "right_adapter"],
                ["right_adapter", "right_barcode", "left_barcode", "left_adapter"],
            ],
            "allowed_distances": {},
        },
        "primer_qc": {
            "enabled": True,
            "sequences": {"mock_primer": MOCK_PRIMER},
            "max_errors": 0,
            "required": "any",
            "warn_on_repeat": True,
        },
    }


def build_mock_reference_bundle() -> ReferenceBundle:
    """Build an in-memory ReferenceBundle with normal and reverse-complement variants."""
    variants: list[ReferenceVariant] = []
    variants.extend(_build_mock_variants(MOCK_LEFT_ADAPTER, "adapter", "left", "left_adapter"))
    variants.extend(_build_mock_variants(MOCK_RIGHT_ADAPTER, "adapter", "right", "right_adapter"))

    for barcode_id, barcode in MOCK_LEFT_BARCODES.items():
        variants.extend(_build_mock_variants(barcode, "barcode_core", "left", barcode_id))
    for barcode_id, barcode in MOCK_RIGHT_BARCODES.items():
        variants.extend(_build_mock_variants(barcode, "barcode_core", "right", barcode_id))

    return ReferenceBundle(
        source_path="in_memory_mock",
        edge="explicit",
        edge_by_side={"left": "explicit", "right": "explicit"},
        adapters={"left": MOCK_LEFT_ADAPTER, "right": MOCK_RIGHT_ADAPTER},
        barcode_cores={
            "left": dict(MOCK_LEFT_BARCODES),
            "right": dict(MOCK_RIGHT_BARCODES),
        },
        variants=variants,
    )


def build_mock_read(
    left_barcode_id: str = "LB1",
    right_barcode_id: str = "RB1",
    amplicon: str = MOCK_AMPLICON,
    orientation: str = "forward",
) -> str:
    """Build a compact read in forward or reverse-complement orientation."""
    read = (
        f"{MOCK_LEFT_ADAPTER}"
        f"{MOCK_LEFT_BARCODES[left_barcode_id]}"
        f"{amplicon}"
        f"{MOCK_RIGHT_BARCODES[right_barcode_id]}"
        f"{MOCK_RIGHT_ADAPTER}"
    )
    if orientation in {"forward", "normal"}:
        return read
    if orientation in {"reverse_complement", "rc"}:
        return reverse_complement(read)
    raise ValueError("orientation must be 'forward' or 'reverse_complement'")


def generate_pipeline_smoke_cases() -> list[MockReadCase]:
    """Generate tiny process_read cases covering the key architecture branches."""
    normal_read = build_mock_read()
    repeated_primer_amplicon = f"{MOCK_PRIMER}TT{MOCK_PRIMER}AA"

    return [
        MockReadCase(
            read_id="normal_forward",
            sequence=normal_read,
            expected_status="Pass",
            expected_trimmed_seq=MOCK_AMPLICON,
            expected_left_barcode_id="LB1",
            expected_right_barcode_id="RB1",
        ),
        MockReadCase(
            read_id="reverse_complement",
            sequence=build_mock_read(orientation="reverse_complement"),
            expected_status="Pass",
            expected_trimmed_seq=reverse_complement(MOCK_AMPLICON),
            expected_left_barcode_id="LB1",
            expected_right_barcode_id="RB1",
        ),
        MockReadCase(
            read_id="adapter_mismatch",
            sequence=introduce_deterministic_error(normal_read, "mismatch", 2),
            expected_status="Pass",
            expected_trimmed_seq=MOCK_AMPLICON,
            expected_left_barcode_id="LB1",
            expected_right_barcode_id="RB1",
            min_penalties={"adapter_penalty": 2},
        ),
        MockReadCase(
            read_id="adapter_insertion",
            sequence=introduce_deterministic_error(normal_read, "insertion", 3),
            expected_status="Pass",
            expected_trimmed_seq=MOCK_AMPLICON,
            expected_left_barcode_id="LB1",
            expected_right_barcode_id="RB1",
            min_penalties={"adapter_penalty": 3},
        ),
        MockReadCase(
            read_id="adapter_deletion",
            sequence=introduce_deterministic_error(normal_read, "deletion", 4),
            expected_status="Pass",
            expected_trimmed_seq=MOCK_AMPLICON,
            expected_left_barcode_id="LB1",
            expected_right_barcode_id="RB1",
            min_penalties={"adapter_penalty": 3},
        ),
        MockReadCase(
            read_id="missing_right_adapter",
            sequence=(
                f"{MOCK_LEFT_ADAPTER}{MOCK_LEFT_BARCODES['LB1']}"
                f"{MOCK_AMPLICON}{MOCK_RIGHT_BARCODES['RB1']}"
            ),
            expected_status="MissingAdapters",
            max_errors=0,
        ),
        MockReadCase(
            read_id="two_left_adapters",
            sequence=(
                f"{MOCK_LEFT_ADAPTER}{MOCK_LEFT_BARCODES['LB1']}"
                f"{MOCK_AMPLICON}{MOCK_LEFT_BARCODES['LB1']}{MOCK_LEFT_ADAPTER}"
            ),
            expected_status="InvalidLayout",
        ),
        MockReadCase(
            read_id="repeat_primer_inside_amplicon",
            sequence=build_mock_read(amplicon=repeated_primer_amplicon),
            expected_status="PrimerWarning",
            expected_trimmed_seq=repeated_primer_amplicon,
            expected_left_barcode_id="LB1",
            expected_right_barcode_id="RB1",
            min_penalties={"primer_penalty": 50},
        ),
        MockReadCase(
            read_id="ambiguous_left_barcode",
            sequence=(
                f"{MOCK_LEFT_ADAPTER}{MOCK_AMBIGUOUS_LEFT_BARCODE}"
                f"{MOCK_AMPLICON}{MOCK_RIGHT_BARCODES['RB1']}{MOCK_RIGHT_ADAPTER}"
            ),
            expected_status="AmbiguousBarcode",
            min_penalties={"barcode_penalty": 2},
        ),
        MockReadCase(
            read_id="leading_tail_before_first_adapter",
            sequence=f"GGGG{normal_read}",
            expected_status="Pass",
            expected_trimmed_seq=MOCK_AMPLICON,
            expected_left_barcode_id="LB1",
            expected_right_barcode_id="RB1",
            min_penalties={"layout_penalty": 4},
        ),
        MockReadCase(
            read_id="right_adapter_at_start_bad_context",
            sequence=(
                f"{MOCK_RIGHT_ADAPTER}{MOCK_LEFT_ADAPTER}"
                f"{MOCK_LEFT_BARCODES['LB1']}{MOCK_AMPLICON}{MOCK_RIGHT_BARCODES['RB1']}"
            ),
            expected_status="MissingBarcode",
        ),
    ]


def _build_mock_variants(
    sequence: str,
    reference_type: str,
    side: str,
    source_id: str,
) -> list[ReferenceVariant]:
    return [
        _make_mock_variant(sequence, reference_type, side, source_id, transform)
        for transform in ("forward", "reverse", "reverse_complement", "complement")
    ]


def _make_mock_variant(
    sequence: str,
    reference_type: str,
    side: str,
    source_id: str,
    transform: str,
) -> ReferenceVariant:
    transformed, strand, sign, direction = _transform_mock_sequence(sequence, transform)
    return ReferenceVariant(
        sequence=transformed,
        reference_type=reference_type,
        side=side,
        source_id=source_id,
        transform=transform,
        strand=strand,
        sign=sign,
        direction=direction,
    )


def _transform_mock_sequence(sequence: str, transform: str) -> tuple[str, str, str, str]:
    if transform == "forward":
        return sequence, "forward", "+", "5_to_3"
    if transform == "reverse":
        return sequence[::-1].upper(), "forward", "+", "3_to_5"
    if transform == "reverse_complement":
        return reverse_complement(sequence), "reverse", "-", "5_to_3"
    if transform == "complement":
        return complement(sequence), "reverse", "-", "3_to_5"
    raise ValueError(f"Unsupported transform: {transform}")

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
