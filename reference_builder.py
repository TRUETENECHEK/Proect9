import os
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Tuple

from data_loader import load_config, load_fasta_dict


TRANSFORMS = ("forward", "reverse", "reverse_complement", "complement")
VALID_EDGES = ("prefix", "suffix", "auto")
DNA_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


@dataclass(frozen=True)
class ReferenceVariant:
    sequence: str
    reference_type: str
    side: str
    source_id: str
    transform: str
    strand: str
    sign: str
    direction: str

    @property
    def metadata(self) -> Dict[str, str]:
        return {
            "side": self.side,
            "source_id": self.source_id,
            "transform": self.transform,
            "strand": self.strand,
            "sign": self.sign,
            "direction": self.direction,
        }


@dataclass(frozen=True)
class ReferenceBundle:
    source_path: str
    edge: str
    edge_by_side: Dict[str, str]
    adapters: Dict[str, str]
    barcode_cores: Dict[str, Dict[str, str]]
    variants: List[ReferenceVariant]


def build_reference_bundle_from_config(config_path: str = "config.yaml") -> ReferenceBundle:
    config = load_config(config_path)
    base_dir = os.path.dirname(os.path.abspath(config_path))
    return build_reference_bundle(config, base_dir=base_dir)


def build_reference_bundle(config: Mapping[str, Any], base_dir: str = ".") -> ReferenceBundle:
    builder_config = config.get("reference_builder", {})
    if not builder_config:
        raise ValueError("Missing reference_builder section in config.yaml")

    fasta_path = builder_config.get(
        "fasta_path",
        os.path.join("Read_file", "barcodes.fasta"),
    )
    fasta_path = _resolve_path(str(fasta_path), base_dir)

    edge = str(builder_config.get("edge", "auto")).lower()
    if edge not in VALID_EDGES:
        raise ValueError(f"reference_builder.edge must be one of {VALID_EDGES}, got {edge!r}")

    groups_config = builder_config.get("groups", {})
    if not groups_config:
        raise ValueError("reference_builder.groups must define at least left and right groups")
    missing_groups = [side for side in ("left", "right") if side not in groups_config]
    if missing_groups:
        raise ValueError(f"reference_builder.groups is missing required groups: {missing_groups}")

    records = {
        record_id: sequence.strip().upper()
        for record_id, sequence in load_fasta_dict(fasta_path).items()
    }
    grouped_records = _split_records_by_prefix(records, groups_config)

    adapters: Dict[str, str] = {}
    barcode_cores: Dict[str, Dict[str, str]] = {}
    edge_by_side: Dict[str, str] = {}

    for side, side_records in grouped_records.items():
        side_edge = str(groups_config[side].get("edge", edge)).lower()
        if side_edge not in VALID_EDGES:
            raise ValueError(f"edge for group {side!r} must be one of {VALID_EDGES}, got {side_edge!r}")

        adapter, cores, selected_edge = _extract_adapter_and_cores(side_records, side_edge)
        adapters[side] = adapter
        barcode_cores[side] = cores
        edge_by_side[side] = selected_edge

    _validate_barcode_cores(barcode_cores)

    variants: List[ReferenceVariant] = []
    for side, adapter in adapters.items():
        variants.extend(_build_variants(adapter, "adapter", side, f"{side}_adapter"))

    for side, cores in barcode_cores.items():
        for source_id, barcode_core in cores.items():
            variants.extend(_build_variants(barcode_core, "barcode_core", side, source_id))

    return ReferenceBundle(
        source_path=fasta_path,
        edge=edge,
        edge_by_side=edge_by_side,
        adapters=adapters,
        barcode_cores=barcode_cores,
        variants=variants,
    )


def _resolve_path(path: str, base_dir: str) -> str:
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(base_dir, path))


def _split_records_by_prefix(
    records: Mapping[str, str],
    groups_config: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Dict[str, str]]:
    grouped: Dict[str, Dict[str, str]] = {side: {} for side in groups_config}
    unmatched: List[str] = []

    for record_id, sequence in records.items():
        matched_sides = [
            side
            for side, side_config in groups_config.items()
            if _matches_prefix(record_id, side_config.get("prefixes", []))
        ]

        if len(matched_sides) > 1:
            raise ValueError(
                f"Record {record_id!r} matches multiple reference groups: {matched_sides}"
            )
        if not matched_sides:
            unmatched.append(record_id)
            continue

        grouped[matched_sides[0]][record_id] = sequence

    empty_groups = [side for side, side_records in grouped.items() if not side_records]
    if empty_groups:
        raise ValueError(f"No FASTA records matched reference groups: {empty_groups}")
    if unmatched:
        raise ValueError(f"FASTA records did not match any reference group: {unmatched}")

    return grouped


def _matches_prefix(record_id: str, prefixes: Iterable[str]) -> bool:
    prefix_list = [str(prefix) for prefix in prefixes]
    if not prefix_list:
        raise ValueError("Each reference group must define at least one prefix")
    return any(record_id.startswith(prefix) for prefix in prefix_list)


def _extract_adapter_and_cores(
    records: Mapping[str, str],
    edge: str,
) -> Tuple[str, Dict[str, str], str]:
    sequences = list(records.values())

    prefix = _longest_common_prefix(sequences)
    suffix = _longest_common_suffix(sequences)

    if edge == "prefix":
        selected_edge = "prefix"
        adapter = prefix
    elif edge == "suffix":
        selected_edge = "suffix"
        adapter = suffix
    else:
        selected_edge = "suffix" if len(suffix) >= len(prefix) else "prefix"
        adapter = suffix if selected_edge == "suffix" else prefix

    if not adapter:
        raise ValueError("Could not find a non-empty common adapter edge with 100% identity")

    if selected_edge == "prefix":
        cores = {
            record_id: sequence[len(adapter):]
            for record_id, sequence in records.items()
        }
    else:
        cores = {
            record_id: sequence[:-len(adapter)]
            for record_id, sequence in records.items()
        }

    return adapter, cores, selected_edge


def _longest_common_prefix(sequences: List[str]) -> str:
    if not sequences:
        return ""

    shortest = min(sequences, key=len)
    for index, base in enumerate(shortest):
        if any(sequence[index] != base for sequence in sequences):
            return shortest[:index]
    return shortest


def _longest_common_suffix(sequences: List[str]) -> str:
    reversed_sequences = [sequence[::-1] for sequence in sequences]
    return _longest_common_prefix(reversed_sequences)[::-1]


def _validate_barcode_cores(barcode_cores: Mapping[str, Mapping[str, str]]) -> None:
    seen: Dict[str, str] = {}

    for side, cores in barcode_cores.items():
        for source_id, barcode_core in cores.items():
            if not barcode_core:
                raise ValueError(f"barcode_core is empty for {side}:{source_id}")
            if barcode_core in seen:
                raise ValueError(
                    "barcode_core values must be unique; "
                    f"{side}:{source_id} duplicates {seen[barcode_core]}"
                )
            seen[barcode_core] = f"{side}:{source_id}"


def _build_variants(
    sequence: str,
    reference_type: str,
    side: str,
    source_id: str,
) -> List[ReferenceVariant]:
    return [
        _make_variant(sequence, reference_type, side, source_id, transform)
        for transform in TRANSFORMS
    ]


def _make_variant(
    sequence: str,
    reference_type: str,
    side: str,
    source_id: str,
    transform: str,
) -> ReferenceVariant:
    transformed, strand, sign, direction = _transform_sequence(sequence, transform)
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


def _transform_sequence(sequence: str, transform: str) -> Tuple[str, str, str, str]:
    if transform == "forward":
        return sequence, "forward", "+", "5_to_3"
    if transform == "reverse":
        return sequence[::-1], "forward", "+", "3_to_5"
    if transform == "reverse_complement":
        return sequence.translate(DNA_COMPLEMENT)[::-1], "reverse", "-", "5_to_3"
    if transform == "complement":
        return sequence.translate(DNA_COMPLEMENT), "reverse", "-", "3_to_5"
    raise ValueError(f"Unsupported transform: {transform}")
