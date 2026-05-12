import os
import re
from collections import defaultdict
from collections.abc import Sequence
from typing import Any

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


PASS_STATUS = "Pass"


def export_demultiplexed_reads(
    reads_data_list,
    output_dir="Demux_Out",
    split_unassigned_by_reason=False,
):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    grouped_reads = defaultdict(list)
    unassigned_reads = defaultdict(list) if split_unassigned_by_reason else []

    for data in reads_data_list:
        status = _clean_text(data.get("status")) or "Discarded"
        sample_id = _clean_text(data.get("sample_id"))
        barcode_id = _clean_text(data.get("barcode_id"))
        output_id = sample_id or barcode_id

        if status == PASS_STATUS and output_id:
            trimmed_seq = data.get("trimmed_seq", "")
            trimmed_quality = data.get("trimmed_quality")
            if _quality_matches_sequence(trimmed_seq, trimmed_quality):
                description_parts = []
                if sample_id:
                    description_parts.append(f"sample={sample_id}")
                if barcode_id:
                    description_parts.append(f"barcode={barcode_id}")

                record = SeqRecord(
                    Seq(trimmed_seq),
                    id=str(data.get("read_id", "")),
                    description=" ".join(description_parts),
                )
                record.letter_annotations["phred_quality"] = list(trimmed_quality)
                grouped_reads[output_id].append(record)
                continue

            _append_unassigned(
                unassigned_reads,
                data,
                reason="QualityLengthMismatch",
                status=status,
                split_by_reason=split_unassigned_by_reason,
            )
            continue

        reason = _unassigned_reason(data, status, output_id)
        _append_unassigned(
            unassigned_reads,
            data,
            reason=reason,
            status=status,
            split_by_reason=split_unassigned_by_reason,
        )

    for output_id, records in grouped_reads.items():
        file_path = os.path.join(output_dir, f"{_safe_filename(output_id)}.fastq")
        SeqIO.write(records, file_path, "fastq")

    if split_unassigned_by_reason:
        for reason, records in unassigned_reads.items():
            if records:
                filename = f"unassigned_{_safe_filename(reason)}.fastq"
                SeqIO.write(records, os.path.join(output_dir, filename), "fastq")
    elif unassigned_reads:
        unassigned_path = os.path.join(output_dir, "unassigned_reads.fastq")
        SeqIO.write(unassigned_reads, unassigned_path, "fastq")


def _append_unassigned(
    unassigned_reads,
    data: dict[str, Any],
    reason: str,
    status: str,
    split_by_reason: bool,
) -> None:
    record = _unassigned_record(data, reason=reason, status=status)
    if split_by_reason:
        unassigned_reads[reason].append(record)
    else:
        unassigned_reads.append(record)


def _unassigned_record(data: dict[str, Any], reason: str, status: str) -> SeqRecord:
    record = SeqRecord(
        Seq(data.get("original_seq", "")),
        id=str(data.get("read_id", "")),
        description=f"unassigned reason={_description_value(reason)} status={_description_value(status)}",
    )
    original_quality = data.get("original_quality")
    if original_quality is not None:
        record.letter_annotations["phred_quality"] = list(original_quality)
    return record


def _quality_matches_sequence(sequence: Any, quality: Any) -> bool:
    if quality is None:
        return False
    if not isinstance(quality, Sequence) or isinstance(quality, (str, bytes, bytearray)):
        return False
    return len(str(sequence)) == len(quality)


def _unassigned_reason(data: dict[str, Any], status: str, output_id: str | None) -> str:
    if status == PASS_STATUS and not output_id:
        return "MissingSampleAndBarcode"
    return (
        _clean_text(data.get("reason"))
        or _clean_text(data.get("rejection_reason"))
        or _clean_text(data.get("diagnostic_message"))
        or status
        or "Unassigned"
    )


def _clean_text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return None
    return text


def _description_value(value: Any) -> str:
    text = _clean_text(value) or "Unassigned"
    return re.sub(r"\s+", "_", text)


def _safe_filename(value: Any) -> str:
    text = _description_value(value)
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._") or "Unknown"
