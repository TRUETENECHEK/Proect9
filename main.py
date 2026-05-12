import csv
import os
from collections.abc import Callable, Iterable
from typing import Any

from data_loader import load_config, load_mapping_excel, read_fastq_generator
from demux_logic import RESULT_FIELDNAMES, process_read
from export_fastq import export_demultiplexed_reads
from reference_builder import build_reference_bundle

RICH_AVAILABLE = False
TQDM_AVAILABLE = False
console = None

try:
    from rich.console import Console
    from rich.progress import (
        Progress,
        SpinnerColumn,
        TextColumn,
        TimeElapsedColumn,
    )
    import rich.traceback

    rich.traceback.install()
    RICH_AVAILABLE = True
    console = Console()
except ImportError:
    try:
        from tqdm import tqdm

        TQDM_AVAILABLE = True
    except ImportError:
        pass


FASTQ_OUTPUT_DIR = "Demux_Out"
MAPPING_PATH = "mapping.xlsx"


def run_pipeline(
    config_path: str,
    fastq_path: str,
    output_csv: str,
    output_dir: str,
    sampling_rate: int | None = None,
) -> int | None:
    if sampling_rate is not None and sampling_rate <= 0:
        raise ValueError("sampling_rate must be a positive integer or None")

    _print_start(sampling_rate)

    try:
        config = load_config(config_path)
        _print_success(f"Loaded config from {config_path}")
    except Exception as exc:
        _print_error(f"Error loading config: {exc}")
        return None

    try:
        base_dir = os.path.dirname(os.path.abspath(config_path))
        reference_bundle = build_reference_bundle(config, base_dir=base_dir)
        _print_success("Built ReferenceBundle")
    except Exception as exc:
        _print_error(f"Error building ReferenceBundle: {exc}")
        return None

    adapters = config.get("sequences", {})
    mapping_dict = load_mapping_excel(MAPPING_PATH)
    export_records: list[dict[str, Any]] = []

    try:
        with open(output_csv, mode="w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDNAMES, extrasaction="ignore")
            writer.writeheader()

            read_records = _sampled_records(read_fastq_generator(fastq_path), sampling_rate)
            processed_count = _process_with_progress(
                records=read_records,
                fastq_path=fastq_path,
                sampling_rate=sampling_rate,
                writer=writer,
                export_records=export_records,
                config=config,
                adapters=adapters,
                mapping_dict=mapping_dict,
                reference_bundle=reference_bundle,
            )
    except Exception as exc:
        _print_error(f"Error reading FASTQ or writing {output_csv}: {exc}")
        return None

    _print_success(f"Finished processing. Results saved to {output_csv}")
    _print_info("Exporting demultiplexed FASTQ files")
    export_demultiplexed_reads(export_records, output_dir=FASTQ_OUTPUT_DIR)

    _print_info(f"Generating report in '{output_dir}' directory")
    try:
        _generate_report(output_csv, output_dir)
        _print_success("Report generation complete")
    except Exception as exc:
        _print_error(f"Error generating report: {exc}")

    return processed_count


def _process_with_progress(
    records: Iterable[tuple[str, str, list[int]]],
    fastq_path: str,
    sampling_rate: int | None,
    writer: csv.DictWriter,
    export_records: list[dict[str, Any]],
    config: dict[str, Any],
    adapters: Any,
    mapping_dict: dict[str, str],
    reference_bundle: Any,
) -> int:
    if RICH_AVAILABLE:
        assert console is not None
        _print_info(_processing_message(fastq_path, sampling_rate))
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            TimeElapsedColumn(),
            console=console,
        ) as progress:
            task = progress.add_task("Demultiplexing")

            def advance(processed_count: int) -> None:
                description = f"Demultiplexing ({processed_count} reads)"
                progress.update(task, advance=1, description=description)

            return _process_records(
                records=records,
                writer=writer,
                export_records=export_records,
                config=config,
                adapters=adapters,
                mapping_dict=mapping_dict,
                reference_bundle=reference_bundle,
                on_processed=advance,
            )

    if TQDM_AVAILABLE:
        wrapped_records = tqdm(  # type: ignore[name-defined]
            records,
            desc="Processing Reads" if sampling_rate is None else "Subsampling Reads",
            unit=" read",
        )
        return _process_records(
            records=wrapped_records,
            writer=writer,
            export_records=export_records,
            config=config,
            adapters=adapters,
            mapping_dict=mapping_dict,
            reference_bundle=reference_bundle,
        )

    _print_info(_processing_message(fastq_path, sampling_rate))
    print_every = 100 if sampling_rate is None else 10

    def report(processed_count: int) -> None:
        if processed_count % print_every == 0:
            suffix = "" if sampling_rate is None else " sampled"
            print(f"Processed {processed_count}{suffix} reads...")

    return _process_records(
        records=records,
        writer=writer,
        export_records=export_records,
        config=config,
        adapters=adapters,
        mapping_dict=mapping_dict,
        reference_bundle=reference_bundle,
        on_processed=report,
    )


def _process_records(
    records: Iterable[tuple[str, str, list[int]]],
    writer: csv.DictWriter,
    export_records: list[dict[str, Any]],
    config: dict[str, Any],
    adapters: Any,
    mapping_dict: dict[str, str],
    reference_bundle: Any,
    on_processed: Callable[[int], None] | None = None,
) -> int:
    processed_count = 0

    for read_id, sequence, phred_quality in records:
        result = _process_one_read(
            read_id=read_id,
            sequence=sequence,
            phred_quality=phred_quality,
            adapters=adapters,
            config=config,
            mapping_dict=mapping_dict,
            reference_bundle=reference_bundle,
        )
        writer.writerow(result)
        export_records.append(_fastq_export_data(result))

        processed_count += 1
        if on_processed is not None:
            on_processed(processed_count)

    return processed_count


def _process_one_read(
    read_id: str,
    sequence: str,
    phred_quality: list[int],
    adapters: Any,
    config: dict[str, Any],
    mapping_dict: dict[str, str],
    reference_bundle: Any,
) -> dict[str, Any]:
    try:
        return process_read(
            read_id=read_id,
            read_seq=sequence,
            adapters=adapters,
            barcodes=None,
            config=config,
            mapping_dict=mapping_dict,
            max_errors=2,
            quality=phred_quality,
            reference_bundle=reference_bundle,
        )
    except Exception as exc:
        return _exception_result(read_id, sequence, phred_quality, exc)


def _exception_result(
    read_id: str,
    sequence: str,
    phred_quality: list[int],
    exc: Exception,
) -> dict[str, Any]:
    return {
        "read_id": read_id,
        "barcode_id": "Error",
        "sample_id": None,
        "score": 0.0,
        "status": "Fail_Exception",
        "trimmed_seq": "",
        "trimmed_quality": None,
        "original_seq": sequence,
        "original_quality": phred_quality,
        "left_barcode_id": None,
        "right_barcode_id": None,
        "adapter_penalty": 0.0,
        "barcode_penalty": 0.0,
        "primer_penalty": 0.0,
        "layout_penalty": 0.0,
        "total_penalty": 0.0,
        "trim_start": None,
        "trim_end": None,
        "diagnostic_message": str(exc),
    }


def _fastq_export_data(result: dict[str, Any]) -> dict[str, Any]:
    return {
        "read_id": result.get("read_id"),
        "status": result.get("status"),
        "barcode_id": result.get("barcode_id"),
        "sample_id": result.get("sample_id"),
        "diagnostic_message": result.get("diagnostic_message"),
        "trimmed_seq": result.get("trimmed_seq", ""),
        "trimmed_quality": result.get("trimmed_quality"),
        "original_seq": result.get("original_seq", ""),
        "original_quality": result.get("original_quality"),
    }


def _sampled_records(
    records: Iterable[tuple[str, str, list[int]]],
    sampling_rate: int | None,
) -> Iterable[tuple[str, str, list[int]]]:
    for index, record in enumerate(records):
        if sampling_rate is None or index % sampling_rate == 0:
            yield record


def _processing_message(fastq_path: str, sampling_rate: int | None) -> str:
    if sampling_rate is None:
        return f"Processing reads from {fastq_path}"
    return f"Processing every {sampling_rate}th read from {fastq_path}"


def _generate_report(output_csv: str, output_dir: str) -> None:
    from visualizer import generate_report

    generate_report(output_csv, output_dir)


def _print_start(sampling_rate: int | None) -> None:
    if sampling_rate is None:
        message = "Starting Demultiplexing Pipeline"
    else:
        message = f"Starting Subsampled Demultiplexing Pipeline (every {sampling_rate}th read)"

    if RICH_AVAILABLE:
        assert console is not None
        console.rule(f"[bold green]{message}[/bold green]")
    else:
        print(f"{message}...")


def _print_success(message: str) -> None:
    if RICH_AVAILABLE:
        assert console is not None
        console.print(f"[green]+[/green] {message}")
    else:
        print(f"[+] {message}")


def _print_info(message: str) -> None:
    if RICH_AVAILABLE:
        assert console is not None
        console.print(f"[cyan]i[/cyan] {message}")
    else:
        print(f"[+] {message}")


def _print_error(message: str) -> None:
    if RICH_AVAILABLE:
        assert console is not None
        console.print(f"[red]-[/red] {message}")
    else:
        print(f"[-] {message}")


def main() -> None:
    run_pipeline(
        config_path="config.yaml",
        fastq_path=os.path.join("Read_file", "test_reads.fastq"),
        output_csv="results.csv",
        output_dir="report_plots",
    )


if __name__ == "__main__":
    main()
