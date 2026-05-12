from main import run_pipeline


SAMPLING_RATE = 1000


def main() -> None:
    run_pipeline(
        config_path="config.yaml",
        fastq_path="Read_file/test_reads.fastq",
        output_csv="results_subsampled.csv",
        output_dir="report_plots_subsampled",
        sampling_rate=SAMPLING_RATE,
    )


if __name__ == "__main__":
    main()
