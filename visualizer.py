import os
import tempfile
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


MAX_BAR_CATEGORIES = 30
PASS_STATUS = "Pass"
PENALTY_COLUMNS = ["adapter_penalty", "barcode_penalty", "primer_penalty"]


def generate_report(csv_log_path, output_dir):
    """
    Reads the demultiplexing CSV log and generates visualization reports.

    The report supports both the current schema with total_penalty and penalty
    breakdown columns, and old results.csv files that only have score/status/barcode_id.
    """
    os.makedirs(output_dir, exist_ok=True)

    try:
        df = pd.read_csv(csv_log_path)
    except FileNotFoundError:
        print(f"Error: File not found at {csv_log_path}")
        return
    except pd.errors.EmptyDataError:
        print(f"Error: The file {csv_log_path} is empty.")
        return

    sns.set_theme(style="whitegrid", palette="muted")

    generated = 0
    generated += _generate_status_counts(df, output_dir)
    generated += _generate_barcode_pair_counts(df, output_dir)
    generated += _generate_penalty_distribution(df, output_dir)
    generated += _generate_rejection_reasons(df, output_dir)
    generated += _generate_penalty_breakdown(df, output_dir)

    print(f"Successfully generated {generated} report files in '{output_dir}'.")


def _generate_status_counts(df: pd.DataFrame, output_dir: str) -> int:
    if "status" not in df.columns:
        return 0

    counts = _counts_table(df["status"], "status")
    if counts.empty:
        return 0

    generated = _write_table(counts, output_dir, "Status_Counts.csv")
    generated += _plot_bar(
        counts,
        x="status",
        y="count",
        title="Read Status Counts",
        xlabel="Status",
        ylabel="Number of Reads",
        output_dir=output_dir,
        filenames=["Status_Counts.png", "Status_Distribution.png"],
        palette="pastel",
    )
    return generated


def _generate_barcode_pair_counts(df: pd.DataFrame, output_dir: str) -> int:
    if not {"left_barcode_id", "right_barcode_id", "barcode_id"} & set(df.columns):
        return 0

    barcode_pairs = _barcode_pair_series(df)
    counts = _counts_table(barcode_pairs, "barcode_pair")
    if counts.empty:
        return 0

    generated = _write_table(counts, output_dir, "Barcode_Pair_Counts.csv")
    generated += _plot_bar(
        counts,
        x="barcode_pair",
        y="count",
        title="Barcode Pair Counts",
        xlabel="Barcode Pair",
        ylabel="Number of Reads",
        output_dir=output_dir,
        filenames=["Barcode_Pair_Counts.png", "Barcode_Distribution.png"],
        palette="viridis",
    )
    return generated


def _generate_penalty_distribution(df: pd.DataFrame, output_dir: str) -> int:
    metric_column = _penalty_metric_column(df)
    if metric_column is None:
        return 0

    values = pd.to_numeric(df[metric_column], errors="coerce").dropna()
    if values.empty:
        return 0

    label = "Total Penalty" if metric_column == "total_penalty" else "Score"
    filename = (
        "Total_Penalty_Distribution.png"
        if metric_column == "total_penalty"
        else "Score_Distribution.png"
    )
    filenames = [filename]
    if metric_column == "total_penalty":
        filenames.append("Score_Distribution.png")

    plt.figure(figsize=(10, 6))
    sns.histplot(values, bins=30, kde=True, color="steelblue")
    plt.title(f"{label} Distribution", fontsize=16, pad=15)
    plt.xlabel(label, fontsize=12)
    plt.ylabel("Number of Reads", fontsize=12)
    plt.tight_layout()
    _save_current_figure(output_dir, filenames)
    plt.close()
    return len(filenames)


def _generate_rejection_reasons(df: pd.DataFrame, output_dir: str) -> int:
    if "status" not in df.columns:
        return 0

    status = df["status"].fillna("").astype(str)
    rejected = df[status != PASS_STATUS].copy()
    if rejected.empty:
        return 0

    reasons = _rejection_reason_series(rejected)
    counts = _counts_table(reasons, "reason")
    if counts.empty:
        return 0

    generated = _write_table(counts, output_dir, "Rejection_Reasons.csv")
    generated += _plot_bar(
        counts,
        x="reason",
        y="count",
        title="Rejection Reasons",
        xlabel="Reason",
        ylabel="Number of Reads",
        output_dir=output_dir,
        filenames=["Rejection_Reasons.png"],
        palette="mako",
    )
    return generated


def _generate_penalty_breakdown(df: pd.DataFrame, output_dir: str) -> int:
    available_columns = [column for column in PENALTY_COLUMNS if column in df.columns]
    if not available_columns:
        return 0

    numeric = df[available_columns].apply(pd.to_numeric, errors="coerce")
    if numeric.dropna(how="all").empty:
        return 0

    summary = pd.DataFrame(
        {
            "component": available_columns,
            "sum": [numeric[column].sum(skipna=True) for column in available_columns],
            "mean": [numeric[column].mean(skipna=True) for column in available_columns],
            "median": [numeric[column].median(skipna=True) for column in available_columns],
            "non_null_count": [numeric[column].count() for column in available_columns],
        }
    )

    generated = _write_table(summary, output_dir, "Penalty_Breakdown.csv")

    plot_data = summary[["component", "sum"]].copy()
    plot_data["component"] = plot_data["component"].str.replace("_penalty", "", regex=False)
    generated += _plot_bar(
        plot_data,
        x="component",
        y="sum",
        title="Penalty Breakdown",
        xlabel="Component",
        ylabel="Total Penalty",
        output_dir=output_dir,
        filenames=["Penalty_Breakdown.png"],
        palette="crest",
    )
    return generated


def _penalty_metric_column(df: pd.DataFrame) -> str | None:
    if "total_penalty" in df.columns:
        return "total_penalty"
    if "score" in df.columns:
        return "score"
    return None


def _barcode_pair_series(df: pd.DataFrame) -> pd.Series:
    left = df["left_barcode_id"].map(_clean_value) if "left_barcode_id" in df.columns else None
    right = df["right_barcode_id"].map(_clean_value) if "right_barcode_id" in df.columns else None
    barcode = df["barcode_id"].map(_clean_value) if "barcode_id" in df.columns else None

    values = []
    for index in df.index:
        left_value = left.loc[index] if left is not None else None
        right_value = right.loc[index] if right is not None else None
        barcode_value = barcode.loc[index] if barcode is not None else None

        if left_value and right_value:
            values.append(f"{left_value}_{right_value}")
        else:
            values.append(barcode_value or "Unknown")

    return pd.Series(values, name="barcode_pair")


def _rejection_reason_series(df: pd.DataFrame) -> pd.Series:
    status = df["status"].map(_clean_value) if "status" in df.columns else None
    reasons = []

    for index in df.index:
        reason = None
        for column in ("reason", "rejection_reason", "diagnostic_message"):
            if column in df.columns:
                reason = _clean_value(df.at[index, column])
                if reason:
                    break
        if not reason and status is not None:
            reason = status.loc[index]
        reasons.append(reason or "Unknown")

    return pd.Series(reasons, name="reason")


def _counts_table(values: pd.Series, column_name: str) -> pd.DataFrame:
    cleaned = values.map(_clean_value).fillna("Unknown")
    counts = cleaned.value_counts(dropna=False).rename_axis(column_name).reset_index(name="count")
    return counts


def _write_table(df: pd.DataFrame, output_dir: str, filename: str) -> int:
    df.to_csv(os.path.join(output_dir, filename), index=False)
    return 1


def _plot_bar(
    counts: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    xlabel: str,
    ylabel: str,
    output_dir: str,
    filenames: list[str],
    palette: str,
) -> int:
    plot_data = counts.head(MAX_BAR_CATEGORIES)
    if plot_data.empty:
        return 0

    width = max(10, min(18, len(plot_data) * 0.65))
    plt.figure(figsize=(width, 6))
    sns.barplot(data=plot_data, x=x, y=y, hue=x, palette=palette, legend=False)
    plt.title(title, fontsize=16, pad=15)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.xticks(rotation=45, ha="right")
    _annotate_bars()
    plt.tight_layout()
    _save_current_figure(output_dir, filenames)
    plt.close()
    return len(filenames)


def _annotate_bars() -> None:
    ax = plt.gca()
    for patch in ax.patches:
        height = patch.get_height()
        if pd.isna(height):
            continue
        label = f"{height:.2f}" if height % 1 else f"{int(height)}"
        ax.annotate(
            label,
            (patch.get_x() + patch.get_width() / 2.0, height),
            ha="center",
            va="center",
            fontsize=10,
            color="black",
            xytext=(0, 5),
            textcoords="offset points",
        )


def _save_current_figure(output_dir: str, filenames: list[str]) -> None:
    for filename in filenames:
        plt.savefig(os.path.join(output_dir, filename), dpi=300)


def _clean_value(value: Any) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return None
    return text
