from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from data_loader import write_primers_fasta


PASS_STATUS = "Pass"
DEFAULT_OUTPUT_PATH = "Read_file/primer.fasta"
DEFAULT_PRIMER_LENGTH = 20
DEFAULT_WINDOW = 40


def discover_primers(
    read_results: Iterable[Mapping[str, Any]],
    config: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """
    Discover primer candidates from confidently parsed reads and write them to FASTA.

    The function is intentionally independent from main.py. It consumes the in-memory
    result dictionaries produced by process_read, uses primer_discovery config values,
    and never changes trimming coordinates or trimmed sequences.
    """
    discovery_config = _discovery_config(config)
    output_path = _config_value(
        discovery_config,
        ("output_path", "fasta_path", "path"),
        DEFAULT_OUTPUT_PATH,
    )

    stats = {
        "total_reads": 0,
        "confident_reads": 0,
        "skipped_reads": 0,
        "candidate_count": 0,
    }

    if discovery_config.get("enabled") is False:
        write_primers_fasta([], str(output_path))
        return {
            "enabled": False,
            "output_path": str(output_path),
            "primers": {},
            "records": [],
            "clusters": {"left": [], "right": []},
            "stats": stats,
        }

    observations_by_side: dict[str, list[dict[str, Any]]] = {"left": [], "right": []}
    for result in read_results:
        stats["total_reads"] += 1
        if not _is_confident_result(result, discovery_config):
            stats["skipped_reads"] += 1
            continue

        read_seq = _read_sequence(result)
        trim_bounds = _trim_bounds(result)
        if read_seq is None or trim_bounds is None:
            stats["skipped_reads"] += 1
            continue

        stats["confident_reads"] += 1
        trim_start, trim_end = trim_bounds
        for side in ("left", "right"):
            observations = _candidate_observations(
                read_seq=read_seq,
                read_id=str(result.get("read_id") or ""),
                trim_start=trim_start,
                trim_end=trim_end,
                side=side,
                discovery_config=discovery_config,
            )
            observations_by_side[side].extend(observations)
            stats["candidate_count"] += len(observations)

    cluster_threshold = int(
        _config_value(
            discovery_config,
            ("edit_distance_threshold", "cluster_distance", "max_edit_distance"),
            1,
        )
    )
    records: list[dict[str, str]] = []
    clusters_for_return: dict[str, list[dict[str, Any]]] = {"left": [], "right": []}

    for side in ("left", "right"):
        clusters = _cluster_observations(
            observations_by_side[side],
            max_distance=max(0, cluster_threshold),
        )
        selected = _select_clusters(clusters, side, discovery_config)
        clusters_for_return[side] = [_cluster_summary(cluster) for cluster in selected]

        multiple = len(selected) > 1
        for index, cluster in enumerate(selected, start=1):
            record_id = f"{side}_primer_{index}" if multiple else f"{side}_primer"
            records.append({"id": record_id, "sequence": cluster["consensus"]})

    write_primers_fasta(records, str(output_path))
    return {
        "enabled": True,
        "output_path": str(output_path),
        "primers": {record["id"]: record["sequence"] for record in records},
        "records": records,
        "clusters": clusters_for_return,
        "stats": stats,
    }


def _discovery_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    if not isinstance(config, Mapping):
        return {}

    if "primer_discovery" in config:
        raw_config = config.get("primer_discovery") or {}
        return dict(raw_config) if isinstance(raw_config, Mapping) else {}

    discovery_keys = {
        "enabled",
        "left_window",
        "right_window",
        "primer_length",
        "length",
        "candidate_lengths",
        "max_total_penalty",
        "edit_distance_threshold",
        "output_path",
    }
    if any(key in config for key in discovery_keys):
        return dict(config)

    return {}


def _is_confident_result(
    result: Mapping[str, Any],
    discovery_config: Mapping[str, Any],
) -> bool:
    required_status = str(
        _config_value(discovery_config, ("status", "required_status"), PASS_STATUS)
    )
    if str(result.get("status")) != required_status:
        return False

    max_total_penalty = _config_value(
        discovery_config,
        ("max_total_penalty", "total_penalty_max", "penalty_threshold"),
        0.0,
    )
    if max_total_penalty is None:
        return True

    total_penalty = _result_total_penalty(result)
    if total_penalty is None:
        return False
    return total_penalty <= float(max_total_penalty)


def _result_total_penalty(result: Mapping[str, Any]) -> float | None:
    if result.get("total_penalty") is not None:
        try:
            return float(result["total_penalty"])
        except (TypeError, ValueError):
            return None

    if result.get("score") is not None:
        try:
            return -float(result["score"])
        except (TypeError, ValueError):
            return None

    return None


def _read_sequence(result: Mapping[str, Any]) -> str | None:
    sequence = (
        result.get("original_seq")
        or result.get("read_seq")
        or result.get("sequence")
        or result.get("seq")
    )
    if not sequence:
        return None
    return str(sequence).upper()


def _trim_bounds(result: Mapping[str, Any]) -> tuple[int, int] | None:
    try:
        trim_start = int(result["trim_start"])
        trim_end = int(result["trim_end"])
    except (KeyError, TypeError, ValueError):
        return None

    if trim_start < 0 or trim_end <= trim_start:
        return None
    return trim_start, trim_end


def _candidate_observations(
    read_seq: str,
    read_id: str,
    trim_start: int,
    trim_end: int,
    side: str,
    discovery_config: Mapping[str, Any],
) -> list[dict[str, Any]]:
    window_size = _side_int(discovery_config, side, "window", DEFAULT_WINDOW)
    candidate_lengths = _candidate_lengths(discovery_config, side)
    allow_ambiguous = bool(discovery_config.get("allow_ambiguous", False))
    require_unique = bool(discovery_config.get("require_unique_in_read", True))

    window_start, window_end = _candidate_window(
        trim_start=trim_start,
        trim_end=trim_end,
        side=side,
        window_size=max(0, window_size),
    )
    observations: list[dict[str, Any]] = []

    for length in candidate_lengths:
        if length <= 0 or window_end - window_start < length:
            continue
        for start in range(window_start, window_end - length + 1):
            end = start + length
            candidate = read_seq[start:end].upper()
            if not _valid_candidate(candidate, allow_ambiguous=allow_ambiguous):
                continue
            if require_unique and _count_occurrences(read_seq, candidate) > 1:
                continue
            observations.append(
                {
                    "sequence": candidate,
                    "read_id": read_id,
                    "side": side,
                    "start": start,
                    "end": end,
                    "offset": _edge_offset(start, end, trim_start, trim_end, side),
                }
            )

    return observations


def _candidate_window(
    trim_start: int,
    trim_end: int,
    side: str,
    window_size: int,
) -> tuple[int, int]:
    if side == "left":
        return trim_start, min(trim_end, trim_start + window_size)
    return max(trim_start, trim_end - window_size), trim_end


def _edge_offset(
    start: int,
    end: int,
    trim_start: int,
    trim_end: int,
    side: str,
) -> int:
    if side == "left":
        return start - trim_start
    return trim_end - end


def _candidate_lengths(
    discovery_config: Mapping[str, Any],
    side: str,
) -> list[int]:
    side_config = discovery_config.get(side)
    side_mapping = side_config if isinstance(side_config, Mapping) else {}
    raw_lengths = _config_value(
        side_mapping,
        ("lengths", "candidate_lengths"),
        None,
    )
    if raw_lengths is None:
        raw_lengths = _config_value(
            discovery_config,
            (
                f"{side}_lengths",
                f"{side}_candidate_lengths",
                "lengths",
                "candidate_lengths",
            ),
            None,
        )
    if raw_lengths is not None:
        return _normalize_lengths(raw_lengths)

    raw_length = _config_value(side_mapping, ("length", "primer_length"), None)
    if raw_length is None:
        raw_length = _config_value(
            discovery_config,
            (f"{side}_length", "primer_length", "length", "candidate_length"),
            DEFAULT_PRIMER_LENGTH,
        )

    return _normalize_lengths(raw_length)


def _normalize_lengths(raw_lengths: Any) -> list[int]:
    if isinstance(raw_lengths, str):
        if "," in raw_lengths:
            items = raw_lengths.split(",")
        else:
            items = [raw_lengths]
    elif isinstance(raw_lengths, Sequence):
        items = list(raw_lengths)
    else:
        items = [raw_lengths]

    lengths = sorted({int(item) for item in items if item is not None})
    return [length for length in lengths if length > 0]


def _valid_candidate(candidate: str, allow_ambiguous: bool) -> bool:
    alphabet = set("ACGTN") if allow_ambiguous else set("ACGT")
    if not candidate or set(candidate) - alphabet:
        return False
    return bool(set(candidate) - {"N"})


def _count_occurrences(sequence: str, needle: str) -> int:
    count = 0
    start = 0
    while True:
        index = sequence.find(needle, start)
        if index < 0:
            return count
        count += 1
        start = index + 1


def _cluster_observations(
    observations: Sequence[Mapping[str, Any]],
    max_distance: int,
) -> list[dict[str, Any]]:
    sequence_stats = _sequence_stats(observations)
    clusters: list[dict[str, Any]] = []

    for sequence, stats in sorted(
        sequence_stats.items(),
        key=lambda item: (
            -len(item[1]["read_ids"]),
            _mean(item[1]["offsets"]),
            -item[1]["count"],
            item[0],
        ),
    ):
        target_cluster = None
        for cluster in clusters:
            if _edit_distance(sequence, cluster["representative"], max_distance) <= max_distance:
                target_cluster = cluster
                break

        if target_cluster is None:
            target_cluster = {
                "representative": sequence,
                "members": Counter(),
                "read_ids": set(),
                "offsets": [],
                "observations": 0,
            }
            clusters.append(target_cluster)

        target_cluster["members"][sequence] += stats["count"]
        target_cluster["read_ids"].update(stats["read_ids"])
        target_cluster["offsets"].extend(stats["offsets"])
        target_cluster["observations"] += stats["count"]
        target_cluster["representative"] = _best_representative(target_cluster)

    for cluster in clusters:
        cluster["consensus"] = _consensus_sequence(cluster["members"])
        cluster["support"] = len(cluster["read_ids"])
        cluster["mean_offset"] = _mean(cluster["offsets"])

    return sorted(clusters, key=_cluster_sort_key)


def _sequence_stats(observations: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    stats: dict[str, dict[str, Any]] = {}
    for observation in observations:
        sequence = str(observation["sequence"])
        record = stats.setdefault(
            sequence,
            {"count": 0, "read_ids": set(), "offsets": []},
        )
        record["count"] += 1
        record["read_ids"].add(str(observation.get("read_id") or ""))
        record["offsets"].append(int(observation.get("offset", 0)))
    return stats


def _select_clusters(
    clusters: Sequence[Mapping[str, Any]],
    side: str,
    discovery_config: Mapping[str, Any],
) -> list[dict[str, Any]]:
    min_support = _side_int(discovery_config, side, "min_support", 1)
    max_per_side = _side_value(discovery_config, side, "max_per_side", 1)
    filtered = [
        dict(cluster)
        for cluster in clusters
        if int(cluster.get("support", 0)) >= min_support
    ]

    if max_per_side is None:
        return filtered

    max_per_side = int(max_per_side)
    if max_per_side < 1:
        return []
    return filtered[:max_per_side]


def _cluster_summary(cluster: Mapping[str, Any]) -> dict[str, Any]:
    members = cluster.get("members", Counter())
    if isinstance(members, Counter):
        member_count = len(members)
    else:
        member_count = len(dict(members))

    return {
        "sequence": cluster.get("consensus"),
        "support": cluster.get("support", 0),
        "observations": cluster.get("observations", 0),
        "mean_offset": cluster.get("mean_offset", 0.0),
        "member_count": member_count,
    }


def _best_representative(cluster: Mapping[str, Any]) -> str:
    members = cluster["members"]
    return sorted(members, key=lambda sequence: (-members[sequence], sequence))[0]


def _consensus_sequence(members: Counter[str]) -> str:
    if not members:
        return ""

    lengths = {len(sequence) for sequence in members}
    if len(lengths) != 1:
        return sorted(members, key=lambda sequence: (-members[sequence], sequence))[0]

    length = next(iter(lengths))
    consensus_bases: list[str] = []
    for index in range(length):
        base_counts: Counter[str] = Counter()
        for sequence, count in members.items():
            base_counts[sequence[index]] += count
        consensus_bases.append(
            sorted(base_counts, key=lambda base: (-base_counts[base], base))[0]
        )
    return "".join(consensus_bases)


def _cluster_sort_key(cluster: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        -int(cluster.get("support", 0)),
        float(cluster.get("mean_offset", 0.0)),
        -int(cluster.get("observations", 0)),
        str(cluster.get("consensus") or ""),
    )


def _edit_distance(left: str, right: str, max_distance: int | None = None) -> int:
    if left == right:
        return 0
    if max_distance is not None and abs(len(left) - len(right)) > max_distance:
        return max_distance + 1

    previous = list(range(len(right) + 1))
    for left_index, left_base in enumerate(left, start=1):
        current = [left_index]
        row_min = current[0]
        for right_index, right_base in enumerate(right, start=1):
            cost = 0 if left_base == right_base else 1
            value = min(
                previous[right_index] + 1,
                current[right_index - 1] + 1,
                previous[right_index - 1] + cost,
            )
            current.append(value)
            row_min = min(row_min, value)
        if max_distance is not None and row_min > max_distance:
            return max_distance + 1
        previous = current
    return previous[-1]


def _side_int(
    config: Mapping[str, Any],
    side: str,
    key: str,
    default: int,
) -> int:
    return int(_side_value(config, side, key, default))


def _side_value(
    config: Mapping[str, Any],
    side: str,
    key: str,
    default: Any,
) -> Any:
    side_config = config.get(side)
    if isinstance(side_config, Mapping) and key in side_config:
        return side_config[key]

    side_key = f"{side}_{key}"
    if side_key in config:
        return config[side_key]

    return config.get(key, default)


def _config_value(
    config: Mapping[str, Any],
    keys: Sequence[str],
    default: Any,
) -> Any:
    for key in keys:
        if key in config:
            return config[key]
    return default


def _mean(values: Sequence[int]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)
