from __future__ import annotations

from dataclasses import dataclass

import regex


_COMPLEMENT_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class MatchRecord:
    element_type: str | None
    ref_id: str | None
    sequence: str
    start: int
    end: int
    errors: tuple[int, int, int]
    penalty: int
    transform: str
    strand: str
    direction: str


def complement(seq):
    """
    Возвращает комплементарную ДНК последовательность.
    """
    return seq.translate(_COMPLEMENT_TABLE).upper()


def reverse_sequence(seq):
    """
    Возвращает последовательность в обратном порядке без комплементации.
    """
    return seq[::-1].upper()


def reverse_complement(seq):
    """
    Возвращает обратное дополнение (Reverse Complement) ДНК последовательности.
    """
    return complement(reverse_sequence(seq))


def get_reverse_complement(seq):
    """
    Старое имя оставлено для обратной совместимости.
    """
    return reverse_complement(seq)


def _penalty_from_errors(errors):
    substitutions, insertions, deletions = errors
    return substitutions + insertions + deletions


def _to_optional_str(value):
    if value is None:
        return None
    return str(value)


def _variant_from_mapping(ref_id, variant):
    sequence = (
        variant.get("sequence")
        or variant.get("seq")
        or variant.get("pattern")
        or variant.get("value")
    )

    if not sequence:
        return None

    return {
        "element_type": _to_optional_str(variant.get("element_type") or variant.get("type")),
        "ref_id": _to_optional_str(
            variant.get("ref_id") or variant.get("id") or variant.get("name") or ref_id
        ),
        "sequence": str(sequence),
    }


def _iter_reference_variants(reference_variants):
    if isinstance(reference_variants, str):
        if reference_variants:
            yield {"element_type": None, "ref_id": None, "sequence": reference_variants}
        return

    if isinstance(reference_variants, dict):
        if any(key in reference_variants for key in ("sequence", "seq", "pattern", "value")):
            variant = _variant_from_mapping(None, reference_variants)
            if variant:
                yield variant
            return

        for ref_id, sequence in reference_variants.items():
            if isinstance(sequence, dict):
                variant = _variant_from_mapping(ref_id, sequence)
                if variant:
                    yield variant
            elif sequence:
                yield {
                    "element_type": None,
                    "ref_id": _to_optional_str(ref_id),
                    "sequence": str(sequence),
                }
        return

    for index, variant in enumerate(reference_variants):
        if isinstance(variant, str):
            if variant:
                yield {"element_type": None, "ref_id": None, "sequence": variant}
            continue

        if isinstance(variant, dict):
            parsed = _variant_from_mapping(None, variant)
            if parsed:
                yield parsed
            continue

        sequence = getattr(variant, "sequence", None) or getattr(variant, "seq", None)
        if sequence:
            yield {
                "element_type": _to_optional_str(getattr(variant, "element_type", None)),
                "ref_id": _to_optional_str(
                    getattr(variant, "ref_id", None) or getattr(variant, "id", None) or index
                ),
                "sequence": str(sequence),
            }


def _find_matches_for_sequence(
    read_seq,
    variant,
    search_sequence,
    max_errors,
    allow_overlaps,
    transform,
    strand,
    direction,
):
    pattern = regex.compile(f"(?e)(?:{regex.escape(search_sequence)}){{e<={max_errors}}}")

    matches = []
    for match in pattern.finditer(read_seq, overlapped=allow_overlaps):
        substitutions, insertions, deletions = match.fuzzy_counts
        errors = (substitutions, insertions, deletions)
        matches.append(
            MatchRecord(
                element_type=variant["element_type"],
                ref_id=variant["ref_id"],
                sequence=variant["sequence"],
                start=match.start(),
                end=match.end(),
                errors=errors,
                penalty=_penalty_from_errors(errors),
                transform=transform,
                strand=strand,
                direction=direction,
            )
        )

    return matches


def find_all_matches(read_seq, reference_variants, max_errors, allow_overlaps=True):
    """
    Возвращает все fuzzy-совпадения reference_variants в read_seq.

    reference_variants может быть строкой, словарем ref_id -> sequence,
    словарем с ключом sequence/seq/pattern/value или списком таких значений.
    Ошибки хранятся явно как (substitutions, insertions, deletions).
    """
    all_matches = []

    for variant in _iter_reference_variants(reference_variants):
        sequence = variant["sequence"]

        all_matches.extend(
            _find_matches_for_sequence(
                read_seq=read_seq,
                variant=variant,
                search_sequence=sequence,
                max_errors=max_errors,
                allow_overlaps=allow_overlaps,
                transform="identity",
                strand="+",
                direction="forward",
            )
        )

        all_matches.extend(
            _find_matches_for_sequence(
                read_seq=read_seq,
                variant=variant,
                search_sequence=reverse_complement(sequence),
                max_errors=max_errors,
                allow_overlaps=allow_overlaps,
                transform="reverse_complement",
                strand="-",
                direction="reverse",
            )
        )

    return sorted(
        all_matches,
        key=lambda record: (
            record.start,
            record.end,
            record.penalty,
            record.element_type or "",
            record.ref_id or "",
            record.transform,
        ),
    )


def find_best_match(read_seq, pattern_seq, max_errors):
    """
    Ищет наилучшее совпадение pattern_seq в read_seq, допуская до max_errors ошибок.
    
    Возвращает:
        tuple (start, end, errors_tuple, is_reverse_complement) 
        - start, end: координаты найденного участка
        - errors_tuple: кортеж (substitutions, insertions, deletions)
        - is_reverse_complement: флаг (True, если найдено на обратной цепи)
    Если совпадений не найдено, возвращает None.
    """
    matches = find_all_matches(read_seq, pattern_seq, max_errors)
    if not matches:
        return None

    best = min(
        matches,
        key=lambda record: (
            record.strand != "+",
            record.penalty,
            record.start,
            record.end,
        ),
    )
    return (best.start, best.end, best.errors, best.strand == "-")


def calculate_element_score(errors, weights_dict):
    """
    Считает скор конкретного найденного элемента на основе словаря весов из конфига.
    """
    subs, ins, dels = errors
    
    mismatch_penalty = weights_dict.get("mismatch", -2)
    insertion_penalty = weights_dict.get("insertion", -3)
    deletion_penalty = weights_dict.get("deletion", -3)
    
    score = (subs * mismatch_penalty) + (ins * insertion_penalty) + (dels * deletion_penalty)
    return score
