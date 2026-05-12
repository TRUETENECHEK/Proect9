from __future__ import annotations

import os
import yaml
from collections.abc import Iterable, Mapping
from typing import Dict, Generator, Tuple, List, Any, TYPE_CHECKING

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None

if TYPE_CHECKING:
    from reference_builder import ReferenceBundle


def _require_seqio():
    if SeqIO is None:
        raise ImportError("Biopython is required for FASTA/FASTQ parsing")
    return SeqIO


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Читает конфигурационный файл в формате YAML и возвращает его содержимое в виде словаря.

    Args:
        config_path (str): Путь к YAML файлу конфигурации.

    Returns:
        Dict[str, Any]: Словарь с параметрами конфигурации.
        
    Raises:
        FileNotFoundError: Если файл конфигурации не найден.
        yaml.YAMLError: Если произошла ошибка при парсинге YAML.
    """
    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    return config if config is not None else {}

def load_fasta_dict(fasta_path: str) -> Dict[str, str]:
    """
    Читает FASTA файл и возвращает словарь, где ключи — это ID последовательностей (например, имена баркодов),
    а значения — сами последовательности в виде строк.

    Args:
        fasta_path (str): Путь к файлу в формате FASTA.

    Returns:
        Dict[str, str]: Словарь вида {"ID_последовательности": "ACGT..."}.
        
    Raises:
        FileNotFoundError: Если FASTA файл не найден.
    """
    fasta_dict: Dict[str, str] = {}
    seqio = _require_seqio()
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for record in seqio.parse(f, "fasta"):
            fasta_dict[record.id] = str(record.seq)
    return fasta_dict


def write_primers_fasta(
    primers: Mapping[str, str] | Iterable[Mapping[str, Any]],
    fasta_path: str = "Read_file/primer.fasta",
) -> None:
    """
    Записывает найденные праймеры в FASTA, всегда перезаписывая файл.

    Args:
        primers: Либо словарь {record_id: sequence}, либо iterable записей с
            ключами id/name и sequence/seq.
        fasta_path: Путь к FASTA файлу для записи.
    """
    output_dir = os.path.dirname(fasta_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    if isinstance(primers, Mapping):
        records = [
            {"id": str(record_id), "sequence": str(sequence)}
            for record_id, sequence in primers.items()
        ]
    else:
        records = list(primers)

    with open(fasta_path, "w", encoding="utf-8") as handle:
        for index, record in enumerate(records, start=1):
            record_id = str(record.get("id") or record.get("name") or f"primer_{index}")
            sequence = str(record.get("sequence") or record.get("seq") or "").upper()
            if not sequence:
                continue
            handle.write(f">{record_id}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(f"{sequence[offset:offset + 80]}\n")


def load_reference_bundle(config_path: str = "config.yaml") -> "ReferenceBundle":
    """
    Загружает YAML-конфиг и FASTA-референсы, затем строит предобработанный набор
    adapter/barcode_core с вариантами ориентации.

    Старая загрузка YAML/FASTA остается доступной через load_config и load_fasta_dict.
    """
    from reference_builder import build_reference_bundle_from_config

    return build_reference_bundle_from_config(config_path)


def read_fastq_generator(fastq_path: str) -> Generator[Tuple[str, str, List[int]], None, None]:
    """
    Генератор для ленивого чтения FASTQ файла. Выдает по одному риду за раз, что позволяет
    обрабатывать большие файлы без переполнения оперативной памяти.

    Args:
        fastq_path (str): Путь к файлу в формате FASTQ.

    Yields:
        Tuple[str, str, List[int]]: Кортеж, содержащий:
            - ID рида (str)
            - Нуклеотидную последовательность (str)
            - Качество (Phred quality scores) в виде списка целых чисел (List[int])
            
    Raises:
        FileNotFoundError: Если FASTQ файл не найден.
    """
    seqio = _require_seqio()
    with open(fastq_path, 'r', encoding='utf-8') as f:
        for record in seqio.parse(f, "fastq"):
            # Извлекаем ID, последовательность в виде строки и оценки качества Phred
            read_id = record.id
            sequence = str(record.seq)
            phred_quality = record.letter_annotations["phred_quality"]
            
            yield (read_id, sequence, phred_quality)


def load_mapping_excel(file_path: str) -> Dict[str, str]:
    import pandas as pd
    if not os.path.exists(file_path): return {}
    try:
        df = pd.read_excel(file_path)
        mapping = {}
        for _, row in df.iterrows():
            bc = str(row.iloc[0]).strip()
            samp = str(row.iloc[1]).strip()
            mapping[bc] = samp
        return mapping
    except Exception as e:
        print(f"Warning: Could not load mapping from {file_path}: {e}")
        return {}
