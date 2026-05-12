import pandas as pd
import subprocess

def load_sample_mapping(excel_path):
    """
    Загружает маппинг баркодов на SampleID из Excel файла.
    
    Ожидается таблица с колонками: Barcode1, Barcode2, SampleID.
    
    Возвращает словарь:
        {(Barcode1, Barcode2): SampleID}
    """
    df = pd.read_excel(excel_path)
    mapping = {}
    
    for _, row in df.iterrows():
        bc1 = str(row['Barcode1']).strip()
        bc2 = str(row['Barcode2']).strip()
        sample_id = str(row['SampleID']).strip()
        mapping[(bc1, bc2)] = sample_id
        
    return mapping


def run_cutadapt_preprocessing(input_fastq, output_fastq, left_adapter_seq, right_adapter_seq):
    """
    Запускает cutadapt для предварительной обрезки внешних общих адаптеров.
    Если cutadapt не установлен, возвращает исходный FASTQ без ошибки.
    """
    try:
        subprocess.run(
            [
                "cutadapt",
                "-a", left_adapter_seq,
                "-A", right_adapter_seq,
                "-o", output_fastq,
                input_fastq,
            ],
            check=True,
        )
    except FileNotFoundError:
        print("Cutadapt не найден, пропускаем этап")
        return input_fastq

    return output_fastq
