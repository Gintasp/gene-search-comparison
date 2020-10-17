from typing import Dict

from Bio.SeqRecord import SeqRecord


class FrequencyTable:
    def __init__(self, freq_table: Dict[str, float], seq_record: SeqRecord):
        self.__freq_table = freq_table
        self.__seq_record = seq_record

    def get_table(self) -> Dict[str, float]:
        return self.__freq_table

    def get_record(self) -> SeqRecord:
        return self.__seq_record
