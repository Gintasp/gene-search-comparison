from typing import List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class SeqFragmentCollection:
    def __init__(self, seq_record: SeqRecord, seq_fragments: List[Seq]):
        self.__seq_record = seq_record
        self.__seq_fragments = seq_fragments

    def get_fragments(self) -> List[Seq]:
        return self.__seq_fragments

    def get_seq_record(self) -> SeqRecord:
        return self.__seq_record
