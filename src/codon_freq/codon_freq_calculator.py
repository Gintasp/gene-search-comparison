import itertools
from typing import Dict, List

from Bio.Seq import Seq

from config import config
from src.model import SeqFragmentCollection
from src.model.frequency_table import FrequencyTable
from src.io import IO
from src.seq_utils import SeqUtils


class CodonFreqCalculator:
    def __init__(self):
        self.__possible_codons = self.__get_all_possible_codons()
        self.__possible_dicodons = self.__get_all_possible_dicodons()

    def get_codon_freq_table(self, seq_fragment_collection: SeqFragmentCollection) -> FrequencyTable:
        freq_table = dict.fromkeys(self.__possible_codons, 0)
        fragments = seq_fragment_collection.get_fragments()

        IO.print('Calculating codon frequencies...')

        for frag in fragments:
            frag_codons = SeqUtils.get_codons(frag)
            for codon in self.__possible_codons:
                freq_table[codon] += frag_codons.count(codon)

        IO.success('Done')

        total_codons = self.__get_total_codons_in_fragments(fragments)
        normalized_table = self.__normalize_freq_table(freq_table, total_codons)

        return FrequencyTable(normalized_table, seq_fragment_collection.get_seq_record())

    def get_dicodon_freq_table(self, seq_fragment_collection: SeqFragmentCollection) -> FrequencyTable:
        freq_table = dict.fromkeys(self.__possible_dicodons, 0)
        fragments = seq_fragment_collection.get_fragments()

        IO.print('Calculating dicodon frequencies...')
        for frag in fragments:
            frag_dicodons = SeqUtils.get_dicodons(frag)
            for dicodon in self.__possible_dicodons:
                freq_table[dicodon] += frag_dicodons.count(dicodon)

        IO.success('Done')

        total_items = self.__get_total_dicodons_in_fragments(fragments)
        normalized_table = self.__normalize_freq_table(freq_table, total_items)

        return FrequencyTable(normalized_table, seq_fragment_collection.get_seq_record())

    @staticmethod
    def __get_all_possible_codons() -> List[str]:
        codons = []

        for comb in itertools.combinations_with_replacement(config.NUCLEOTIDES, 3):
            for perm in itertools.permutations(comb):
                codons.append(''.join(perm))

        return list(dict.fromkeys(codons))

    @staticmethod
    def __get_all_possible_dicodons() -> List[str]:
        dicodons = []

        for comb in itertools.combinations_with_replacement(config.NUCLEOTIDES, 6):
            for perm in itertools.permutations(comb):
                dicodons.append(''.join(perm))

        return list(dict.fromkeys(dicodons))

    @staticmethod
    def __get_total_codons_in_fragments(seq_fragments: List[Seq]) -> int:
        result = 0
        for frag in seq_fragments:
            result += len(SeqUtils.get_codons(frag))

        return result

    @staticmethod
    def __get_total_dicodons_in_fragments(seq_fragments: List[Seq]) -> int:
        result = 0
        for frag in seq_fragments:
            result += len(SeqUtils.get_dicodons(frag))

        return result

    @staticmethod
    def __normalize_freq_table(freq_table: Dict[str, float], total_items: int) -> Dict[str, float]:
        for key in freq_table:
            freq_table[key] = freq_table[key] / total_items

        return freq_table
