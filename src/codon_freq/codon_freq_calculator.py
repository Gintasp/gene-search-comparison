import itertools
from typing import Dict, List

from Bio.Seq import Seq

from config import config
from src.model import SeqFragmentCollection
from src.model.frequency_table import FrequencyTable
from src.io import IO
from src.seq_utils import SeqUtils


class CodonFreqCalculator:
    def get_codon_freq_table(self, seq_fragment_collection: SeqFragmentCollection) -> FrequencyTable:
        freq_table = {}
        fragments = seq_fragment_collection.get_fragments()

        IO.print('Calculating codon frequencies...')
        for codon in self.__get_all_possible_codons():
            for frag in fragments:
                for frag_cod in SeqUtils.get_codons(frag):
                    if codon == frag_cod:
                        if codon in freq_table:
                            freq_table[codon] = freq_table[codon] + 1
                        else:
                            freq_table[codon] = 1

        IO.success('Done')

        normalized_table = self.__normalize_freq_table(freq_table, fragments)

        return FrequencyTable(normalized_table, seq_fragment_collection.get_seq_record())

    @staticmethod
    def __get_all_possible_codons() -> List[str]:
        codons = []

        for comb in itertools.combinations_with_replacement(config.NUCLEOTIDES, 3):
            for perm in itertools.permutations(comb):
                codons.append(''.join(perm))

        return list(dict.fromkeys(codons))

    @staticmethod
    def __get_total_codons_in_fragments(seq_fragments: List[Seq]) -> int:
        result = 0
        for frag in seq_fragments:
            result += len(SeqUtils.get_codons(frag))

        return result

    def __normalize_freq_table(self, freq_table: Dict[str, float], seq_fragments: List[Seq]) -> Dict[str, float]:
        total_codons = self.__get_total_codons_in_fragments(seq_fragments)

        for key in freq_table:
            freq_table[key] = freq_table[key] / total_codons

        return freq_table
