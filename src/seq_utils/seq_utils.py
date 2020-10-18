from typing import List

from Bio import SeqIO
from Bio.Seq import Seq

from config import config
from src.io import IO
from src.model import SeqFragmentCollection


class SeqUtils:
    def __init__(self, filename: str):
        self.seq_record = SeqIO.read(filename, config.SEQ_FILE_FORMAT)
        IO.success('Loaded file ' + filename)

    def find_start_stop_fragments(self) -> SeqFragmentCollection:
        start_codons = config.START_CODONS
        stop_codons = config.STOP_CODONS
        start_stop_fragments = []
        nucleotides = self.seq_record.seq

        IO.print('Parsing start - stop codon sequences...')
        for nuc_strand in [nucleotides, nucleotides.reverse_complement()]:
            for frame in range(0, 2):
                frame_codons = self.get_codons(nuc_strand, frame)
                fragment = Seq('')
                is_read_started = False
                for index, codon in enumerate(frame_codons):
                    if codon in start_codons:
                        fragment += codon
                        is_read_started = True
                    elif is_read_started and codon not in stop_codons:
                        fragment += codon
                    elif is_read_started and codon in stop_codons:
                        fragment += codon
                        start_stop_fragments.append(fragment)
                        fragment = ''
                        is_read_started = False

        IO.success('Done')

        return SeqFragmentCollection(
            self.seq_record,
            self.__filter_out_shorter_than(start_stop_fragments, config.MIN_START_STOP_FRAGMENT_LENGTH)
        )

    @staticmethod
    def get_codons(seq: Seq, frame=0) -> List[str]:
        return SeqUtils.__chunk_seq(seq, frame, 3)

    @staticmethod
    def get_dicodons(seq: Seq, frame=0) -> List[str]:
        return SeqUtils.__chunk_seq(seq, frame, 6)

    @staticmethod
    def __chunk_seq(seq: Seq, start: int, length: int) -> List[str]:
        chunks = []

        for index in range(start, len(seq), length):
            chunk = seq[index:index + length]

            if len(chunk) < length:
                continue

            chunks.append(str(chunk))

        return chunks

    @staticmethod
    def __filter_out_shorter_than(fragments: List[Seq], min_length: int) -> List[Seq]:
        filtered = []
        IO.print(f'Filtering out sequences shorter than {min_length} symbols...')
        for fragment in fragments:
            if len(fragment) > min_length:
                filtered.append(fragment)

        IO.success('Done')

        return filtered
