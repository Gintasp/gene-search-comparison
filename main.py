from os import listdir
from time import time

from src.codon_freq import CodonFreqCalculator
from src.io import IO
from src.matrix import PhylipMatrix
from src.seq_utils import SeqUtils
from src.utils import project_root


def get_file_paths():
    paths = []

    dna_files = listdir(f'{project_root()}/dna-data')
    dna_files.sort()

    for filename in dna_files:
        paths.append(f'{project_root()}/dna-data/{filename}')

    return paths


if __name__ == '__main__':
    start_time = time()
    freqCalculator = CodonFreqCalculator()
    codon_freq_tables = []
    dicodon_freq_tables = []

    for file in get_file_paths():
        seq_fragment_collection = (SeqUtils(file)).find_start_stop_fragments()
        codon_freq_tables.append(freqCalculator.get_codon_freq_table(seq_fragment_collection))
        dicodon_freq_tables.append(freqCalculator.get_dicodon_freq_table(seq_fragment_collection))
        print('')

    IO.print('Generating Phylip matrix for codon frequencies...')
    codon_matrix = PhylipMatrix(codon_freq_tables)
    codon_matrix.print()

    IO.print('Generating Phylip matrix for dicodon frequencies...')
    dicodon_matrix = PhylipMatrix(dicodon_freq_tables)
    dicodon_matrix.print()

    end_time = time()

    print(f'\nAll done, took {round(end_time - start_time, 2)}s')
