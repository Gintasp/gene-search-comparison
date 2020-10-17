from os import listdir
from time import time

from src.codon_freq import CodonFreqCalculator
from src.matrix import PhylipMatrix
from src.seq_utils import SeqUtils
from src.utils import project_root


def get_filenames():
    files = []

    for filename in listdir(f'{project_root()}/dna-data'):
        files.append(f'{project_root()}/dna-data/{filename}')

    return files


if __name__ == '__main__':
    start_time = time()
    freqCalculator = CodonFreqCalculator()
    freq_tables = []

    for file in get_filenames():
        seq_fragment_collection = (SeqUtils(file)).find_start_stop_fragments()
        freq_tables.append(freqCalculator.get_codon_freq_table(seq_fragment_collection))
        print('')

    matrix = PhylipMatrix(freq_tables)
    matrix.print()

    end_time = time()

    print(f'\nAll done, took {round(end_time - start_time, 2)}s')
