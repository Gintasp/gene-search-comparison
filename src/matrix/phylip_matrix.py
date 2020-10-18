from typing import List

from src.model import FrequencyTable

from config import config


class PhylipMatrix:
    def __init__(self, freq_tables: List[FrequencyTable]):
        self.__freq_tables = freq_tables

    def print(self) -> None:
        self.__print_header()
        for freq_table in self.__freq_tables:
            print(freq_table.get_record().id, end=' ')
            for freq_table_inner in self.__freq_tables:
                distance = self.__calculate_distance(freq_table, freq_table_inner)
                print(f'{format(distance, f".{config.PHYLIP_MATRIX_PRECISION}f")}', end=' ')
            print('')
        print('')

    def __calculate_distance(self, freq_table_1: FrequencyTable, freq_table_2: FrequencyTable) -> float:
        codon_distances = []
        table_1 = freq_table_1.get_table()
        table_2 = freq_table_2.get_table()

        for key in table_1:
            distance = self.__root_mean_square_deviation(table_1[key], table_2[key])
            codon_distances.append(distance)

        return sum(codon_distances)

    @staticmethod
    def __root_mean_square_deviation(x: float, y: float) -> float:
        return (x - y) ** 2

    def __print_header(self) -> None:
        print(f'\n{len(self.__freq_tables)}')
