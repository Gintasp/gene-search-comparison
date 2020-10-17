from datetime import datetime

from config import config


class IO:
    @staticmethod
    def success(content) -> None:
        print(f'[OK]: {content}')

    @staticmethod
    def failure(content) -> None:
        print(f'[FAILED]: {content}')

    @staticmethod
    def print(content) -> None:
        print(f'[{datetime.now().strftime(config.LOGGER_TIME_FORMAT)}]: {content}')
