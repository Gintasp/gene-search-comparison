from pathlib import Path


def project_root() -> str:
    return Path(__file__).parent.parent.parent
