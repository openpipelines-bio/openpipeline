import importlib
import pytest
from pathlib import Path


def pytest_collect_file(file_path: Path, parent):
    if file_path.name == ".viash_script.sh":
        # Allow file ending in .sh to be imported
        importlib.machinery.SOURCE_SUFFIXES.append(".viash_script.sh")
        return pytest.Module.from_parent(parent, path=file_path)


def pytest_collection_finish(session):
    importlib.machinery.SOURCE_SUFFIXES.remove(".viash_script.sh")
