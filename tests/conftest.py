import pytest
from pathlib import Path

TEST_FILES_PATH = Path("test_files")


@pytest.fixture(scope="session")
def test_files_path():
    return TEST_FILES_PATH


def get_sam_lines():
    with open(TEST_FILES_PATH/"sam_50.sam", "r") as handle:
        return handle.read().split("\n")[0:-1]


@pytest.fixture(params=get_sam_lines(), ids=lambda x: x.split("\t")[0])
def sam_line(request):
    return request.param.split("\t")
