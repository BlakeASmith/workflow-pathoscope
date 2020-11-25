from pathlib import Path

TEST_FILES = Path(__file__).parent / "test_files"

BEST_HIT_PATH = TEST_FILES / "best_hit"
RESULTS_PATH = TEST_FILES / "results.json"
EM_PATH = TEST_FILES / "em"
ISOLATES_VTA_PATH = TEST_FILES / "to_isolates.vta"
MATRIX_PATH = TEST_FILES / "ps_matrix"
REF_LENGTHS_PATH = TEST_FILES / "ref_lengths.json"
SAM_PATH = TEST_FILES / "test_al.sam"
SCORES = TEST_FILES / "scores"
TO_SUBTRACTION_PATH = TEST_FILES / "to_subtraction.json"
UNU_PATH = TEST_FILES / "unu"
VTA_PATH = TEST_FILES / "test.vta"
INDEXES_PATH = TEST_FILES /"index"
FASTQ_PATH = TEST_FILES / "test.fq"
HOST_PATH = INDEXES_PATH / "host"

def test_map_default_isolates():
    assert False


def test_generate_isolate_fasta():
    assert False


def test_build_isolate_index():
    assert False


def test_map_isolates():
    assert False


def test_map_subtraction():
    assert False


def test_subtract_mapping():
    assert False


def test_pathoscope():
    assert False


def test_workflow():
    assert False
