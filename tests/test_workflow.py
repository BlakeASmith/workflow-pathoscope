import pytest
from pathlib import Path
from virtool_workflow.fixtures.scope import WorkflowFixtureScope
from virtool_workflow_runtime.db import VirtoolDatabase
from virtool_workflow_runtime.config.configuration import db_connection_string, db_name

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


@pytest.fixture(scope="session")
def otu_resource():
    maps, otus = {}, {}

    with open(VTA_PATH, "r") as vta:
        for line in vta:
            ref_id = line.split(",")[1]
            otu_id = f"otu_{ref_id}"

            maps[ref_id] = otu_id

            otus[otu_id] = {
                "otu": otu_id,
                "version": 2
            }

    return maps, otus


@pytest.fixture(scope="module")
def fixture_scope(otu_resource):
    db = VirtoolDatabase(db_name(), db_connection_string())
    jobs, samples, analyses, indexes = db["jobs"], db["samples"], db["analyses"], db["indexes"]

    await analyses.insert_one({
        "_id": "baz",
        "workflow": "pathoscope_bowtie",
        "ready": False,
        "sample": {
            "id": "foobar"
        },
        "subtraction": {
            "id": "Prunus persica"
        }
    })

    await jobs.insert_one({
        "_id": "foobar",
        "task": "pathoscope_bowtie",
        "args": {
            "sample_id": "foobar",
            "analysis_id": "baz",
            "ref_id": "original",
            "index_id": "index3"
        },
        "proc": 2,
        "mem": 8
    })

    await indexes.insert_one({
        "_id": "index3",
        "manifest": {
            "foobar": 10,
            "reo": 5,
            "baz": 6
        },
        "sequence_otu_map": otu_resource
    })

    await samples.insert_one({
        "_id": "foobar",
        "paired": False,
        "library_type": "normal",
        "quality": {
            "count": 1337,
            "length": [78, 101]
        },
        "subtraction": {
            "id": "Arabidopsis thaliana"
        }
    })

    scope = WorkflowFixtureScope()

    scope["job_id"] = "foobar"
    scope["manifest"] = {
        "foobar": 10,
        "reo": 5,
        "baz": 6,
    }

    return scope



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
