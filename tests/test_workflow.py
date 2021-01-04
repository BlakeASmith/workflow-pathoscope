import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

import workflow
from virtool_workflow.fixtures.scope import WorkflowFixtureScope
from virtool_workflow_runtime.config.configuration import db_connection_string, db_name
from virtool_workflow_runtime.db import VirtoolDatabase
from virtool_workflow.storage.paths import temp_path
from virtool_workflow.analysis.analysis_info import AnalysisInfo

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


@pytest.fixture()
async def fixture_scope(otu_resource):
    db = VirtoolDatabase(db_name(), db_connection_string())
    jobs, samples, analyses, indexes = db["jobs"], db["samples"], db["analyses"], db["indexes"]

    scope = WorkflowFixtureScope()

    scope["analysis_document"] = {
        "_id": "baz",
        "workflow": "pathoscope_bowtie",
        "ready": False,
        "sample": {
            "id": "foobar"
        },
        "subtraction": {
            "id": "Prunus persica"
        }
    }

    scope["job_id"] = "foobar"
    scope["job_document"] = {
        "_id": scope["job_id"],
        "task": "pathoscope_bowtie",
        "args": {
            "sample_id": "foobar",
            "analysis_id": "baz",
            "ref_id": "original",
            "index_id": "index3"
        },
        "proc": 2,
        "mem": 8
    }

    scope["analysis_info"] = AnalysisInfo(scope["job_document"]["sample_id"])

    scope["index_document"] = {
        "_id": "index3",
        "manifest": {
            "foobar": 10,
            "reo": 5,
            "baz": 6
        },
        "sequence_otu_map": otu_resource
    }

    scope["sample_document"] = {
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
    }

    scope["manifest"] = {
        "foobar": 10,
        "reo": 5,
        "baz": 6,
    }

    await scope.get_or_instantiate("temp_path")
    scope["sample_path"] = scope["temp_path"]/"samples/foobar"
    scope["sample_path"].mkdir(parents=True)

    await scope.get_or_instantiate("analysis_args")

    return scope


async def test_map_default_isolates(fixture_scope: WorkflowFixtureScope):
    shutil.copyfile(FASTQ_PATH, fixture_scope["sample_path"]/"reads_1.fq")

    fixture_scope["reads"] = SimpleNamespace(paths=fixture_scope["analysis_args"].read_paths)

    bound = await fixture_scope.bind(workflow.map_default_isolates)
    await bound()

    assert sorted(fixture_scope["intermediate"]["to_otus"]) == sorted([
        "NC_013110",
        "NC_017938",
        "NC_006057",
        "NC_007448",
        "JQ080272",
        "NC_001836",
        "NC_003347",
        "NC_016509",
        "NC_017939",
        "NC_006056",
        "NC_003623",
        "KX109927",
        "NC_016416",
        "NC_001948",
        "NC_021148",
        "NC_003615",
        "NC_004006"
    ])


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
