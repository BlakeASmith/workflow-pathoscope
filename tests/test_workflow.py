import json
import os
import shutil
import workflow

from pathlib import Path
import pytest

from virtool_workflow import WorkflowFixtureScope
from virtool_workflow_runtime.db.db import VirtoolDatabase
from virtool_workflow_runtime.config.configuration import db_name, db_connection_string

TEST_FILES_PATH = Path("test_files")
BEST_HIT_PATH = TEST_FILES_PATH / "best_hit"
RESULTS_PATH = TEST_FILES_PATH/"results.json"
EM_PATH = TEST_FILES_PATH/"em"
ISOLATES_VTA_PATH = TEST_FILES_PATH/"to_isolates.vta"
MATRIX_PATH = TEST_FILES_PATH/"ps_matrix"
REF_LENGTHS_PATH = TEST_FILES_PATH/"ref_lengths.json"
SAM_PATH = TEST_FILES_PATH/"test_al.sam"
SCORES = TEST_FILES_PATH/"scores"
TO_SUBTRACTION_PATH = TEST_FILES_PATH/"to_subtraction.json"
UNU_PATH = TEST_FILES_PATH/"unu"
VTA_PATH = TEST_FILES_PATH/"test.vta"
INDEXES_PATH = TEST_FILES_PATH/"index"
FASTQ_PATH = TEST_FILES_PATH/"test.fq"
HOST_PATH = str(INDEXES_PATH/"host")


@pytest.fixture(scope="session")
def scope():
    with WorkflowFixtureScope() as _scope:
        yield _scope


@pytest.fixture(scope="session")
def otu_resource():
    map_dict, otus = {}, {}

    with VTA_PATH.open("r") as handle:
        for line in handle:
            ref_id = line.split(",")[1]

            otu_id = "otu_{}".format(ref_id)

            map_dict[ref_id] = otu_id

            otus[otu_id] = {
                "otu": otu_id,
                "version": 2
            }

    return map_dict, otus


@pytest.fixture(scope="session", autouse=True)
def db(otu_resource):
    db = VirtoolDatabase(db_name(), db_connection_string())

    sequence_otu_map, _ = otu_resource

    await db["analyses"].insert_one({
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

    await db["jobs"].insert_one({
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

    await db["indexes"].insert_one({
        "_id": "index3",
        "manifest": {
            "foobar": 10,
            "reo": 5,
            "baz": 6
        },
        "sequence_otu_map": sequence_otu_map
    })

    await db["indexes"].insert_one({
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

    return db


@pytest.fixture
async def sample_path(scope) -> Path:
    return await scope.get_or_instantiate("sample_path")


@pytest.fixture
async def analysis_work_path(scope) -> Path:
    return await scope.get_or_instantiate("analysis_work_path")


async def test_map_default_isolates(scope: WorkflowFixtureScope, sample_path: Path):
    shutil.copyfile(FASTQ_PATH, sample_path/"reads_1.fq")

    map_default_isolates = await scope.bind(workflow.map_default_isolates)
    await map_default_isolates()

    assert sorted(scope["otu_ids"]) == sorted([
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


async def test_map_isolates(snapshot, scope, sample_path: Path):
    shutil.copyfile(FASTQ_PATH, sample_path/"reads_1.fq")

    analysis_work_path = await scope.get_or_instantiate("analysis_work_path")

    for filename in os.listdir(INDEXES_PATH):
        if "reference" in filename:
            shutil.copyfile(
                Path(INDEXES_PATH, filename),
                Path(analysis_work_path, filename.replace("reference", "isolates"))
            )

    scope["number_of_processes"] = 2

    map_isolates = await scope.bind(workflow.map_isolates)
    await map_isolates()

    vta_path = analysis_work_path/"to_isolates.vta"

    with vta_path.open("r") as f:
        data = sorted([line.rstrip() for line in f])
        snapshot.assert_match(data, "isolates")


async def test_map_subtraction(snapshot, scope):
    scope["number_of_processes"] = True
    scope["subtraction_path"] = HOST_PATH

    analysis_work_path: Path = await scope.get_or_instantiate("analysis_work_path")

    shutil.copyfile(FASTQ_PATH, analysis_work_path/"mapped.fastq")

    map_subtraction = await scope.bind(workflow.map_subtraction)
    map_subtraction()

    sorted_lines = sorted(scope["subtraction_ids"])

    snapshot.assert_match(sorted_lines, "subtraction")


async def test_subtract_mapping(scope, analysis_work_path: Path):
    with TO_SUBTRACTION_PATH.open("r") as handle:
        scope["subtraction_ids"] = json.load(handle)

    shutil.copyfile(VTA_PATH, analysis_work_path/"to_isolates.vta")

    subtract_mapping = await scope.bind(workflow.subtract_mapping)
    await subtract_mapping()

    assert scope["results"]["subtracted_count"] == 4


async def test_pathoscope(snapshot, analysis_work_path: Path, scope):
    with open(REF_LENGTHS_PATH, "r") as handle:
        scope["intermediate"]["ref_lengths"] = json.load(handle)

    shutil.copyfile(
        VTA_PATH,
        analysis_work_path / "to_isolates.vta"
    )

    scope["sequence_otu_map"] = {
        "NC_016509": "foobar",
        "NC_001948": "foobar",
        "13TF149_Reovirus_TF1_Seg06": "reo",
        "13TF149_Reovirus_TF1_Seg03": "reo",
        "13TF149_Reovirus_TF1_Seg07": "reo",
        "13TF149_Reovirus_TF1_Seg02": "reo",
        "13TF149_Reovirus_TF1_Seg08": "reo",
        "13TF149_Reovirus_TF1_Seg11": "reo",
        "13TF149_Reovirus_TF1_Seg04": "reo",
        "NC_004667": "foobar",
        "NC_003347": "foobar",
        "NC_003615": "foobar",
        "NC_003689": "foobar",
        "NC_011552": "foobar",
        "KX109927": "baz",
        "NC_008039": "foobar",
        "NC_015782": "foobar",
        "NC_016416": "foobar",
        "NC_003623": "foobar",
        "NC_008038": "foobar",
        "NC_001836": "foobar",
        "JQ080272": "baz",
        "NC_017938": "foobar",
        "NC_008037": "foobar",
        "NC_007448": "foobar"
    }

    pathoscope = await scope.bind(workflow.pathoscope)
    await pathoscope()

    with (analysis_work_path/"reassigned.vta").open("r") as f:
        data = sorted([line.rstrip() for line in f])
        snapshot.assert_match(data)

    with (analysis_work_path/"report.tsv").open("r") as f:
        data = sorted([line.rstrip() for line in f])
        snapshot.assert_match(data)

    snapshot.assert_match(scope["results"])

