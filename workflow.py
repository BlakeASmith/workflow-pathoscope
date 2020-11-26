from pathlib import Path
from typing import Dict, Any, Callable

import pathoscope
from virtool_workflow import fixture
from virtool_workflow import step

# TODO: use `Reads` class that will come in the next virtool_workflow release.
Reads = Any


@fixture
def intermediate() -> Dict[str, Any]:
    return {}


@step
async def map_default_isolates(
        proc: int,
        index_path: Path,
        reads: Reads,
        run_subprocess: Callable,
        intermediate: Dict[str, Any]
):
    command = [
        "bowtie2",
        "-p", str(proc),
        "--no-unal",
        "--local",
        "--score-min", "L,20,1.0",
        "-N", "0",
        "-L", "15",
        "-x", index_path,
        "-U", ",".join(str(path) for path in reads.paths)
    ]

    to_otus = set()

    async def _stdout_handler(line: str):
        nonlocal to_otus

        if line[0] == "#" or line[0] == "@":
            return

        fields = line.split("\t")

        # Bitwise FLAG - 0x4: segment unmapped
        if int(fields[1]) & 0x4 == 4:
            return

        ref_id = fields[2]

        if ref_id == "*":
            return

        # Skip if the p_score does not meet the minimum cutoff.
        if pathoscope.find_sam_align_score(fields) < 0.01:
            return

        to_otus.add(ref_id)

    await run_subprocess(command, stdout_handler=_stdout_handler)

    intermediate["to_otus"] = to_otus
    return "Finished mapping default isolates."


@step
def generate_isolate_fasta(
        temp_analysis_path: Path,
        sequence: Dict
):
    fasta_path = temp_analysis_path/"isolate_index.fa"

    ref_lengths = {}




@step
def build_isolate_index():
    ...


@step
def map_isolates():
    ...


@step
def map_subtraction():
    print("1")
    ...


@step
def subtract_mapping():
    ...


@step
def pathoscope():
    ...


