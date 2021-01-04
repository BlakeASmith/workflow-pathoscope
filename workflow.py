import aiofiles
from pathlib import Path
from typing import Dict, Any, Callable

import pathoscope
import virtool_core.history.db
from virtool_workflow import fixture, step
from virtool_workflow.analysis.reads.reads import Reads
from virtool_workflow.execution.run_subprocess import RunSubprocess




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

        segment_unmapped = int(fields[1]) & 0x4 == 4
        if segment_unmapped:
            return

        ref_id = fields[2]

        if ref_id == "*":
            return

        # Skip if the p_score does not meet the minimum cutoff.
        p_score = pathoscope.find_sam_align_score(fields)
        if p_score < 0.01:
            return

        to_otus.add(ref_id)

    await run_subprocess(command, stdout_handler=_stdout_handler)

    intermediate["to_otus"] = to_otus
    return "Finished mapping default isolates."


@step
def generate_isolate_fasta(
        temp_analysis_path: Path,
        sequence_otu_map: Dict,
        intermediate: Dict,
        database: Dict,
        data_path: Path,
):
    fasta_path = temp_analysis_path/"isolate_index.fa"

    # The ids of OTUs whose default sequences had mappings.
    otu_ids = {sequence_otu_map[sequence_id] for sequence_id in intermediate["to_otus"]}

    # Get the database documents for the sequences
    async with aiofiles.open(fasta_path, "w") as f:
        # Iterate through each otu id referenced by the hit sequence ids.
        for otu_id in otu_ids:
            # TODO: get manifest as fixture
            otu_version = ...
            #otu_version = job.params["manifest"][otu_id]

            # TODO: create utility for `patch_do_version`
            _, patched, _ = await virtool_core.history.db.patch_to_version(
                database,
                str(data_path),
                otu_id,
                otu_version
            )

            ref_lengths = {}
            for isolate in patched["isolates"]:
                for sequence in isolate["sequences"]:
                    await f.write(f">{sequence['_id']}\n{sequence['sequence']}\n")
                    ref_lengths[sequence["_id"]] = len(sequence["sequence"])

    del intermediate["to_otus"]

    intermediate["ref_lengths"] = ref_lengths


@step
def build_isolate_index(run_subprocess: RunSubprocess,
                        number_of_processes: int,
                        temp_analysis_path: Path):
    command = [
        "bowtie2-build",
        "--threads", str(number_of_processes),
        temp_analysis_path/"isolate_index.fa",
        temp_analysis_path/"isolates",
    ]

    await run_subprocess(command)


@step
def map_isolates(
        number_of_processes: int,
        temp_analysis_path: Path,
        reads: Reads,
        run_subprocess: RunSubprocess
):
    command = [
        "bowtie2",
        "-p", str(number_of_processes - 1),
        "--no-unal",
        "--local",
        "--score-min", "L,20,1.0",
        "-N", "0",
        "-L", "15",
        "-k", "100",
        "--al", temp_analysis_path/"mapped.fastq",
        "-x", temp_analysis_path/"isolates",
        "-U", ",".join(str(path) for path in reads.paths)
    ]

    async with aiofiles.open(temp_analysis_path/"to_isolates.vta", "w") as f:
        async def stdout_handler(line, p_score_cutoff=0.01):
            line = line.decode()

            if line[0] == "@" or line == "#":
                return

            fields = line.split("\t")

            # Bitwise FLAG - 0x4 : segment unmapped
            if int(fields[1]) & 0x4 == 4:
                return

            ref_id = fields[2]

            if ref_id == "*":
                return

            p_score = pathoscope.find_sam_align_score(fields)

            # Skip if the p_score does not meet the minimum cutoff.
            if p_score < p_score_cutoff:
                return

            read_id = fields[0]
            read_pos = fields[3]
            read_length = len(fields[9])

            read_info = [
                read_id,
                ref_id,
                read_pos,
                str(read_length),
                str(p_score)
            ]

            await f.write(f"{','.join(read_info)}\n")

        await run_subprocess(command, stdout_handler=stdout_handler)



@step
def map_subtraction():
    ...


@step
def subtract_mapping():
    ...


@step
def pathoscope():
    ...


