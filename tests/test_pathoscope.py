import os
import pytest
import shutil
import pickle
import filecmp
import difflib
from pathlib import Path

import pathoscope


@pytest.fixture(scope="session")
def expected_em(test_files_path):
    with open(test_files_path/"em", "rb") as handle:
        em_dict = pickle.load(handle)

    return em_dict


@pytest.fixture(scope="session")
def expected_scores(test_files_path):
    with open(test_files_path/"scores", "rb") as handle:
        scores = pickle.load(handle)

    return scores


def test_find_sam_align_score(sam_line, expected_scores):
    assert pathoscope.find_sam_align_score(sam_line) == expected_scores["".join(sam_line)]


def test_build_matrix(tmpdir, test_files_path):
    shutil.copy(test_files_path/"test.vta", str(tmpdir))
    vta_path = os.path.join(str(tmpdir), "test.vta")

    with open(test_files_path/"ps_matrix", "rb") as handle:
        expected = pickle.load(handle)

    actual = pathoscope.build_matrix(vta_path, 0.01)

    assert expected[0] == actual[0]
    assert expected[1] == actual[1]

    assert sorted(expected[2]) == sorted(actual[2])
    assert sorted(expected[3]) == sorted(actual[3])

@pytest.mark.parametrize("theta_prior", [0, 1e-5])
@pytest.mark.parametrize("pi_prior", [0, 1e-5])
@pytest.mark.parametrize("epsilon", [1e-6, 1e-7, 1e-8])
@pytest.mark.parametrize("max_iter", [5, 10, 20, 30])
def test_em(tmpdir, theta_prior, pi_prior, epsilon, max_iter, expected_em, test_files_path):
    shutil.copy(test_files_path/"test.vta", str(tmpdir))
    vta_path = os.path.join(str(tmpdir), "test.vta")

    u, nu, refs, _ = pathoscope.build_matrix(vta_path, 0.01)

    result = pathoscope.em(u, nu, refs, max_iter, epsilon, pi_prior, theta_prior)

    file_string = "_".join([str(i) for i in ["em", theta_prior, pi_prior, epsilon, max_iter]])

    for i in [0, 1, 2]:
        assert sorted(result[i]) == sorted(expected_em[file_string][i])

    assert result[3] == expected_em[file_string][3]


def test_compute_best_hit(test_files_path):
    """
    Test that :meth:`compute_best_hit` gives the expected result given some input data.

    """
    with open(test_files_path/"ps_matrix", "rb") as handle:
        matrix_tuple = pickle.load(handle)

    with open(test_files_path/"best_hit", "rb") as handle:
        assert pickle.load(handle) == pathoscope.compute_best_hit(*matrix_tuple)


def test_rewrite_align(tmpdir, test_files_path):
    with open(test_files_path/"unu", "rb") as f:
        u, nu = pickle.load(f)

    shutil.copy(test_files_path/"test.vta", str(tmpdir))
    vta_path = os.path.join(str(tmpdir), "test.vta")

    rewrite_path = os.path.join(str(tmpdir), "rewrite.vta")

    pathoscope.rewrite_align(u, nu, vta_path, 0.01, rewrite_path)

    assert filecmp.cmp(test_files_path/"updated.vta", rewrite_path)
    assert not filecmp.cmp(vta_path, rewrite_path)


def test_calculate_coverage(tmpdir, test_files_path):
    ref_lengths = dict()

    shutil.copy(test_files_path/"test.vta", str(tmpdir))
    vta_path = os.path.join(str(tmpdir), "test.vta")

    with open(test_files_path/"test_al.sam", "r") as handle:
        for line in handle:
            if line[0:3] == "@SQ":
                ref_id = None
                length = None

                for field in line.rstrip().split("\t"):
                    if "SN:" in field:
                        ref_id = field.split(":")[1]
                    if "LN:" in field:
                        length = int(field.split(":")[1])

                assert ref_id and length
                assert ref_id not in ref_lengths

                ref_lengths[ref_id] = length

    pathoscope.calculate_coverage(vta_path, ref_lengths)


def test_write_report(tmpdir, test_files_path):
    shutil.copy(test_files_path/"test.vta", str(tmpdir))
    vta_path = os.path.join(str(tmpdir), "test.vta")

    u, nu, refs, reads = pathoscope.build_matrix(vta_path, 0.01)

    with open(test_files_path/"best_hit", "rb") as handle:
        best_hit_initial_reads, best_hit_initial, level_1_initial, level_2_initial = pickle.load(handle)

    init_pi, pi, _, nu = pathoscope.em(u, nu, refs, 30, 1e-7, 0, 0)

    best_hit_final_reads, best_hit_final, level_1_final, level_2_final = pathoscope.compute_best_hit(
        u,
        nu,
        refs,
        reads
    )

    report_path = os.path.join(str(tmpdir), "report.tsv")

    pathoscope.write_report(
        report_path,
        pi,
        refs,
        len(reads),
        init_pi,
        best_hit_initial,
        best_hit_initial_reads,
        best_hit_final,
        best_hit_final_reads,
        level_1_initial,
        level_2_initial,
        level_1_final,
        level_2_final
    )

    with open(report_path, "r") as report:
        with open(test_files_path/"report.tsv") as test_report:
            diff = difflib.unified_diff(report.readlines(), test_report.readlines())
            assert len(list(diff)) == 0



