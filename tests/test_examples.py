import math
import re
import subprocess
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_command(args, cwd=ROOT):
    return subprocess.run(
        args,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
        timeout=30,
    )


def parse_md_output(output):
    rows = []
    for line in output.splitlines():
        parts = line.split()
        if not parts:
            continue
        if len(parts) != 5:
            raise AssertionError(f"Unexpected MD output row: {line!r}")
        step, atom = int(parts[0]), int(parts[1])
        coords = tuple(float(value) for value in parts[2:])
        rows.append((step, atom, coords))
    return rows


def assert_md_output_is_stable(testcase, output):
    rows = parse_md_output(output)
    testcase.assertEqual(len(rows), 150)

    by_step = defaultdict(dict)
    for step, atom, coords in rows:
        testcase.assertTrue(all(math.isfinite(value) for value in coords))
        by_step[step][atom] = coords

    testcase.assertEqual(sorted(by_step), list(range(0, 1000, 100)))

    for atoms in by_step.values():
        testcase.assertEqual(len(atoms), 15)
        for base in range(0, 15, 3):
            oxygen = atoms[base]
            hydrogen1 = atoms[base + 1]
            hydrogen2 = atoms[base + 2]
            testcase.assertAlmostEqual(math.dist(oxygen, hydrogen1), 0.96, delta=0.02)
            testcase.assertAlmostEqual(math.dist(oxygen, hydrogen2), 0.96, delta=0.02)


def assert_docking_output_is_stable(testcase, output, trajectory_file):
    match = re.search(r"Best Docking Score:\s+([0-9.]+)", output)
    testcase.assertIsNotNone(match)
    testcase.assertAlmostEqual(float(match.group(1)), 5.0)
    testcase.assertIn("Best Docking Rotation: 0 degrees", output)
    testcase.assertIn("Best Docking Translation: X=-3", output)

    lines = trajectory_file.read_text().splitlines()
    testcase.assertEqual(len(lines), 8 * 11 * 11 * 12)
    testcase.assertEqual(lines[0], "10")
    testcase.assertTrue(lines[1].startswith("Frame 0:"))


class ExampleRegressionTests(unittest.TestCase):
    def test_perl_md_output_shape_and_bond_lengths(self):
        result = run_command(["perl", "md_simulation/perl/md_simulation_h20.pl"])
        assert_md_output_is_stable(self, result.stdout)

    def test_cpp_md_output_shape_and_bond_lengths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            binary = Path(tmpdir) / "md_simulation_h20"
            run_command(
                [
                    "g++",
                    "-std=c++11",
                    "-O2",
                    "-Wall",
                    "-Wextra",
                    "-pedantic",
                    "-o",
                    str(binary),
                    str(ROOT / "md_simulation/cpp/md_simulation_h20.cpp"),
                ]
            )
            result = run_command([str(binary)])
        assert_md_output_is_stable(self, result.stdout)

    def test_perl_docking_best_score_and_trajectory_shape(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            result = run_command(
                ["perl", str(ROOT / "docking/perl/peptide_docking.pl")],
                cwd=tmp_path,
            )
            assert_docking_output_is_stable(
                self,
                result.stdout,
                tmp_path / "docking_trajectory.xyz",
            )

    def test_cpp_docking_best_score_and_trajectory_shape(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            binary = tmp_path / "peptide_docking"
            run_command(
                [
                    "g++",
                    "-std=c++11",
                    "-O2",
                    "-Wall",
                    "-Wextra",
                    "-pedantic",
                    "-o",
                    str(binary),
                    str(ROOT / "docking/cpp/peptide_docking.cpp"),
                ]
            )
            result = run_command([str(binary)], cwd=tmp_path)
            assert_docking_output_is_stable(
                self,
                result.stdout,
                tmp_path / "docking_trajectory.xyz",
            )


if __name__ == "__main__":
    unittest.main()
