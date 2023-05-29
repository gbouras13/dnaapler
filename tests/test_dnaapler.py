"""
Unit tests for dnaapler.

Usage: pytest

"""

# import
import unittest
import os
import click
from pathlib import Path
import pandas as pd
import pytest


# import functions
from dnaapler import validation
from dnaapler import processing
from dnaapler import external_tools
from dnaapler import dnaapler

# move to folder with mock files. First try Github structure, then try pulled repository structure

test_data = Path("tests/test_data")


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


class TestValidateFasta(unittest.TestCase):
    """Tests of Fasta validation functions"""

    def test_non_fasta_input(self):
        with self.assertRaises(click.exceptions.Exit):
            non_fasta_file = os.path.join(test_data, "non_fasta.txt")
            ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
            validation.validate_fasta(non_fasta_file, ctx)

    def test_real_fasta_input(self):
        non_fasta_file = os.path.join(test_data, "nucl_test.fna")
        ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
        validation.validate_fasta(non_fasta_file, ctx)
        # checks the ctx is the same, no error
        assert ctx == ctx

    def test_non_fasta_custom_input(self):
        with self.assertRaises(click.exceptions.Exit):
            nucleotide_fasta_file = os.path.join(test_data, "non_fasta.txt")
            ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
            validation.validate_custom_db_fasta(nucleotide_fasta_file, ctx)

    def test_non_aa_custom_fasta_input(self):
        with self.assertRaises(click.exceptions.Exit):
            nucleotide_fasta_file = os.path.join(test_data, "nucl_test.fna")
            ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
            validation.validate_custom_db_fasta(nucleotide_fasta_file, ctx)


class TestReorientSequence(unittest.TestCase):
    """Tests for reorient_sequence"""

    def test_reorient_sequence_outside_range(self):
        # Test scenario where the row is outside of range
        with self.assertRaises(KeyError):
            blast_file = os.path.join(test_data, "SAOMS1_blast_output_correct.txt")
            col_list = [
                "qseqid",
                "qlen",
                "sseqid",
                "slen",
                "length",
                "qstart",
                "qend",
                "sstart",
                "send",
                "pident",
                "nident",
                "gaps",
                "mismatch",
                "evalue",
                "bitscore",
                "qseq",
                "sseq",
            ]
            blast_df = pd.read_csv(
                blast_file, delimiter="\t", index_col=False, names=col_list
            )
            input = os.path.join(test_data, "SAOMS1.fasta")
            out_file = os.path.join(test_data, "fake_reoriented.fasta")
            ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
            gene = "terL"
            i = int(3)
            processing.reorient_sequence(blast_df, input, out_file, gene, i)

    def test_reorient_sequence_correct(self):
        # test where it works as expected
        blast_file = os.path.join(test_data, "SAOMS1_blast_output_correct.txt")
        col_list = [
            "qseqid",
            "qlen",
            "sseqid",
            "slen",
            "length",
            "qstart",
            "qend",
            "sstart",
            "send",
            "pident",
            "nident",
            "gaps",
            "mismatch",
            "evalue",
            "bitscore",
            "qseq",
            "sseq",
        ]
        blast_df = pd.read_csv(
            blast_file, delimiter="\t", index_col=False, names=col_list
        )
        input = os.path.join(test_data, "SAOMS1.fasta")
        out_file = os.path.join(test_data, "fake_reoriented.fasta")
        ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
        gene = "terL"
        i = int(1)
        processing.reorient_sequence(blast_df, input, out_file, gene, i)
        assert ctx == ctx


class TestReorientSequenceRandom(unittest.TestCase):
    """Tests for reorient_sequence_random"""

    def test_reorient_sequence_random_bad_strand(self):
        # Test scenario where the strand is outside 1 or -1
        with self.assertRaises(UnboundLocalError):
            input = os.path.join(test_data, "SAOMS1.fasta")
            out_file = os.path.join(test_data, "fake_reoriented.fasta")
            start = 1000
            strand = 24
            processing.reorient_sequence_random(input, out_file, start, strand)


class TestBlastOutput(unittest.TestCase):
    """Tests for process_blast_output_and_reorient"""

    def test_process_blast_output_and_reorient_invalid_blast_file(self):
        # Test scenario where the blast input is invalud
        with self.assertRaises(click.exceptions.Exit):
            blast_file = pd.DataFrame({"qstart": [1]})
            input = os.path.join(test_data, "SAOMS1.fasta")
            output = os.path.join(test_data, "fake_reoriented.fasta")
            ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
            gene = "terL"
            processing.process_blast_output_and_reorient(
                input, blast_file, output, ctx, gene
            )

    def test_process_blast_output_and_reorient_already_oriented(self):
        # Test scenario where the blast output suggests the contig is already oriented correctly
        with self.assertRaises(click.exceptions.Exit):
            blast_file = os.path.join(
                test_data, "SAOMS1_blast_output_already_oriented.txt"
            )
            input = os.path.join(test_data, "SAOMS1.fasta")
            output = os.path.join(test_data, "fake_reoriented.fasta")
            ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
            gene = "terL"
            processing.process_blast_output_and_reorient(
                input, blast_file, output, ctx, gene
            )

    def test_process_blast_output_and_reorient_wrong_start_codon(self):
        # Test scenario where the best BLAST hit has no valid start codon
        with self.assertRaises(click.exceptions.Exit):
            blast_file = os.path.join(
                test_data, "SAOMS1_blast_output_wrong_start_codon.txt"
            )
            input = os.path.join(test_data, "SAOMS1.fasta")
            output = os.path.join(test_data, "fake_reoriented.fasta")
            ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
            gene = "terL"
            processing.process_blast_output_and_reorient(
                input, blast_file, output, ctx, gene
            )

    def test_process_blast_output_and_reorient_no_one(self):
        # Test scenario where the no BLAST hit begins with 1 (start of gene)
        with self.assertRaises(click.exceptions.Exit):
            blast_file = os.path.join(test_data, "SAOMS1_blast_output_no_one.txt")
            input = os.path.join(test_data, "SAOMS1.fasta")
            output = os.path.join(test_data, "fake_reoriented.fasta")
            ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
            gene = "terL"
            processing.process_blast_output_and_reorient(
                input, blast_file, output, ctx, gene
            )

    def test_process_blast_output_and_reorient_correct(self):
        # Test scenario where the no BLAST hit begins with 1 (start of gene)
        blast_file = os.path.join(test_data, "SAOMS1_blast_output_correct.txt")
        input = os.path.join(test_data, "SAOMS1.fasta")
        output = os.path.join(test_data, "fake_reoriented.fasta")
        ctx = click.Context(click.Command("cmd"), obj={"prop": "A Context"})
        gene = "terL"
        processing.process_blast_output_and_reorient(
            input, blast_file, output, ctx, gene
        )
        # checks the ctx is the same, no error
        assert ctx == ctx

    def test_begin_dnaapler(self):
        # Test begin
        input = os.path.join(test_data, "SAOMS1.fasta")
        threads = str(8)
        gene = "terL"
        tmp = 1
        outdir = os.path.join(test_data, "bad_dir")
        dnaapler.begin_dnaapler(input, outdir, threads, gene)
        assert tmp == 1

    def test_end_dnaapler(self):
        # Test scenario where the no BLAST hit begins with 1 (start of gene)
        time = 2324.0
        tmp = 1
        dnaapler.end_dnaapler(time)
        assert tmp == 1


import glob
import sys
from pathlib import Path
from unittest.mock import patch

from dnaapler.constants import repo_root
from dnaapler.external_tools import ExternalTool

### external tools

"""
taken from tbpore
"""


class TestExternalTools:
    @patch.object(
        ExternalTool,
        ExternalTool._build_command.__name__,
        return_value=["mocked", "command", "arg"],
    )
    @patch.object(Path, Path.mkdir.__name__)
    def test___constructor(self, mkdir_mock, build_command_mock):
        logdir = Path("logs")

        external_tool = external_tools.ExternalTool(
            "tool", "input", "output", "params", logdir
        )

        assert external_tool.command == ["mocked", "command", "arg"]
        assert external_tool.command_as_str == "mocked command arg"
        assert (
            external_tool.out_log
            == "logs/tool_c238863b32d18040bbf255fa3bf0dc91e9afa268335b56f51abe3c6d1fd83261.out"
        )
        assert (
            external_tool.err_log
            == "logs/tool_c238863b32d18040bbf255fa3bf0dc91e9afa268335b56f51abe3c6d1fd83261.err"
        )

        build_command_mock.assert_called_once_with("tool", "input", "output", "params")
        mkdir_mock.assert_called_once_with(parents=True, exist_ok=True)

    def test___build_command___simple_command(self):
        expected_escaped_command = ["tool", "param1", "param2", "-o", "out", "-i", "in"]
        actual_escaped_command = ExternalTool._build_command(
            "tool", "-i in", "-o out", "param1 param2"
        )
        assert expected_escaped_command == actual_escaped_command

    def test___build_command___single_quote_escaped(self):
        expected_escaped_command = [
            "tool",
            "params",
            "with",
            "escaped arg",
            "-o",
            "escaped out",
            "-i",
            "escaped in",
        ]
        actual_escaped_command = ExternalTool._build_command(
            "tool", "-i 'escaped in'", "-o 'escaped out'", "params with 'escaped arg'"
        )
        assert expected_escaped_command == actual_escaped_command

    def test___build_command___double_quote_escaped(self):
        expected_escaped_command = [
            "tool",
            "params",
            "with",
            "escaped arg",
            "-o",
            "escaped out",
            "-i",
            "escaped in",
        ]
        actual_escaped_command = ExternalTool._build_command(
            "tool", '-i "escaped in"', '-o "escaped out"', 'params with "escaped arg"'
        )
        assert expected_escaped_command == actual_escaped_command

    def test___run(self):
        logsdir = repo_root / "tests/helpers/logs"
        logsdir.mkdir(parents=True, exist_ok=True)
        for file in logsdir.iterdir():
            file.unlink()

        python_script = str(repo_root / "tests/helpers/run_test.py")
        external_tool = ExternalTool(
            sys.executable,
            "input",
            "output",
            python_script,
            logsdir,
        )

        external_tool.run()

        out_file = glob.glob(f"{logsdir}/*.out")[0]
        with open(out_file) as out_file_fh:
            lines = out_file_fh.readlines()
            assert lines == ["out\n"]

        err_file = glob.glob(f"{logsdir}/*.err")[0]
        with open(err_file) as err_file_fh:
            lines = err_file_fh.readlines()
            assert lines == [
                "err\n",
                f"Command line: {sys.executable} {python_script} output input\n",
            ]
