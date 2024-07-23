"""
Unit tests for dnaapler.

Usage: pytest

"""

import glob
import os
import sys

# import
import unittest
from pathlib import Path
from unittest.mock import patch

import click
import pandas as pd
import pytest
from loguru import logger

from src.dnaapler.utils.constants import repo_root
from src.dnaapler.utils.external_tools import ExternalTool
from src.dnaapler.utils.processing import (
    process_blast_output_and_reorient,
    reorient_sequence,
    reorient_sequence_random,
)
from src.dnaapler.utils.util import begin_dnaapler, end_dnaapler
from src.dnaapler.utils.validation import (
    check_evalue,
    validate_choice_autocomplete,
    validate_custom_db_fasta,
    validate_fasta,
)

# import functions


# test data
test_data = Path("tests/test_data")
overall_inputs_test_data = Path("tests/test_data/overall_inputs")
logger.add(lambda _: sys.exit(1), level="ERROR")


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


class TestValidateFasta(unittest.TestCase):
    """Tests of Fasta validation functions"""

    def test_non_fasta_input(self):
        with self.assertRaises(SystemExit):
            non_fasta_file = os.path.join(test_data, "non_fasta.txt")
            validate_fasta(non_fasta_file)

    def test_real_fasta_input(self):
        non_fasta_file = os.path.join(test_data, "nucl_test.fna")
        validate_fasta(non_fasta_file)
        # checks the ctx is the same, no error

    def test_non_fasta_custom_input(self):
        with self.assertRaises(SystemExit):
            nucleotide_fasta_file = os.path.join(test_data, "non_fasta.txt")
            validate_custom_db_fasta(nucleotide_fasta_file)

    def test_non_aa_custom_fasta_input(self):
        with self.assertRaises(SystemExit):
            nucleotide_fasta_file = os.path.join(test_data, "nucl_test.fna")
            validate_custom_db_fasta(nucleotide_fasta_file)


class TestReorientSequence(unittest.TestCase):
    """Tests for reorient_sequence"""

    def test_reorient_sequence_top_hit_no_start_codon(self):
        # Test scenario where the top blast hit doesn't have a start codon
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
        gene = "terL"
        overlapping_orf = True
        reorient_sequence(blast_df, input, out_file, gene, overlapping_orf)

    def test_reorient_sequence_top_hit_with_start_codon(self):
        # Test scenario where the top blast hit does have a start codon
        blast_file = os.path.join(test_data, "NC_007458_dnaapler_blast_output.txt")
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
        input = os.path.join(overall_inputs_test_data, "NC_007458.fasta")
        out_file = os.path.join(
            overall_inputs_test_data, "fake_reoriented_NC_007458.fasta"
        )
        gene = "terL"
        overlapping_orf = False
        reorient_sequence(blast_df, input, out_file, gene, overlapping_orf)


class TestReorientSequenceRandom(unittest.TestCase):
    """Tests for reorient_sequence_random"""

    def test_reorient_sequence_random_bad_strand(self):
        # Test scenario where the strand is outside 1 or -1
        with self.assertRaises(UnboundLocalError):
            input = os.path.join(test_data, "SAOMS1.fasta")
            out_file = os.path.join(test_data, "fake_reoriented.fasta")
            start = 1000
            strand = 24
            reorient_sequence_random(input, out_file, start, strand)


class TestBlastOutput(unittest.TestCase):
    """Tests for process_blast_output_and_reorient"""

    def test_process_blast_output_and_reorient_invalid_blast_file(self):
        # Test scenario where the blast input is invalud
        with self.assertRaises(SystemExit):
            blast_file = pd.DataFrame({"qstart": [1]})
            input = os.path.join(test_data, "SAOMS1.fasta")
            output = os.path.join(test_data, "fake_reoriented.fasta")
            gene = "terL"
            process_blast_output_and_reorient(input, blast_file, output, gene)

    def test_process_blast_output_and_reorient_already_oriented(self):
        # Test scenario where the blast output suggests the contig is already oriented correctly
        blast_file = os.path.join(test_data, "SAOMS1_blast_output_already_oriented.txt")
        input = os.path.join(test_data, "SAOMS1.fasta")
        output = os.path.join(test_data, "fake_reoriented.fasta")
        gene = "terL"
        blast_success = process_blast_output_and_reorient(
            input, blast_file, output, gene
        )
        assert blast_success == True

    def test_process_blast_output_and_reorient_wrong_start_codon(self):
        # Test scenario where the best BLAST hit has no valid start codon
        blast_file = os.path.join(
            test_data, "SAOMS1_blast_output_wrong_start_codon.txt"
        )
        input = os.path.join(test_data, "SAOMS1.fasta")
        output = os.path.join(test_data, "fake_reoriented.fasta")
        gene = "terL"
        blast_success = process_blast_output_and_reorient(
            input, blast_file, output, gene
        )
        assert blast_success == True

    def test_process_blast_output_and_reorient_correct(self):
        # Test scenario where the no BLAST hit begins with 1 (start of gene)
        blast_file = os.path.join(test_data, "SAOMS1_blast_output_correct.txt")
        input = os.path.join(test_data, "SAOMS1.fasta")
        output = os.path.join(test_data, "fake_reoriented.fasta")
        gene = "terL"
        blast_success = process_blast_output_and_reorient(
            input, blast_file, output, gene
        )
        assert blast_success == True

    def test_begin_dnaapler(self):
        # Test begin
        input = os.path.join(test_data, "SAOMS1.fasta")
        threads = str(8)
        gene = "terL"
        outdir = os.path.join(test_data, "bad_dir")
        params = {
        "--input": input,
        "--output": outdir,
        "--threads": threads,
    }
        begin_dnaapler(input, outdir, threads, gene, params)

    def test_end_dnaapler(self):
        time = 2324.0
        end_dnaapler(time)


class TestEValue(unittest.TestCase):
    """Tests of Evalue"""

    def test_evalue_char(self):
        with self.assertRaises(SystemExit):
            evalue = "sfsd"
            check_evalue(evalue)

    def test_evalue_char_mix(self):
        with self.assertRaises(SystemExit):
            evalue = "1e-10t"
            check_evalue(evalue)

    def test_evalue_char_int(self):
        evalue = "5"
        check_evalue(evalue)

    def test_evalue_int(self):
        evalue = 5
        check_evalue(evalue)

    def test_evalue_sci(self):
        evalue = "1e-10"
        check_evalue(evalue)


class TestChoiceAutocomplete(unittest.TestCase):
    """Tests of Choice Autocomplete"""

    def test_evalue_bad_char(self):
        # fake values
        ctx = "1"
        param = "2"
        value = "sfsd"
        with self.assertRaises(click.BadParameter):
            val = validate_choice_autocomplete(ctx, param, value)

    def test_evalue_none(self):
        value = "none"
        ctx = "1"
        param = "2"
        val = validate_choice_autocomplete(ctx, param, value)

    def test_evalue_mys(self):
        value = "mystery"
        ctx = "1"
        param = "2"
        val = validate_choice_autocomplete(ctx, param, value)

    def test_evalue_nearest(self):
        value = "nearest"
        ctx = "1"
        param = "2"
        val = validate_choice_autocomplete(ctx, param, value)


# external tools

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

        external_tool = ExternalTool("tool", "input", "output", "params", logdir)

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
        logsdir = repo_root.parent.parent / "tests/helpers/logs"
        logsdir.mkdir(parents=True, exist_ok=True)
        for file in logsdir.iterdir():
            file.unlink()

        python_script = str(repo_root.parent.parent / "tests/helpers/run_test.py")
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
            assert lines == ["err\n"]


class TestFailExternal(unittest.TestCase):
    """Fail Extenral Tool Test"""

    def test___run_exit(self):
        with self.assertRaises(FileNotFoundError):
            logsdir = repo_root.parent.parent / "tests/helpers/logs"
            logsdir.mkdir(parents=True, exist_ok=True)
            for file in logsdir.iterdir():
                file.unlink()

            python_script = str(repo_root.parent.parent / "tests/helpers/run_test.py")
            external_tool = ExternalTool(
                "break_here",
                "input",
                "output",
                python_script,
                logsdir,
            )

            external_tool.run()
