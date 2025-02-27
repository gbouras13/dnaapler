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
    process_MMseqs2_output_and_reorient,
    reorient_sequence,
    reorient_sequence_random,
)
from src.dnaapler.utils.util import begin_dnaapler, end_dnaapler, check_duplicate_headers
from src.dnaapler.utils.validation import (
    check_evalue,
    validate_choice_autocomplete,
    validate_custom_db_fasta,
    validate_input,
)

# import functions


# test data
test_data = Path("tests/test_data")
overall_inputs_test_data = Path("tests/test_data/overall_inputs")
logger.add(lambda _: sys.exit(1), level="ERROR")


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


class TestValidateInput(unittest.TestCase):
    """Tests of input validation functions"""

    def test_non_fasta_input(self):
        with self.assertRaises(SystemExit):
            non_fasta_file = os.path.join(test_data, "non_fasta.txt")
            validate_input(non_fasta_file)

    def test_real_fasta_input(self):
        fasta_file = os.path.join(test_data, "nucl_test.fna")
        validate_input(fasta_file)
        # checks the ctx is the same, no error

    def test_circular_gfa_input(self):
        gfa_file = os.path.join(test_data, "nucl_test_circular.gfa")
        validate_input(gfa_file)

    def test_linear_gfa_input(self):
        gfa_file = os.path.join(test_data, "nucl_test_linear.gfa")
        validate_input(gfa_file)

    def test_non_fasta_custom_input(self):
        with self.assertRaises(SystemExit):
            nucleotide_fasta_file = os.path.join(test_data, "non_fasta.txt")
            validate_custom_db_fasta(nucleotide_fasta_file)

    def test_non_aa_custom_fasta_input(self):
        with self.assertRaises(SystemExit):
            nucleotide_fasta_file = os.path.join(test_data, "nucl_test.fna")
            validate_custom_db_fasta(nucleotide_fasta_file)

    def test_no_duplicate_fasta(self):
        fasta_file = os.path.join(test_data, "nucl_test.fna")
        check_duplicate_headers(fasta_file)

    def test_duplicate_fasta(self):
        with self.assertRaises(SystemExit):
            fasta_file = os.path.join(test_data, "duplicate_names.fna")
            check_duplicate_headers(fasta_file)

    def test_no_duplicate_gfa(self):
        gfa_file = os.path.join(test_data, "nucl_test_circular.gfa")
        check_duplicate_headers(gfa_file)

    def test_duplicate_gfa(self):
        with self.assertRaises(SystemExit):
            gfa_file = os.path.join(test_data, "duplicate_names.gfa")
            check_duplicate_headers(gfa_file)


class TestReorientSequence(unittest.TestCase):
    """Tests for reorient_sequence"""

    def test_reorient_sequence_top_hit_no_start_codon(self):
        # Test scenario where the top MMseqs2 hit doesn't have a start codon
        MMseqs2_file = os.path.join(test_data, "SAOMS1_MMseqs2_output_correct.txt")
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
        MMseqs2_df = pd.read_csv(
            MMseqs2_file, delimiter="\t", index_col=False, names=col_list
        )
        input = os.path.join(test_data, "SAOMS1.fasta")
        out_file = os.path.join(test_data, "fake_reoriented.fasta")
        gene = "terL"
        overlapping_orf = True
        reorient_sequence(MMseqs2_df, input, out_file, gene, overlapping_orf)

    def test_reorient_sequence_top_hit_with_start_codon(self):
        # Test scenario where the top MMseqs2 hit does have a start codon
        MMseqs2_file = os.path.join(test_data, "NC_007458_dnaapler_MMseqs2_output.txt")
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
        MMseqs2_df = pd.read_csv(
            MMseqs2_file, delimiter="\t", index_col=False, names=col_list
        )
        input = os.path.join(overall_inputs_test_data, "NC_007458.fasta")
        out_file = os.path.join(
            overall_inputs_test_data, "fake_reoriented_NC_007458.fasta"
        )
        gene = "terL"
        overlapping_orf = False
        reorient_sequence(MMseqs2_df, input, out_file, gene, overlapping_orf)


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
    """Tests for process_MMseqs2_output_and_reorient"""

    def test_process_MMseqs2_output_and_reorient_invalid_MMseqs2_file(self):
        # Test scenario where the MMseqs2 input is invalud
        with self.assertRaises(SystemExit):
            MMseqs2_file = pd.DataFrame({"qstart": [1]})
            input = os.path.join(test_data, "SAOMS1.fasta")
            output = os.path.join(test_data, "fake_reoriented.fasta")
            gene = "terL"
            process_MMseqs2_output_and_reorient(input, MMseqs2_file, output, gene)

    def test_process_MMseqs2_output_and_reorient_already_oriented(self):
        # Test scenario where the MMseqs2 output suggests the contig is already oriented correctly
        MMseqs2_file = os.path.join(
            test_data, "SAOMS1_MMseqs2_output_already_oriented.txt"
        )
        input = os.path.join(test_data, "SAOMS1.fasta")
        output = os.path.join(test_data, "fake_reoriented.fasta")
        gene = "terL"
        MMseqs2_success = process_MMseqs2_output_and_reorient(
            input, MMseqs2_file, output, gene
        )
        assert MMseqs2_success == True

    def test_process_MMseqs2_output_and_reorient_wrong_start_codon(self):
        # Test scenario where the best BLAST hit has no valid start codon
        MMseqs2_file = os.path.join(
            test_data, "SAOMS1_MMseqs2_output_wrong_start_codon.txt"
        )
        input = os.path.join(test_data, "SAOMS1.fasta")
        output = os.path.join(test_data, "fake_reoriented.fasta")
        gene = "terL"
        MMseqs2_success = process_MMseqs2_output_and_reorient(
            input, MMseqs2_file, output, gene
        )
        assert MMseqs2_success == True

    def test_process_MMseqs2_output_and_reorient_correct(self):
        # Test scenario where the no BLAST hit begins with 1 (start of gene)
        MMseqs2_file = os.path.join(test_data, "SAOMS1_MMseqs2_output_correct.txt")
        input = os.path.join(test_data, "SAOMS1.fasta")
        output = os.path.join(test_data, "fake_reoriented.fasta")
        gene = "terL"
        MMseqs2_success = process_MMseqs2_output_and_reorient(
            input, MMseqs2_file, output, gene
        )
        assert MMseqs2_success == True

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
