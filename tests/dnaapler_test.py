'''
Unit tests for dnaapler.

Usage: pytest

'''

# import
import unittest
import os
import logging
import click
from pathlib import Path


# import functions
from dnaapler import validation

# move to folder with mock files. First try Github structure, then try pulled repository structure

test_data = Path("tests/test_data")


class TestValidateFasta(unittest.TestCase):
    """ Test for the function that examines the cutoffs given for core and low-frequency genes"""

    def test_non_fasta_input(self):
        with self.assertRaises(click.exceptions.Exit):
            non_fasta_file =  os.path.join(test_data, 'non_fasta.txt')
            ctx = click.Context(click.Command('cmd'), obj={'prop': 'A Context'})
            validation.validate_fasta(non_fasta_file, ctx)
