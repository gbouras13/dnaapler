#!/usr/bin/env python3
"""
Test script to verify the process_ignore_input functionality works with both
comma-separated lists, file paths, and stdin.
"""

import os
import sys
import tempfile
import unittest
from io import StringIO
from unittest.mock import patch

# Add the src directory to the path so we can import dnaapler
sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
)

from dnaapler.utils.validation import process_ignore_input  # noqa: E402


class TestProcessIgnoreInput(unittest.TestCase):

    def test_empty_input(self):
        """Test empty string and None input"""
        self.assertEqual(process_ignore_input(""), [])
        self.assertEqual(process_ignore_input(None), [])

    def test_comma_separated_basic(self):
        """Test basic comma-separated list functionality"""
        result = process_ignore_input("chr1,chr2,chr3")
        expected = ["chr1", "chr2", "chr3"]
        self.assertEqual(result, expected)

    def test_comma_separated_with_spaces(self):
        """Test comma-separated with spaces around commas"""
        result = process_ignore_input("chr1, chr2 , chr3")
        expected = ["chr1", "chr2", "chr3"]
        self.assertEqual(result, expected)

    def test_single_chromosome(self):
        """Test single chromosome name"""
        result = process_ignore_input("chr1")
        expected = ["chr1"]
        self.assertEqual(result, expected)

    def test_comma_separated_vs_path_detection(self):
        """Test that comma-separated lists without path characteristics work"""
        # These should be treated as comma-separated lists (no path separators or extensions)
        result = process_ignore_input("chr1,chr2,chr3")
        expected = ["chr1", "chr2", "chr3"]
        self.assertEqual(result, expected)

        result = process_ignore_input("contig_1,contig_2")
        expected = ["contig_1", "contig_2"]
        self.assertEqual(result, expected)

    def test_file_path_with_multiple_entries(self):
        """Test file path with multiple entries per line (space-separated)"""
        # Create a temporary file with test data (similar to existing test data)
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as tmp:
            tmp.write("MW460250_1_subset\n")
            tmp.write(
                "1 length=181436 plasmid_copy_number_short=1.0x plasmid_copy_number_long=1.04x circular=true\n"
            )
            tmp_path = tmp.name

        try:
            result = process_ignore_input(tmp_path)
            expected_set = {
                "MW460250_1_subset",
                "1",
            }  # Should split by space and take first part
            result_set = set(result)
            self.assertEqual(result_set, expected_set)
            self.assertEqual(len(result), 2)
        finally:
            os.unlink(tmp_path)

    def test_file_path_simple(self):
        """Test file path with simple entries (one per line)"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as tmp:
            tmp.write("chr1\n")
            tmp.write("chr2\n")
            tmp.write("chr3\n")
            tmp_path = tmp.name

        try:
            result = process_ignore_input(tmp_path)
            expected_set = {"chr1", "chr2", "chr3"}
            result_set = set(result)
            self.assertEqual(result_set, expected_set)
            self.assertEqual(len(result), 3)
        finally:
            os.unlink(tmp_path)

    def test_nonexistent_file_path(self):
        """Test handling of nonexistent file path (should exit with error)"""
        with self.assertRaises(SystemExit):
            process_ignore_input("/path/that/does/not/exist.txt")

    def test_actual_test_data_file(self):
        """Test with the actual ignore.txt file from test data"""
        ignore_file = os.path.join(
            os.path.dirname(__file__), "test_data", "overall_inputs", "ignore.txt"
        )
        if os.path.exists(ignore_file):
            result = process_ignore_input(ignore_file)
            # Should contain the entries from the file
            self.assertGreater(len(result), 0)
            self.assertIn("MW460250_1_subset", result)
            self.assertIn("1", result)
        else:
            self.skipTest("Actual test file not found")

    def test_path_detection(self):
        """Test that nonexistent paths cause system exit"""
        # These should be treated as file paths and cause system exit since files don't exist
        test_paths = [
            "/path/to/file.txt",
            "relative/path/file.tsv",
            "file.csv",
            "C:\\Windows\\path\\file.list",
        ]

        for path in test_paths:
            with self.assertRaises(
                SystemExit, msg=f"Path {path} should cause system exit"
            ):
                process_ignore_input(path)

    def test_stdin_input(self):
        """Test reading from stdin using '-' as input"""
        test_input = "chr1\nchr2\nchr3\n"
        with patch("sys.stdin", StringIO(test_input)):
            result = process_ignore_input("-")
            expected_set = {"chr1", "chr2", "chr3"}
            result_set = set(result)
            self.assertEqual(result_set, expected_set)
            self.assertEqual(len(result), 3)

    def test_stdin_with_spaces(self):
        """Test stdin input with space-separated content (should take first part)"""
        test_input = (
            "MW460250_1_subset\n1 length=181436 plasmid_copy_number_short=1.0x\n"
        )
        with patch("sys.stdin", StringIO(test_input)):
            result = process_ignore_input("-")
            expected_set = {"MW460250_1_subset", "1"}
            result_set = set(result)
            self.assertEqual(result_set, expected_set)
            self.assertEqual(len(result), 2)

    def test_stdin_empty(self):
        """Test stdin with empty input"""
        test_input = ""
        with patch("sys.stdin", StringIO(test_input)):
            result = process_ignore_input("-")
            self.assertEqual(result, [])


if __name__ == "__main__":
    unittest.main()
