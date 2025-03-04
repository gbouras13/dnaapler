"""
Unit tests for dnaapler overall

Usage: pytest .

"""

# import
import os
import shutil

# import functions
import subprocess
import unittest
from pathlib import Path

import pytest

test_data = Path("tests/test_data")
overall_test_data = Path(f"{test_data}/overall_inputs")
contig_end_data = Path(f"{test_data}/contig_ends")


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_phage(tmp_dir):
    """test phage"""
    input_fasta: Path = f"{overall_test_data}/NC_007458.fasta"
    cmd = f"dnaapler phage -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_phage_start_codon_not_found(tmp_dir):
    """test phage where the tophit has no start codon in alignment"""
    input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
    cmd = f"dnaapler phage -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_chrom(tmp_dir):
    """test chrom"""
    input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
    cmd = f"dnaapler chromosome -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_all_archaea(tmp_dir):
    """test all archaea"""
    input_fasta: Path = f"{overall_test_data}/CP001742.1_archaea.fasta"
    cmd = f"dnaapler all -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_archaea(tmp_dir):
    """test archaea"""
    input_fasta: Path = f"{overall_test_data}/CP001742.1_archaea.fasta"
    cmd = f"dnaapler archaea -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_chrom_start_codon_not_found(tmp_dir):
    """test chrom"""
    input_fasta: Path = f"{overall_test_data}/chromosome_top_hit_no_start_codon.fasta"
    cmd = f"dnaapler chromosome -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_chrom_diff_eval(tmp_dir):
    """test chrom with different e value"""
    input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
    cmd = f"dnaapler chromosome -i {input_fasta} -o {tmp_dir} -t 1 -e 0.1 -f"
    exec_command(cmd)


def test_chrom_mystery_autocomplete(tmp_dir):
    """test chrom - with phage (so no hit)"""
    input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
    cmd = f"dnaapler chromosome -i {input_fasta} -o {tmp_dir} -t 1 -f -a mystery"
    exec_command(cmd)


def test_chrom_nearest_autocomplete(tmp_dir):
    """test chrom - with phage (so no hit)"""
    input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
    cmd = f"dnaapler chromosome -i {input_fasta} -o {tmp_dir} -t 1 -f -a nearest"
    exec_command(cmd)


def test_plas(tmp_dir):
    """test plas"""
    input_fasta: Path = f"{overall_test_data}/plasmid.fasta"
    cmd = f"dnaapler plasmid -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_plas_top_hit_no_start_codon(tmp_dir):
    """test plas tophit no start codon"""
    input_fasta: Path = f"{overall_test_data}/plasmid_top_hit_no_start_codon.fasta"
    cmd = f"dnaapler plasmid -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_plas_nearest_autocomplete(tmp_dir):
    """test plas - with phage (so no hit)"""
    input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
    cmd = f"dnaapler plasmid -i {input_fasta} -o {tmp_dir} -t 1 -f -a nearest"
    exec_command(cmd)


def test_plas_mystery_autocomplete(tmp_dir):
    """test plas - with phage (so no hit)"""
    input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
    cmd = f"dnaapler plasmid -i {input_fasta} -o {tmp_dir} -t 1 -f -a mystery"
    exec_command(cmd)


def test_phage_mystery(tmp_dir):
    """test phage mystery - with plasmid (so no hit)"""
    input_fasta: Path = f"{overall_test_data}/plasmid.fasta"
    cmd = f"dnaapler phage -i {input_fasta} -o {tmp_dir} -t 1 -f -a mystery"
    exec_command(cmd)


def test_phage_nearest(tmp_dir):
    """test phage nearest - with plasmid (so no hit)"""
    input_fasta: Path = f"{overall_test_data}/plasmid.fasta"
    cmd = f"dnaapler phage -i {input_fasta} -o {tmp_dir} -t 1 -f -a nearest"
    exec_command(cmd)


def test_mys(tmp_dir):
    """test mystery"""
    input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
    cmd = f"dnaapler mystery -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_nearest(tmp_dir):
    """test nearest"""
    input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
    cmd = f"dnaapler nearest -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_largest(tmp_dir):
    """test largest"""
    input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
    cmd = f"dnaapler largest -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_custom(tmp_dir):
    """test custom"""
    input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
    custom: Path = f"{test_data}/fake_custom.faa"
    cmd = f"dnaapler custom -i {input_fasta} -o {tmp_dir} -c {custom} -t 1 -f"
    exec_command(cmd)


def test_bulk_chromosome(tmp_dir):
    """test bulk chromosome"""
    input_fasta: Path = f"{overall_test_data}/bulk_chromosome.fasta"
    cmd = f"dnaapler bulk -m chromosome -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_bulk_phage(tmp_dir):
    """test bulk phage"""
    input_fasta: Path = f"{overall_test_data}/bulk_phage.fasta"
    cmd = f"dnaapler bulk -m phage -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_bulk_plasmid(tmp_dir):
    """test bulk plasmid"""
    input_fasta: Path = f"{overall_test_data}/bulk_plasmid.fasta"
    cmd = f"dnaapler bulk -m plasmid -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_bulk_archaea(tmp_dir):
    """test bulk archaea"""
    input_fasta: Path = f"{overall_test_data}/bulk_archaea.fasta"
    cmd = f"dnaapler bulk -m archaea -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_bulk_custom(tmp_dir):
    """test bulk custom"""
    input_fasta: Path = f"{overall_test_data}/bulk_chromosome.fasta"
    custom: Path = f"{test_data}/fake_custom.faa"
    cmd = f"dnaapler bulk -m custom -i {input_fasta} -o {tmp_dir} -c {custom} -t 1 -f"
    exec_command(cmd)


def test_all_custom(tmp_dir):
    """test all"""
    input_fasta: Path = f"{overall_test_data}/all_test.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_custom(tmp_dir):
    """test custom"""
    input_fasta: Path = f"{overall_test_data}/all_test.fasta"
    custom: Path = f"{test_data}/fake_custom.faa"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f -c {custom}"
    exec_command(cmd)


def test_no_overlap(tmp_dir):
    """test no overlap"""
    input_fasta: Path = f"{overall_test_data}/PP04977_pilon_subset.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_all_dnaa_repa(tmp_dir):
    """test all dnaa repa"""
    input_fasta: Path = f"{overall_test_data}/all_test.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f --db dnaa,repa"
    exec_command(cmd)


def test_gfa_all(tmp_dir):
    """test all with a GFA input - seq 3 contains overlap and seq 5 in non-circular"""
    input_fasta: Path = f"{overall_test_data}/all_test.gfa"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)
    output_gfa_lines = open(f"{tmp_dir}/dnaapler_reoriented.gfa").read().splitlines()
    seq_1 = output_gfa_lines[1].split("\t")[2]
    seq_2 = output_gfa_lines[4].split("\t")[2]
    seq_3 = output_gfa_lines[7].split("\t")[2]
    seq_4 = output_gfa_lines[10].split("\t")[2]
    seq_5 = output_gfa_lines[13].split("\t")[2]
    assert seq_1.startswith("GTGTCACTTT") and len(seq_1) == 398227
    assert seq_2.startswith("GTGACTGATC") and len(seq_2) == 186230
    assert seq_3.startswith("ATGATCGTAG") and len(seq_3) == 1240
    assert seq_4.startswith("ATGTCAGAAG") and len(seq_4) == 2088
    assert seq_5.startswith("GTTCTAGCAT") and len(seq_5) == 1000


"""
this one is for hybracter
"""


def test_all_dnaa_repa_cog1474(tmp_dir):
    """test all dnaa repa cog1474 - for hybracter"""
    input_fasta: Path = f"{overall_test_data}/all_test.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f --db dnaa,cog1474,repa"
    exec_command(cmd)


def test_all_dnaa_terl(tmp_dir):
    """test all dnaa terl"""
    input_fasta: Path = f"{overall_test_data}/all_reorient_and_no_reorient.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f --db dnaa,terl"
    exec_command(cmd)


def test_all_repa_terl(tmp_dir):
    """test all repa terl"""
    input_fasta: Path = f"{overall_test_data}/all_reorient_and_no_reorient.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f --db repa,terl"
    exec_command(cmd)


def test_all_ignore(tmp_dir):
    """test all"""
    input_fasta: Path = f"{overall_test_data}/all_test.fasta"
    ignore_file: Path = f"{overall_test_data}/ignore.txt"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f --ignore {ignore_file}"
    exec_command(cmd)


def test_all_no_reorientation(tmp_dir):
    """test all where there is already orientation or no orientation"""
    input_fasta: Path = f"{overall_test_data}/all_reorient_and_no_reorient.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_all_autocomplete_mystery(tmp_dir):
    """test all where autcompletion is required mystery"""
    input_fasta: Path = f"{overall_test_data}/all_test_autocomplete.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f -a mystery"
    exec_command(cmd)


def test_all_autocomplete_nearest(tmp_dir):
    """test all where autcompletion is required nearest"""
    input_fasta: Path = f"{overall_test_data}/all_test_autocomplete.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f -a nearest"
    exec_command(cmd)


def test_all_autocomplete_largest(tmp_dir):
    """test all where autcompletion is required largest"""
    input_fasta: Path = f"{overall_test_data}/all_test_autocomplete.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f -a largest"
    exec_command(cmd)


"""
new contig end reorientation
"""


def test_all_chromosome_contig_end(tmp_dir):
    """test all where dnaA found overlapping contig end - taken from https://github.com/gbouras13/dnaapler/issues/90"""
    input_fasta: Path = f"{contig_end_data}/GN04821.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_all_phage_contig_end(tmp_dir):
    """test all where terL found overlapping contig end"""
    input_fasta: Path = f"{contig_end_data}/C333_sa3int_phage.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_all_archaea_contig_end(tmp_dir):
    """test all where cog1474 found overlapping contig end"""
    input_fasta: Path = f"{contig_end_data}/CP001742.1_archaea.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_all_plasmid_contig_end(tmp_dir):
    """test all where repA found overlapping contig end"""
    input_fasta: Path = f"{contig_end_data}/plasmid.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_all_mixed_contig_end(tmp_dir):
    """test all where all found overlapping contig ends plus some contigs without end reorientation"""
    input_fasta: Path = f"{contig_end_data}/mixed.fasta"
    cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_citation():
    """test citation"""
    cmd = f"dnaapler citation"
    exec_command(cmd)


class TestExits(unittest.TestCase):
    """Tests of End to End common failures"""

    def test_chrom_double(self):
        """test chrom double to test force is working"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
            outdir: Path = f"{overall_test_data}/chrom_out"
            cmd = f"dnaapler chromosome -i {input_fasta} -o {outdir} -t 1"
            exec_command(cmd)
            exec_command(cmd)

    def test_phage_chrom(self):
        """test chrom with phage input - test where nothing found"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
            outdir: Path = f"{overall_test_data}/phage_out"
            cmd = f"dnaapler chromosome -i {input_fasta} -o {outdir} -t 1 -f"
            exec_command(cmd)

    def test_chrom_bad_autocomplete(self):
        """test chrom with bad autcomplete"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
            outdir: Path = f"{overall_test_data}/phage_out"
            cmd = (
                f"dnaapler chromosome -i {input_fasta} -o {outdir} -t 1 -f -a bad_auto"
            )
            exec_command(cmd)

    def test_chrom_already_oriented(self):
        """test chrom with already oriented input"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{overall_test_data}/SAOMS1_reoriented.fasta"
            outdir: Path = f"{overall_test_data}/chrom_out"
            cmd = f"dnaapler chromosome -i {input_fasta} -o {outdir} -t 1 -f"
            exec_command(cmd)

    def test_bulk_single_genome(self):
        """test bulk with single genome"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{overall_test_data}/SAOMS1_reoriented.fasta"
            outdir: Path = f"{overall_test_data}/chrom_out"
            cmd = f"dnaapler bulk -m phage -i {input_fasta} -o {outdir} -t 1 -f"
            exec_command(cmd)

    def test_bulk_bad_mode(self):
        """test bulk with incorrect -m value"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{overall_test_data}/SAOMS1_reoriented.fasta"
            outdir: Path = f"{overall_test_data}/chrom_out"
            cmd = f"dnaapler bulk -m bad -i {input_fasta} -o {outdir} -t 1 -f"
            exec_command(cmd)

    def test_bulk_custom_no_db(self):
        """test bulk custom with no db"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{overall_test_data}/bulk_chromosome.fasta"
            outdir: Path = f"{overall_test_data}/chrom_out"
            cmd = f"dnaapler bulk -m custom -i {input_fasta} -o {outdir} -t 1 -f"
            exec_command(cmd)

    def test_bulk_phage_with_chrom(self):
        """test bulk chromosome with phage genomes"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{overall_test_data}/bulk_phage.fasta"
            outdir: Path = f"{overall_test_data}/chrom_out"
            cmd = f"dnaapler bulk -m chromosome -i {input_fasta} -o {outdir} -t 1 -f"
            exec_command(cmd)

    def test_all_autocomplete_mystery_too_small(self):
        """test all where the autocompletion mystery fails as the contig has < 4 CDS"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{test_data}/no_hit_plasmid.fasta"
            outdir: Path = f"{overall_test_data}/plas_out"
            cmd = f"dnaapler all -i {input_fasta} -o {outdir} -t 1 -f -a mystery"
            exec_command(cmd)

    def test_all_autocomplete_largest_too_small(self):
        """test all where the autocompletion largest fails as the contig has < 4 CDS"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{test_data}/no_hit_plasmid.fasta"
            outdir: Path = f"{overall_test_data}/plas_out"
            cmd = f"dnaapler all -i {input_fasta} -o {outdir} -t 1 -f -a largest"
            exec_command(cmd)

    def test_all_autocomplete_nearest_no_hits(self):
        """test all with autcomplete no hits"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{test_data}/nucl_test.fna"
            outdir: Path = f"{overall_test_data}/bulk_out"
            cmd = f"dnaapler all -i {input_fasta} -o {outdir} -t 1 -f -a nearest "
            exec_command(cmd)

    def test_all_dupe_header(self):
        """test all with autcomplete no hits"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{test_data}/dupe_header.fasta"
            outdir: Path = f"{overall_test_data}/bulk_out"
            cmd = f"dnaapler all -i {input_fasta} -o {outdir} -t 1 -f  "
            exec_command(cmd)

    def test_bulk_dupe_header(self):
        """test bulk with autcomplete no hits"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{test_data}/dupe_header.fasta"
            outdir: Path = f"{overall_test_data}/bulk_out"
            cmd = f"dnaapler bulk -i {input_fasta} -o {outdir} -t 1 -f  "
            exec_command(cmd)

    def test_gfa_no_circular(self):
        """test all with a GFA input that contains no circular sequences"""
        with self.assertRaises(RuntimeError):
            input_gfa: Path = f"{overall_test_data}/no_circular.gfa"
            outdir: Path = f"{overall_test_data}/all_out"
            cmd = f"dnaapler all -i {input_gfa} -o {outdir} -t 1 -f  "
            exec_command(cmd)

    def test_gfa_inexact_overlap(self):
        """test all with a GFA input that contains an inexact overlap"""
        with self.assertRaises(RuntimeError):
            input_gfa: Path = f"{overall_test_data}/inexact_overlap.gfa"
            outdir: Path = f"{overall_test_data}/all_out"
            cmd = f"dnaapler all -i {input_gfa} -o {outdir} -t 1 -f  "
            exec_command(cmd)


remove_directory(f"{overall_test_data}/phage_out")
remove_directory(f"{overall_test_data}/chrom_out")
remove_directory(f"{overall_test_data}/bulk_out")
remove_directory(f"{overall_test_data}/plas_out")
remove_directory(f"{overall_test_data}/all_out")
