"""
Contains functions related to GFA processing.
Uses some code from https://github.com/rrwick/Circular-Contig-Extractor
"""

import os
import re
import shutil
from pathlib import Path

from Bio import SeqIO
from loguru import logger


def is_gfa(input_file):
    """
    Check if the file is in GFA format.
    """
    with open(input_file, "r") as handle:
        first_char = handle.read(1)
        if first_char not in {"H", "S", "L"}:
            return False
    return gfa_sequence_count(input_file) > 0


def prep_gfa(input_file, output_dir, prefix, summary_type=None):
    """
    If the input file given to Dnaapler is a GFA file, this function is run early in Dnaapler's
    pipeline. It saves a temporary FASTA file which contains the circular sequences from the GFA,
    which is what Dnaapler reorients.

    If the GFA contains no circular sequences there is nothing to reorient, so the input GFA is
    instead copied to the output and all sequences are written out as a linear FASTA, and the caller
    is told to skip reorientation.

    Args:
    * summary_type (str | None): for a GFA with no circular sequences, controls whether a
    *      reorientation summary is also written. "all" and "bulk" write the corresponding summary
    *      (with every contig marked as not reoriented); None (the single-gene commands) writes none.

    Returns:
    * bool: whether or not the input was GFA format
    * str: FASTA input file to reorient (if the input was GFA, this is the temp FASTA file, but if
    *      the input was FASTA this is just the same FASTA)
    * str: GFA input file (if the input was FASTA this is None)
    * bool: whether reorientation should be skipped (True only for a GFA with no circular sequences,
    *       in which case the output GFA and FASTA have already been written)
    """
    if not is_gfa(input_file):
        return False, input_file, None, False

    contigs, links = load_gfa(input_file)
    circular_contigs = find_circular_contigs(contigs, links)

    if not circular_contigs:
        logger.warning(
            f"{input_file} contains no circular sequences. No contigs will be reoriented; the "
            "input GFA will be copied to the output and all sequences written out as a linear FASTA."
        )
        write_gfa_passthrough(input_file, output_dir, prefix, summary_type)
        return True, input_file, input_file, True

    temp_input_fasta = os.path.join(output_dir, "input.fasta")
    circular_contigs = trim_overlaps(circular_contigs)
    write_fasta(circular_contigs, temp_input_fasta)
    logger.info(
        f"number of circular sequences in {input_file}: {len(circular_contigs)}"
    )
    return True, temp_input_fasta, input_file, False


def finalise_gfa(temp_input_fasta, gfa_input_file, output_fasta):
    """
    If the input file given to Dnaapler is a GFA file, this function is run at the end of Dnaapler's
    pipeline. It creates the output GFA and rewrites the output FASTA so that it contains all
    contigs from the GFA (circular contigs reoriented, non-circular contigs passed through
    unchanged), making it suitable for downstream tools such as polishers.
    """
    remove_file(Path(temp_input_fasta))
    # save_reoriented_gfa reads the circular-only output_fasta to build the GFA, so it must run
    # before we overwrite output_fasta with the complete (all-contigs) version below.
    reoriented_gfa = save_reoriented_gfa(gfa_input_file, output_fasta)
    gfa_to_fasta(reoriented_gfa, output_fasta)


def remove_file(file_path: Path):
    if file_path.exists():
        file_path.unlink()


def write_gfa_passthrough(input_gfa, output_dir, prefix, summary_type=None):
    """
    Handles a GFA input that contains no circular sequences: copies the input GFA to the output
    {prefix}_reoriented.gfa and writes all of its sequences out, unchanged, as a linear
    {prefix}_reoriented.fasta. Nothing is reoriented.

    For the `all` and `bulk` commands (summary_type "all"/"bulk"), a reorientation summary is also
    written, with every contig marked as not reoriented.
    """
    reoriented_gfa = os.path.join(output_dir, f"{prefix}_reoriented.gfa")
    reoriented_fasta = os.path.join(output_dir, f"{prefix}_reoriented.fasta")
    logger.info(f"copying input GFA to {reoriented_gfa}")
    shutil.copyfile(input_gfa, reoriented_gfa)
    gfa_to_fasta(input_gfa, reoriented_fasta)
    if summary_type is not None:
        write_no_reorientation_summary(
            gfa_sequence_names(input_gfa), output_dir, prefix, summary_type
        )


def write_no_reorientation_summary(contig_names, output_dir, prefix, summary_type):
    """
    Writes a reorientation summary TSV marking every contig as not reoriented. Used for the `all`
    and `bulk` commands when a GFA input has no circular sequences. The columns mirror the summaries
    written by dnaapler.utils.all and dnaapler.utils.bulk respectively.
    """
    if summary_type == "all":
        columns = [
            "Contig",
            "Gene_Reoriented",
            "Start",
            "Strand",
            "Top_Hit",
            "Top_Hit_Length",
            "Covered_Length",
            "Coverage",
            "Identical_AAs",
            "Identity_Percentage",
            "Overlapping_Contig_End",
        ]
        summary_file = os.path.join(
            output_dir, f"{prefix}_all_reorientation_summary.tsv"
        )
    elif summary_type == "bulk":
        columns = [
            "Contig",
            "Start",
            "Strand",
            "Top_Hit",
            "Top_Hit_Length",
            "Covered_Length",
            "Coverage",
            "Identical_AAs",
            "Identity_Percentage",
        ]
        summary_file = os.path.join(
            output_dir, f"{prefix}_bulk_reorientation_summary.tsv"
        )
    else:
        return

    logger.info(f"writing reorientation summary to {summary_file}")
    with open(summary_file, "wt") as f:
        f.write("\t".join(columns) + "\n")
        for name in contig_names:
            f.write(
                "\t".join([name] + ["No_reorientation"] * (len(columns) - 1)) + "\n"
            )


def load_gfa(filename):
    contigs, links = [], []
    with open(filename, "rt") as gfa_file:
        for line in gfa_file:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "S":
                contigs.append(parts[1:3])
            if parts[0] == "L":
                links.append(parts[1:6])
    return contigs, links


def gfa_sequence_count(filename):
    return len(load_gfa(filename)[0])


def gfa_sequence_names(filename):
    return [contig[0] for contig in load_gfa(filename)[0]]


def find_circular_contigs(contigs, links):
    """
    Returns a list of contigs with a simple circular structure (one circularising link and no other
    links). The return list contains tuples of (name, sequence, cigar).
    """
    circular_links = {}
    for seg_a, strand_a, seg_b, strand_b, cigar in links:
        if seg_a == seg_b and strand_a == strand_b:
            circular_links[seg_a] = cigar
    for seg_a, strand_a, seg_b, strand_b, _ in links:
        if seg_a != seg_b or strand_a != strand_b:
            circular_links.pop(seg_a, None)
            circular_links.pop(seg_b, None)
    circular_contigs = []
    for name, seq in contigs:
        if name in circular_links:
            circular_contigs.append((name, seq, circular_links[name]))
    return circular_contigs


def trim_overlaps(contigs):
    trimmed_contigs = []
    for name, seq, cigar in contigs:
        overlap = get_overlap_from_cigar(cigar)
        if overlap is None:
            logger.error(f"Error: cannot determine overlap from CIGAR string {cigar}")
        trimmed_contigs.append((name, trim_seq(seq, overlap)))
    return trimmed_contigs


def get_overlap_from_cigar(cigar):
    match = re.match(r"^(\d+)M$", cigar)
    return int(match.group(1)) if match else None


def trim_seq(seq, trim_amount):
    if trim_amount is None or trim_amount == 0:
        return seq
    else:
        return seq[:-trim_amount]


def write_fasta(contigs, filename):
    with open(filename, "wt") as f:
        for name, seq in contigs:
            f.write(f">{name}\n")
            f.write(f"{seq}\n")


def save_reoriented_gfa(original_gfa, reoriented_fasta):
    """
    Copies the original input GFA file to an output GFA file, replacing sequences with their
    reoriented versions when possible. Since the reoriented sequences may have had overlap trimmed
    off, it will also modify the CIGAR strings for circularising links.
    """
    assert reoriented_fasta.endswith("_reoriented.fasta")
    reoriented_gfa = reoriented_fasta[:-17] + "_reoriented.gfa"
    logger.info(f"saving reoriented sequences to GFA format in {reoriented_gfa}")
    reoriented_seqs, reoriented_genes = load_reoriented_fasta(reoriented_fasta)
    with open(original_gfa, "rt") as in_gfa, open(reoriented_gfa, "wt") as out_gfa:
        for line in in_gfa:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "S" and parts[1] in reoriented_seqs:
                parts[2] = reoriented_seqs[parts[1]]
                if parts[1] in reoriented_genes:
                    parts.append(f"RT:z:{reoriented_genes[parts[1]]}")
                line = "\t".join(parts) + "\n"
            elif (
                parts[0] == "L" and parts[1] == parts[3] and parts[1] in reoriented_seqs
            ):
                parts[5] = "0M"
                line = "\t".join(parts) + "\n"
            out_gfa.write(line)
    return reoriented_gfa


def gfa_to_fasta(gfa_file, fasta_file):
    """
    Writes all sequences (the S lines) from a GFA file to a FASTA file. Contigs that were reoriented
    carry an "RT:z:<gene>" tag in the GFA (added by save_reoriented_gfa); these are annotated in the
    FASTA header with "rotated=True rotated_gene=<gene>" to match the `dnaapler all` convention.
    Non-circular contigs are written out unchanged with a plain header.
    """
    logger.info(f"saving reoriented sequences to FASTA format in {fasta_file}")
    with open(gfa_file, "rt") as in_gfa, open(fasta_file, "wt") as out_fasta:
        for line in in_gfa:
            parts = line.rstrip("\n").split("\t")
            if parts[0] != "S":
                continue
            name, seq = parts[1], parts[2]
            gene = None
            for field in parts[3:]:
                if field.startswith("RT:z:"):
                    gene = field[len("RT:z:") :]
                    break
            if gene is not None:
                header = f">{name} rotated=True rotated_gene={gene}"
            else:
                header = f">{name}"
            out_fasta.write(f"{header}\n{seq}\n")


def load_reoriented_fasta(reoriented_fasta):
    """
    Reads the Dnaapler reoriented FASTA file and returns two dictionaries:
    1. names -> sequences
    2. names -> rotated genes (if any)
    """
    seq_dict = {}
    gene_dict = {}
    for record in SeqIO.parse(reoriented_fasta, "fasta"):
        seq_dict[record.id] = str(record.seq)
        match = re.search(r"rotated_gene=([^ \t]+)", record.description)
        if match:
            gene_dict[record.id] = match.group(1)
    return seq_dict, gene_dict
