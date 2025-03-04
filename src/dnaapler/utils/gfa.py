"""
Contains functions related to GFA processing.
Uses some code from https://github.com/rrwick/Circular-Contig-Extractor
"""

import os
import re
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


def prep_gfa(input_file, output_dir):
    """
    If the input file given to Dnaapler is a GFA file, this function is run early in Dnaapler's
    pipeline. It will save a temporary FASTA file which contains the circular sequences from the
    GFA.

    Returns:
    * bool: whether or not the input was GFA format
    * str: FASTA input file to reorient (if the input was GFA, this is the temp FASTA file, but if
    *      the input was FASTA this is just the same FASTA)
    * str: GFA input file (if the input was FASTA this is None)
    """
    if is_gfa(input_file):
        temp_input_fasta = os.path.join(output_dir, "input.fasta")
        save_circular_sequences_as_fasta(input_file, temp_input_fasta)
        return True, temp_input_fasta, input_file
    else:
        return False, input_file, None


def finalise_gfa(temp_input_fasta, gfa_input_file, output_fasta):
    """
    If the input file given to Dnaapler is a GFA file, this function is run at the end of Dnaapler's
    pipeline. It will create the output GFA and remove the output FASTA (because non-circular
    sequences will be missing from the FASTA).
    """
    remove_file(Path(temp_input_fasta))
    save_reoriented_gfa(gfa_input_file, output_fasta)
    remove_file(Path(output_fasta))


def remove_file(file_path: Path):
    if file_path.exists():
        file_path.unlink()


def save_circular_sequences_as_fasta(gfa_file, fasta_file):
    """
    Identifies the circular sequences in the GFA file and saves them in FASTA format.
    """
    contigs, links = load_gfa(gfa_file)
    contigs = find_circular_contigs(contigs, links)
    if not contigs:
        logger.error(f"Error: {gfa_file} file contains no circular sequences.")
    contigs = trim_overlaps(contigs)
    write_fasta(contigs, fasta_file)
    logger.info(f"number of circular sequences in {gfa_file}: {len(contigs)}")


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
