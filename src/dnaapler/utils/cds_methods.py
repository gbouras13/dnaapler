import os
import random
from pathlib import Path

import pyrodigal
from Bio import SeqIO
from loguru import logger

from dnaapler.utils.constants import DNAAPLER_DB
from dnaapler.utils.external_tools import ExternalTool
from dnaapler.utils.processing import (
    process_blast_output_and_reorient,
    reorient_sequence_random,
)


def run_mystery(ctx, input: Path, seed_value: int, output: Path, prefix: str) -> None:
    """
    Reorients a DNA sequence with a random Coding DNA Sequence (CDS) from the input file.

    Args:
        ctx (object): A context or environment object, possibly from a command-line tool.
        input (str): Path to the input DNA sequence file in FASTA format.
        seed_value (int): Seed value for randomization.
        output (str): Path to the output directory where the reoriented sequence will be saved.
        prefix (str): Prefix for the output file name.

    Returns:
        None

    This function searches for CDS using Pyrodigal in the input DNA sequence, then randomly selects
    one CDS (not the first or last) and reorients the DNA sequence with respect to that CDS.
    It saves the reoriented sequence in the specified output directory with the given prefix.
    """

    # get number of records of input
    logger.info("Searching for CDS with pyrodigal")
    orf_finder = pyrodigal.GeneFinder(meta=True)

    # set seed
    random.seed(int(seed_value))

    # there will only be 1 record
    for i, record in enumerate(SeqIO.parse(input, "fasta")):
        genes = orf_finder.find_genes(str(record.seq))
        # get number of genes
        gene_count = len(genes)

        # ensure has > 3 genes
        if gene_count < 4:
            logger.error(
                f"{input} has less than 4 CDS. You probably shouldn't be using dnaapler mystery!"
            )
            ctx.exit(2)

        logger.info("Reorienting with a random CDS (that is not the first or last).")

        # ensure not first or last gene
        reorient_gene_number = random.randint(2, gene_count - 1)

        logger.info(f"Gene number {reorient_gene_number} was selected.")
        start = genes[reorient_gene_number].begin
        strand = genes[reorient_gene_number].strand

        if strand == 1:
            strand_eng = "forward"
        else:
            strand_eng = "negative"

        logger.info(f"Your random CDS has a start coordinate of {start}.")
        logger.info(f"Your random CDS is on the {strand_eng} strand.")

        output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
        reorient_sequence_random(input, output_processed_file, start, strand)


def run_nearest(ctx, input: Path, output: Path, prefix: str) -> None:
    """
    Reorients a DNA sequence to begin with the first Coding DNA Sequence (CDS) found in the input file.

    Args:
        ctx (object): A context or environment object, possibly from a command-line tool.
        input (str): Path to the input DNA sequence file in FASTA format.
        output (str): Path to the output directory where the reoriented sequence will be saved.
        prefix (str): Prefix for the output file name.

    Returns:
        None

    This function searches for CDS using Pyrodigal in the input DNA sequence, then reorients the DNA
    sequence to begin with the first CDS found. The reoriented sequence is saved in the specified
    output directory with the given prefix.
    """
    # get number of records of input
    logger.info("Searching for CDS with pyrodigal")
    orf_finder = pyrodigal.GeneFinder(meta=True)

    # there will only be 1 record
    for i, record in enumerate(SeqIO.parse(input, "fasta")):
        genes = orf_finder.find_genes(str(record.seq))
        # get number of genes
        gene_count = len(genes)

        # ensure has > 1 genes
        if gene_count < 2:
            logger.error(
                f"{input} has less than 2 CDS. You probably shouldn't be using dnaapler nearest!"
            )
            ctx.exit(2)

        logger.info("Reorienting to begin with the first CDS.")

        reorient_gene_number = 1

        start = genes[reorient_gene_number].begin
        strand = genes[reorient_gene_number].strand

        if strand == 1:
            strand_eng = "forward"
        else:
            strand_eng = "negative"

        logger.info(f"Your first CDS has a start coordinate of {start}.")
        logger.info(f"Your first CDS is on the {strand_eng} strand.")

        output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
        reorient_sequence_random(input, output_processed_file, start, strand)


def run_largest(ctx, input: Path, output: Path, prefix: str) -> None:
    """
    Reorients a DNA sequence to begin with the largest Coding DNA Sequence (CDS) found in the input file.

    Args:
        ctx (object): A context or environment object, possibly from a command-line tool.
        input (Path): Path to the input DNA sequence file in FASTA format.
        seed_value (int): Seed value for randomization.
        output (Path): Path to the output directory where the reoriented sequence will be saved.
        prefix (str): Prefix for the output file name.

    Returns:
        None

    This function searches for CDS using Pyrodigal in the input DNA sequence, then reorients the DNA
    sequence to begin with the largest CDS found based on the size of the coding region. The reoriented
    sequence is saved in the specified output directory with the given prefix.
    """
    # get number of records of input
    logger.info("Searching for CDS with pyrodigal")
    orf_finder = pyrodigal.GeneFinder(meta=True)

    # there will only be 1 record
    for i, record in enumerate(SeqIO.parse(input, "fasta")):
        genes = orf_finder.find_genes(str(record.seq))
        # get number of genes
        gene_count = len(genes)

        # ensure has > 3 genes
        if gene_count < 4:
            logger.error(
                f"{input} has less than 4 CDS. You probably shouldn't be using dnaapler mystery!"
            )

        logger.info("Reorienting with the largest CDS.")

        # iterate over dict

        size_dict = {}
        gene_index = 0

        for gene in genes:
            size = abs(gene.end - gene.begin)
            size_dict[gene_index] = size
            gene_index += 1

        # Find the gene with the max overlap
        largest_gene_index = max(size_dict, key=lambda key: size_dict[key])

        start = genes[largest_gene_index].begin
        strand = genes[largest_gene_index].strand
        max_size = size_dict[largest_gene_index] / 3

        if strand == 1:
            strand_eng = "forward"
        else:
            strand_eng = "negative"

        logger.info(
            f"Your largest CDS has a start coordinate of {start} and has a size of {max_size} AAs."
        )
        logger.info(f"Your largest CDS is on the {strand_eng} strand.")

        output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
        reorient_sequence_random(input, output_processed_file, start, strand)


def run_blast_based_method(
    ctx, input: Path, output: Path, prefix: str, gene: str, evalue: float, threads: int
) -> bool:
    """
    Run a BLAST-based approach to reorient a DNA sequence based on a specific gene.

    Args:
        ctx (object): A context or environment object, possibly from a command-line tool.
        input (Path): Path to the input DNA sequence file in FASTA format.
        output (Path): Path to the output directory where the reoriented sequence will be saved.
        prefix (str): Prefix for the output file name.
        gene (str): Name of the gene used for BLAST search (options: 'dnaA', 'repA', 'terL', 'custom').
        evalue (float): E-value threshold for BLAST search.
        threads (int): Number of threads for BLAST search.

    Returns:
        bool: True if the BLAST-based approach succeeded, False otherwise.
    """

    # sets DB directory based of the gene name

    # defines db name
    db_name = "dnaA_db"
    if gene == "dnaA":
        db_name = "dnaA_db"
    elif gene == "repA":
        db_name = "repA_db"
    elif gene == "terL":
        db_name = "terL_db"

    # chromosome path
    # blast
    logdir = Path(f"{output}/logs")
    blast_output = os.path.join(output, f"{prefix}_blast_output.txt")

    db = os.path.join(DNAAPLER_DB, db_name)
    if gene == "custom":
        db = os.path.join(output, "custom_db", "custom_db")
    blast = ExternalTool(
        tool="blastx",
        input=f"-query {input}",
        output=f"-out {blast_output}",
        params=f'-db {db} -evalue  {evalue} -num_threads {str(threads)} -outfmt " 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq "',
        logdir=logdir,
    )

    ExternalTool.run_tool(blast, ctx)

    # reorient the genome based on the BLAST hit
    output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
    blast_success = process_blast_output_and_reorient(
        input, blast_output, output_processed_file, gene
    )

    return blast_success


def run_blast_based_method_bulk(
    ctx, input: Path, output: Path, prefix: str, gene: str, evalue: float, threads: int
) -> bool:
    """
    Run a bulk BLAST-based approach to reorient multiple DNA sequences based on a specific gene.

    Args:
        ctx (object): A context or environment object, possibly from a command-line tool.
        input (Path): Path to the directory containing input DNA sequence files in FASTA format.
        output (Path): Path to the output directory where the reoriented sequences will be saved.
        prefix (str): Prefix for the output file names.
        gene (str): Name of the gene used for BLAST search (options: 'dnaA', 'repA', 'terL', 'custom').
        evalue (float): E-value threshold for BLAST search.
        threads (int): Number of threads for BLAST search.

    Returns:
        bool: True if the bulk BLAST-based approach succeeded for all sequences, False otherwise.
    """

    # sets DB directory based of the gene name

    # defines db name
    db_name = "dnaA_db"
    if gene == "dnaA":
        db_name = "dnaA_db"
    elif gene == "repA":
        db_name = "repA_db"
    elif gene == "terL":
        db_name = "terL_db"

    # chromosome path
    # blast
    logdir = Path(f"{output}/logs")
    blast_output = os.path.join(output, f"{prefix}_blast_output.txt")

    db = os.path.join(DNAAPLER_DB, db_name)
    if gene == "custom":
        db = os.path.join(output, "custom_db", "custom_db")
    blast = ExternalTool(
        tool="blastx",
        input=f"-query {input}",
        output=f"-out {blast_output}",
        params=f'-db {db} -evalue  {evalue} -num_threads {str(threads)} -outfmt " 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq "',
        logdir=logdir,
    )

    ExternalTool.run_tool(blast, ctx)

    # reorient the genome based on the BLAST hit
    output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
    blast_success = process_blast_output_and_reorient(
        input, blast_output, output_processed_file, gene
    )

    return blast_success
