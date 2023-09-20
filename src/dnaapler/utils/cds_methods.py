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


def run_mystery(ctx, input, seed_value, output, prefix):
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


def run_nearest(ctx, input, output, prefix):
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


def run_blast_based_method(ctx, input, output, prefix, gene, evalue, threads):
    """
    returns: bool -  blast_success, whether or not the BLAST based approach succeeded
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
        params=f'-db {db} -evalue  {evalue} -num_threads {threads} -outfmt " 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq "',
        logdir=logdir,
    )

    ExternalTool.run_tool(blast, ctx)

    # reorient the genome based on the BLAST hit
    output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
    blast_success = process_blast_output_and_reorient(
        input, blast_output, output_processed_file, gene
    )

    return blast_success


def run_blast_based_method_bulk(ctx, input, output, prefix, gene, evalue, threads):
    """
    returns: bool -  blast_success, whether or not the BLAST based approach succeeded
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
        params=f'-db {db} -evalue  {evalue} -num_threads {threads} -outfmt " 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq "',
        logdir=logdir,
    )

    ExternalTool.run_tool(blast, ctx)

    # reorient the genome based on the BLAST hit
    output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
    blast_success = process_blast_output_and_reorient(
        input, blast_output, output_processed_file, gene
    )

    return blast_success
