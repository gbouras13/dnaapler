import os
import random
import shutil
from pathlib import Path

import random
import pyrodigal
from Bio import SeqIO
from loguru import logger

from dnaapler.utils.processing import (
    reorient_sequence_random,
)

def run_mystery(ctx, input, seed_value, output, prefix ):
# get number of records of input
    orf_finder = pyrodigal.OrfFinder(meta=True)

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