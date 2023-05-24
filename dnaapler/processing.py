import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from loguru import logger
import click


def process_blast_output_and_reorient(input, blast_file, out_file, ctx: click.Context, gene: str):
    col_list = ["qseqid", "qlen", "sseqid", "slen", "length", "qstart", "qend", "sstart", "send", "pident", "nident", "gaps", "mismatch", "evalue", "bitscore", "qseq", "sseq"] 
    blast_df = pd.read_csv(blast_file, delimiter= '\t', index_col=False , names=col_list) 

    # take the top hit - blast sorts by bitscore
    # if the start or end coordinate is 1
    if blast_df['qstart'][0] == 1 or blast_df['qend'][0] == 1:
        logger.error(f"Based on the BLAST output, your input is already oriented to begin with {gene}.")
        ctx.exit(2)

    # prokaryotes can use AUG M, GUG V or UUG L as start codons - for Pseudomonas aeruginosa PA01  dnaA actually starts with V
    # Therefore, require the start codon to be the 1st in the searched sequence match - unlikely to not be but nonetheless
    elif blast_df['qseq'][0][0] in ['M', 'V', 'L'] and (blast_df['sstart'][0] == 1 or blast_df['send'][0] == 1):
        logger.info(f"{gene} identified. Re-orienting.")
        dnaa_start_on_chromosome = blast_df['qstart'][0]
        strand = "fwd"
        if blast_df['qstart'][0] > blast_df['qend'][0]:
            strand = "rev"

        ####################
        # reorientation
        ####################
        record = SeqIO.read(input, "fasta")
        # length of chromosome
        length = len(record.seq)

        # reorient to start at the terminase  
        if strand == "fwd":
            left = record.seq[(int(dnaa_start_on_chromosome) - 1):length]
            right = record.seq[0:int(dnaa_start_on_chromosome)-1]
            total_seq = left + right

        # revese compliment if the strand is negative
        if strand == "rev":
            record.seq = record.seq.reverse_complement()
            left = record.seq[(length - int(dnaa_start_on_chromosome)):length]
            right = record.seq[0:(length - int(dnaa_start_on_chromosome))]
            total_seq = left + right

        # updates the sequence
        record.seq = total_seq

        # writes to file
        with open(out_file, 'w') as out_fa:
            SeqIO.write(record, out_fa, 'fasta')

    else:
        logger.info(f"{gene} start not identified. Please check your input file or try dnaapler custom. If you have assembled an unusual species, this may also cause this error.")
        ctx.exit(2)
    



# function to touch create a file 
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch_file(path):
    with open(path, 'a'):
        os.utime(path, None)

# to create empty fasta file files
def touch_output_files(out_dir, prefix):
    touch_file(os.path.join(out_dir, prefix + "_reoriented.fasta"))

