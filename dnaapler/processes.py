import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from loguru import logger



def blast(input, out_dir, prefix, bin_path, threads, logger):
    blast_db = os.path.join(bin_path, "..", 'db', "dnaA_db") 
    output_file = os.path.join(out_dir, prefix + "_blast_output.txt") 
    try:
        blast = sp.Popen(["blastx", "-query", input, "-db",blast_db, "-evalue",  "1e-20", "-num_threads", str(threads),   "-out", output_file, "-outfmt", "6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq"  ], stdout=sp.PIPE, stderr=sp.PIPE) 
        write_to_log(blast.stderr, logger)
    except:
        sys.exit("Error with blast\n")  


def process_output(input, out_dir, prefix, logger):
    blast_file =  os.path.join(out_dir, prefix + "_blast_output.txt") 
    col_list = ["qseqid", "qlen", "sseqid", "slen", "length", "qstart", "qend", "sstart", "send", "pident", "nident", "gaps", "mismatch", "evalue", "bitscore", "qseq", "sseq"] 
    blast_df = pd.read_csv(blast_file, delimiter= '\t', index_col=False , names=col_list) 

    # take the top hit - blast sorts by bitscore

    if blast_df['qstart'][0] == 1 or blast_df['qend'][0] == 1:
        logger.info("Your input chromosome is already oriented to begin with dnaA.")
        sys.exit("Your input chromosome is already oriented to begin with dnaA.")

    # prokaryotes can use AUG M, GUG V or UUG L as start codons - PA01  dnaA actually starts with V
    # also require the start codon to be the 1st in the searched sequence match - unlikely to not be but nonetheless
    elif blast_df['qseq'][0][0] in ['M', 'V', 'L'] and (blast_df['sstart'][0] == 1 or blast_df['send'][0] == 1):
        print("dnaA identified. Re-orienting chromosome.")
        logger.info("dnaA identified. Re-orienting chromosome.")

        dnaa_start_on_chromosome = blast_df['qstart'][0]

        # 
        strand = "fwd"
        if blast_df['qstart'][0] > blast_df['qend'][0]:
            strand = "rev"

        reorient_chromosome(input, out_dir, prefix,  dnaa_start_on_chromosome,  strand)

    else:
        logger.info("dnaA start not identified. Please check your input file. If you have assembled an unusual species, this may also cause this error.")
        sys.exit('dnaA start not identified. Please check your input file. If you have assembled an unusual species, this may also cause this error.')
    


def reorient_chromosome(input, out_dir, prefix, dnaa_start, strand):

    record = SeqIO.read(input, "fasta")
    # length of chromosome
    length = len(record.seq)

    # reorient to start at the terminase  
    if strand == "fwd":
        left = record.seq[(int(dnaa_start) - 1):length]
        right = record.seq[0:int(dnaa_start)-1]
        total_seq = left + right

    # revese compliment if the strand is negative
    if strand == "rev":
        record.seq = record.seq.reverse_complement()
        left = record.seq[(length - int(dnaa_start)):length]
        right = record.seq[0:(length - int(dnaa_start))]
        total_seq = left + right

    record.seq = total_seq

    outfile = os.path.join(out_dir, prefix + "_reoriented.fasta")  

    with open(outfile, 'w') as out_fa:
        SeqIO.write(record, out_fa, 'fasta')


# function to touch create a file 
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch_file(path):
    with open(path, 'a'):
        os.utime(path, None)

# to create empty fasta file files
def touch_output_files(out_dir, prefix):
    touch_file(os.path.join(out_dir, prefix + "_reoriented.fasta"))

