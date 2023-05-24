import time
from loguru import logger
import shutil
from Bio import SeqIO
from pathlib import Path
import os
import click
import sys
import pyrodigal
import random

from .external_tools import ExternalTool

from .constants import (
    DNAAPLER_DB
)

from .util import (
    dnaapler_base,
    get_version,
    OrderedCommands,
    print_citation,
)

from .validation import (
    instantiate_dirs,
    validate_fasta,
    validate_custom_db_fasta
)


from .processing import (
    process_blast_output_and_reorient,
    reorient_sequence_random
)





"""
some code taken/adapted from tbpore https://github.com/mbhall88/tbpore
"""

log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

def begin_dnaapler(input, output, threads, gene):
    """
    begins dnaapler
    returns start time
    """
    # get start time
    start_time = time.time()
    # initial logging stuff
    log_file = os.path.join(output, f"dnaapler_{start_time}.log")
    # adds log file
    logger.add(log_file)
    logger.info(f"You are using dnaapler version {get_version()}")
    logger.info(f"Repository homepage is https://github.com/gbouras13/dnaapler")
    logger.info(f"Written by George Bouras: george.bouras@adelaide.edu.au")
    logger.info(f"Your input FASTA is {input}")
    logger.info(f"Your output directory  is {output}")
    logger.info(f"You have specified {threads} threads to use with blastx")
    logger.info(f"You have specified {gene} gene to reoirent your sequence")
    return start_time

def end_dnaapler(start_time):
    """
    finishes dnaapler
    """

   # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("dnaapler has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
           '-i', "--input",
           help="Path to assembly chromosome/plasmid file in FASTA format",
           type=click.Path(),
           required=True,
        ),
        click.option(
            '-o', "--output",
            default="output.dnaapler",
            show_default=True,
            type=click.Path(),
            help="Output directory ",
        ),
        click.option(
            '-t',
            "--threads", 
            help="Number of threads to use with BLAST.",
            default=1, 
            show_default=True
        ),
        click.option(
            '-p',
            "--prefix",
            default="dnaapler",
            help="Prefix for output files. Not required.",
            show_default=True,
        ),
        click.option(
            "-f",
            "--force",
            is_flag=True,
            help="Force overwrites the output directory.",
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
def main_cli():
    logger.info(f"Starting dnaapler")



"""
Chromosome command
"""

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
def chromosome(
    ctx,
    input,
    output,
    threads,
    prefix,
    force,
    **kwargs
):
    """Reorients your sequence to begin with the dnaA chromosomal replication initiation gene"""

    ### validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force, ctx)

    # defines gene
    gene = "dnaA"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta 
    validate_fasta(input, ctx)
    ##### use external_tools.py 

    # chromosome path
    # blast
    logdir = Path(f"{output}/logs")
    blast_output = os.path.join(output, f"{prefix}_blast_output.txt")
    # dnaA da
    db = os.path.join(DNAAPLER_DB, f"dnaA_db")
    blast = ExternalTool(
        tool="blastx",
        input=f"-query {input}",
        output=f"-out {blast_output}",
        params=f"-db {db} -evalue  1e-10 -num_threads {threads} -outfmt \" 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq \"",
        logdir=logdir,
    )

    # tools_to_run = (
    #     blast,
    #     blast
    # )
    #ExternalTool.run_tools(tools_to_run, ctx)
    #  makeblastdb -in dnaA.fasta -dbtype prot -out dnaA_db
    #  makeblastdb -in terL.fasta -dbtype prot -out terL_db
    #  makeblastdb -in repA.fasta -dbtype prot -out repA_db

    # think I onle need 
    # phr
    # pin 
    # psq

    ExternalTool.run_tool(blast, ctx)

    # reorient the 
    output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
    process_blast_output_and_reorient(input, blast_output,output_processed_file,  ctx, gene)

    # end dnaapler
    end_dnaapler(start_time)


"""
Plasmid command
"""

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
def plasmid(
    ctx,
    input,
    output,
    threads,
    prefix,
    force,
    **kwargs
):
    """Reorients your sequence to begin with the repA replication initiation gene"""

    ### validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force, ctx)

    # defines gene
    gene = "repA"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta 
    validate_fasta(input, ctx)
    ##### use external_tools.py 

    # chromosome path
    # blast
    logdir = Path(f"{output}/logs")
    blast_output = os.path.join(output, f"{prefix}_blast_output.txt")
    # dnaA da
    db = os.path.join(DNAAPLER_DB, f"repA_db")
    blast = ExternalTool(
        tool="blastx",
        input=f"-query {input}",
        output=f"-out {blast_output}",
        params=f"-db {db} -evalue  1e-10 -num_threads {threads} -outfmt \" 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq \"",
        logdir=logdir,
    )

    ExternalTool.run_tool(blast, ctx)

    # reorient the 
    output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
    process_blast_output_and_reorient(input, blast_output,output_processed_file,  ctx, gene)

    # end dnaapler
    end_dnaapler(start_time)

"""
Phage command
"""

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
def phage(
    ctx,
    input,
    output,
    threads,
    prefix,
    force,
    **kwargs
):
    """Reorients your sequence to begin with the terL large terminase subunit"""

    ### validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force, ctx)

    # defines gene
    gene = "terL"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta 
    validate_fasta(input, ctx)
    ##### use external_tools.py 

    # chromosome path
    # blast
    logdir = Path(f"{output}/logs")
    blast_output = os.path.join(output, f"{prefix}_blast_output.txt")
    # dnaA da
    db = os.path.join(DNAAPLER_DB, f"terL_db")
    blast = ExternalTool(
        tool="blastx",
        input=f"-query {input}",
        output=f"-out {blast_output}",
        params=f"-db {db} -evalue  1e-10 -num_threads {threads} -outfmt \" 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq \"",
        logdir=logdir,
    )

    ExternalTool.run_tool(blast, ctx)

    # reorient the 
    output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
    process_blast_output_and_reorient(input, blast_output,output_processed_file,  ctx, gene)

    # end dnaapler
    end_dnaapler(start_time)



"""
custom command
"""

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-c",
    "--custom_db",
    help="FASTA file with amino acids that will be used as a custom blast database to reorient your sequence however you want.",
    type=click.Path(),
    required=True
)
def custom(
    ctx,
    input,
    output,
    threads,
    prefix,
    force,
    custom_db,
    **kwargs
):
    """Reorients your sequence with a custom database"""

    ### validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force, ctx)

    # defines gene
    gene = "custom"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta 
    validate_fasta(input, ctx)

    # validates custom fasta input for database 
    validate_custom_db_fasta(custom_db, ctx)



    ##### use external_tools.py 

    # make db

    db_dir = os.path.join(output, f"custom_db")
    Path(db_dir).mkdir(parents=True, exist_ok=True)
    custom_db_fasta = os.path.join(db_dir, "custom_db.faa" )
    shutil.copy2(custom_db, custom_db_fasta)

    logdir = Path(f"{output}/logs")
    # cstom db
    # make custom db
    custom_database = os.path.join(db_dir, "custom_db" )
    makeblastdb = ExternalTool(
        tool="makeblastdb",
        input=f"-in {custom_db_fasta}",
        output=f"-out {custom_database}",
        params=f"-dbtype prot ",
        logdir=logdir
    )

    # chromosome path
    # blast
    logdir = Path(f"{output}/logs")
    blast_output = os.path.join(output, f"{prefix}_blast_output.txt")
    # dnaA da
    db = os.path.join(db_dir, f"custom_db")
    blast = ExternalTool(
        tool="blastx",
        input=f"-query {input}",
        output=f"-out {blast_output}",
        params=f"-db {db} -evalue  1e-10 -num_threads {threads} -outfmt \" 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq \"",
        logdir=logdir,
    )

    tools_to_run = (
        makeblastdb,
        blast
    )
    ExternalTool.run_tools(tools_to_run, ctx)

    # think I only need 
    # phr
    # pin 
    # psq

    # reorient the 
    output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
    process_blast_output_and_reorient(input, blast_output,output_processed_file,  ctx, gene)

    # end dnaapler
    end_dnaapler(start_time)


"""
mystery command
"""

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "--seed_value",
    help="Seed to ensure reproducibility.",
    type=int,
    default=13, 
    show_default=True
)
def mystery(
    ctx,
    input,
    output,
    threads,
    prefix,
    seed_value,
    force,
    **kwargs
):
    """Reorients your sequence with a mystery gene"""

    ### validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force, ctx)

    # defines gene
    gene = "mystery"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta 
    validate_fasta(input, ctx)

    logger.info(f"Searching for genes with pyrodigal")

    # get number of records of input 
    orf_finder = pyrodigal.OrfFinder(meta=True)

    # set seed
    random.seed(int(seed_value))

    # there will only be 1 record
    for i, record in enumerate(SeqIO.parse(input, 'fasta' )):
        genes = orf_finder.find_genes(str(record.seq))
        # get number of genes
        gene_count = len(genes)

         # ensure has > 3 genes
        if gene_count < 4:
            logger.error(f"{input} has less than 4 genes. You probably shouldn't be using dnaapler random!")
            ctx.exit(2)

        logger.info(f"Reorienting with a random gene (that is not the first or last).")

        # ensure not first or last gene
        reorient_gene_number = random.randint(2, gene_count-1)

        logger.info(f"Gene number {reorient_gene_number} was selected.")
        start = genes[reorient_gene_number].begin
        strand = genes[reorient_gene_number].strand

        if strand == 1:
            strand_eng = "forward"
        else:
            strand_eng = "negative"

        logger.info(f"Your random gene has a start coordinate of {start}.")
        logger.info(f"Your random gene is on the {strand_eng} strand.")

        output_processed_file = os.path.join(output, f"{prefix}_reoriented.fasta")
        reorient_sequence_random(input, output_processed_file,  start, strand)

    # finish dnaapler
    end_dnaapler(start_time)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


#main_cli.add_command(run)
main_cli.add_command(citation)


