import time
from loguru import logger
import shutil
from Bio import SeqIO
from pathlib import Path
import os
import click
import sys

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
    validate_fasta
)


from .processing import (
    process_blast_output_and_reorient
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







@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


#main_cli.add_command(run)
main_cli.add_command(citation)













# if __name__ == "__main__":





#     print("Running Blast.")
#     logger.info("Running Blast.")

    

#     # run blast
#     absolute_path = os.path.dirname(__file__)
#     processes.blast(args.chromosome, args.outdir, prefix, absolute_path, args.threads, logger)

#     # process output
#     hits = processes.process_output(args.chromosome, args.outdir, prefix, logger) 
    


#     # # extract phage if 2 hits, otherwise create empty file (for snakemake etc)
#     # if hits == 2:
#     #     print("Extracting hlb disrupting sequence.")
#     #     logger.info("Extracting hlb disrupting sequence.")
#     #     processes.extract_prophage(args.chromosome, args.outdir, prefix, logger)
#     # else:
#     #     processes.touch_output_files(args.outdir, prefix)

#     # Determine elapsed time
#     elapsed_time = time.time() - start_time
#     elapsed_time = round(elapsed_time, 2)

#     # Show elapsed time for the process
#     logger.info("dnaapler has finished")
#     logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

#     print("dnaapler has finished")
#     print("Elapsed time: "+str(elapsed_time)+" seconds")

    




