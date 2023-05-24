#!/usr/bin/env python3
# import input_commands
# import processes
# import os
# import logging
import time
from loguru import logger
import time
import shutil
from Bio import SeqIO

from .util import (
    dnaapler_base,
    get_version,
    OrderedCommands,
    print_citation,
)

import os
import click
import sys


"""
code taken from tbpore https://github.com/mbhall88/tbpore
"""

log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)



"""
code adapted from Snaketool https://github.com/beardymcjohnface/Snaketool
"""

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
            '-t,',
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
Main run command
"""

@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
def run(
    ctx,
    input,
    output,
    threads,
    prefix,
    force,
    **kwargs
):
    """Run dnaapler"""

    # get start time
    start_time = time.time()
    # initial logging stuff
    logger.add(os.path.join(output, "dnaapler.log"), rotation="500 MB", level="DEBUG")
    logger.info(f"You are using dnaapler version {get_version()}")
    logger.info(f"Repository homepage is https://github.com/gbouras13/dnaapler")
    logger.info(f"Written by George Bouras: george.bouras@adelaide.edu.au")
    logger.info(f"Your input FASTA specified is {input}")
    logger.info(f"Your output directory specified  is {output}")
    logger.info(f"You have specified {threads} threads")

	# Checks the output directory
	# remove outdir on force
    logger.info(f"Checking the out directory")
    if force == True:
        if os.path.isdir(output) == True:
            shutil.rmtree(output)
        else:
            logger.info(f"\n--force was specified even though the output directory does not already exist. Continuing. \n")
    else:
        if os.path.isdir(output) == True:
            logger.error(f"Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory.")
            ctx.exit(2) 
	# instantiate outdir
    if os.path.isdir(output) == False:
        os.mkdir(output)

###################################
    # checks input FASTA
###################################
    logger.info(f"Checking the input FASTA")
    # to get extension
    with open(input, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta):
            logger.info(f"{input} file checked")
        else:
            logger.error("Error: {input} file is not in the FASTA format. Please check your input file")
            ctx.exit(2) 

    ##### use external_tools.py from tbpore














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

    




