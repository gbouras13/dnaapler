"""MISC FUNCTIONS
You shouldn't need to tweak these much if at all
"""

import os
import subprocess as sp
import sys
import time

import click
from loguru import logger


class OrderedCommands(click.Group):
    """This class will preserve the order of subcommands, which is useful when printing --help"""

    def list_commands(self, ctx: click.Context):
        return list(self.commands)


def dnaapler_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(dnaapler_base("VERSION"), "r") as f:
        version = f.readline()
    return version


def echo_click(msg, log=None):
    click.echo(msg, nl=False, err=True)
    if log:
        with open(log, "a") as l:
            l.write(msg)


def print_citation():
    with open(dnaapler_base("CITATION"), "r") as f:
        for line in f:
            echo_click(line)


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)

"""
begin and end functions
"""


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
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"You are using dnaapler version {get_version()}")
    logger.info(f"Repository homepage is https://github.com/gbouras13/dnaapler")
    logger.info(f"Written by George Bouras: george.bouras@adelaide.edu.au")
    logger.info(f"Your input FASTA is {input}")
    logger.info(f"Your output directory  is {output}")
    logger.info(f"You have specified {threads} threads to use with blastx")
    logger.info(f"You have specified {gene} gene to reoirent your sequence")
    blast_version = check_blast_version()
    logger.info(f"blastx version {blast_version} found")
    return start_time


def check_blast_version():
    # blast
    try:
        process = sp.Popen(["blastx", "-version"], stdout=sp.PIPE, stderr=sp.STDOUT)
        blast_out, _ = process.communicate()
        blast_out = blast_out.decode().strip()
        blast_out = blast_out.split("\n")[0]
        blast_version = blast_out.split(" ")[1]
    except:
        logger.error("blastx not found. Please install blast")

    return blast_version


def end_dnaapler(start_time):
    """
    finishes dnaapler
    """

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("dnaapler has finished")
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")
