"""MISC FUNCTIONS
You shouldn't need to tweak these much if at all
"""

import os
import shutil
import subprocess as sp
import sys
import time
from pathlib import Path

import click
import pyrodigal
from Bio import SeqIO
from loguru import logger

from dnaapler.utils.cds_methods import run_largest, run_mystery, run_nearest


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
        with open(log, "a") as lo:
            lo.write(msg)


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


def begin_dnaapler(input, output, threads, gene, params):
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
    logger.info("Repository homepage is https://github.com/gbouras13/dnaapler")
    logger.info("Written by George Bouras: george.bouras@adelaide.edu.au")
    logger.info(f"Your input FASTA is {input}")
    logger.info(f"Your output directory  is {output}")
    logger.info(f"You have specified {threads} threads to use with MMseqs2")
    logger.info(f"You have specified {gene} gene(s) to reorient your sequence")
    # check MMseqs2 version
    check_mmseqs2_version()
    check_pyrodigal_version()
    for key, value in params.items():
        logger.info(f"Parameter: {key} {value}")
    return start_time


def check_pyrodigal_version():
    """
    checks the pyrodigal version
    """
    # MMseqs
    message = "Checking pyrodigal installation."
    logger.info(message)

    try:
        pyrodigal_version = pyrodigal.__version__
        pyrodigal_major_version = int(pyrodigal_version.split(".")[0])

        if pyrodigal_major_version < 3:
            logger.error(
                "Pyrodigal is the wrong version. It needs to be v3.0.0 or higher. Please reinstall Dnaapler."
            )

        logger.info(f"Pyrodigal version is v{pyrodigal_version}")
        logger.info(f"Pyrodigal version is ok.")

    except Exception:
        message = "Pyrodigal not found."
        logger.error(message)


def check_mmseqs2_version():
    """
    checks the MMseqs2 version
    """

    # 13.45111 to match pharokka!

    message = "Checking MMseqs2 installation."
    logger.info(message)
    try:
        process = sp.Popen(["mmseqs"], stdout=sp.PIPE, stderr=sp.STDOUT)
        mmseqs_out, _ = process.communicate()
        mmseqs_out = mmseqs_out.decode()
        for line in mmseqs_out.split("\n"):
            if "Version" in line:
                mmseqs_version = line.split(" ")[2]
                break
        else:
            raise ValueError("MMseqs2 version not found")

        mmseqs_major_version = int(mmseqs_version.split(".")[0])
        mmseqs_minor_version = mmseqs_version.split(".")[1]

        logger.info(
            f"MMseqs2 version found is v{mmseqs_major_version}.{mmseqs_minor_version}"
        )

        if mmseqs_major_version != 13:
            logger.error("MMseqs2 is the wrong version. Please install v13.45111")
        if mmseqs_minor_version != "45111":
            logger.error("MMseqs2 is the wrong version. Please install v13.45111")

        logger.info("MMseqs2 version is ok.")

    except Exception:
        message = "MMseqs2 not found. Please install MMseqs2 v 13.45111"
        logger.error(message)


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


def run_autocomplete(
    MMseqs_success, autocomplete, ctx, input, seed_value, output, prefix
):
    """Processes
    :param: MMseqs_success: bool - whether a BLAST hit with a valid start codon was identified
    :param: autocomplete: str - either "none" "mystery" or "nearest"
    :return:
    """

    # if there was
    if MMseqs_success == False:
        if autocomplete == "none":
            logger.error(
                "BLAST based reorientation failed.\n"
                f"Because you chose the {autocomplete} autocomplete strategy as a backup strategy for reorientation, Dnaapler will now exit.\n"
                "Please check your input file, or choose mystery or nearest with the --autocomplete flag."
            )
        elif autocomplete == "mystery":
            logger.info(
                "BLAST based reorientation failed.\n"
                f"You chose the {autocomplete} autocomplete strategy as a backup strategy for reorientation.\n "
                f"Running {autocomplete}."
            )
            run_mystery(ctx, input, seed_value, output, prefix)
        elif autocomplete == "nearest":
            logger.info(
                "BLAST based reorientation failed.\n"
                f"You chose the {autocomplete} autocomplete strategy as a backup strategy for reorientation.\n "
                f"Running {autocomplete}."
            )
            run_nearest(ctx, input, output, prefix)
        elif autocomplete == "largest":
            logger.info(
                "BLAST based reorientation failed.\n"
                f"You chose the {autocomplete} autocomplete strategy as a backup strategy for reorientation.\n "
                f"Running {autocomplete}."
            )
            run_largest(ctx, input, output, prefix)


"""
for all and bulk
"""


def check_duplicate_headers(fasta_file: Path) -> None:
    """
    checks if there are duplicated in the FASTA header
    https://github.com/gbouras13/pharokka/issues/293
    """
    header_set = set()

    # Iterate through the FASTA file and check for duplicate headers
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        if header in header_set:
            logger.error(
                f"Duplicate FASTA header {header} found in the input file {fasta_file}."
            )  # errors if duplicate header found
        else:
            header_set.add(header)
    # if it finished it will be fine


def remove_directory(dir_path: Path) -> None:
    """
    Remove a directory and all its contents if it exists.

    Parameters:
        dir_path (Path): Path to the directory to remove.

    Returns:
        None
    """
    if dir_path.exists():
        shutil.rmtree(dir_path)
