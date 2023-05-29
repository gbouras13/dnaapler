#!/usr/bin/env python3

"""graphbin: Refined binning of metagenomic contigs using assembly graphs."""

import logging
import os
import sys

import click

from graphbin.utils import (
    graphbin_Canu,
    graphbin_Flye,
    graphbin_MEGAHIT,
    graphbin_Miniasm,
    graphbin_SGA,
    graphbin_SPAdes,
)


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.7.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Production"


class ArgsObj:
    def __init__(
        self,
        assembler,
        graph,
        contigs,
        paths,
        binned,
        output,
        prefix,
        max_iteration,
        diff_threshold,
        delimiter,
    ):
        self.assembler = assembler
        self.graph = graph
        self.contigs = contigs
        self.paths = paths
        self.binned = binned
        self.output = output
        self.prefix = prefix
        self.max_iteration = max_iteration
        self.diff_threshold = diff_threshold
        self.delimiter = delimiter


@click.command()
@click.option(
    "--assembler",
    help="name of the assembler used (SPAdes, SGA or MEGAHIT). GraphBin supports Flye, Canu and Miniasm long-read assemblies as well.",
    type=click.Choice(
        ["spades", "sga", "megahit", "flye", "canu", "miniasm"], case_sensitive=False
    ),
    required=True,
)
@click.option(
    "--graph",
    help="path to the assembly graph file",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--contigs",
    help="path to the contigs file",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--paths",
    help="path to the contigs.paths (metaSPAdes) or assembly.info (metaFlye) file",
    type=click.Path(exists=True),
    required=False,
)
@click.option(
    "--binned",
    help="path to the .csv file with the initial binning output from an existing tool",
    type=click.Path(exists=True),
    required=True,
)
@click.option(
    "--output",
    help="path to the output folder",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    required=True,
)
@click.option(
    "--prefix",
    help="prefix for the output file",
    type=str,
    required=False,
)
@click.option(
    "--max_iteration",
    help="maximum number of iterations for label propagation algorithm",
    type=int,
    default=100,
    show_default=True,
    required=False,
)
@click.option(
    "--diff_threshold",
    help="difference threshold for label propagation algorithm",
    type=click.FloatRange(0, 1),
    default=0.1,
    show_default=True,
    required=False,
)
@click.option(
    "--delimiter",
    help="delimiter for input/output results. Supports a comma (,), a semicolon (;), a tab ($'\\t'), a space (\" \") and a pipe (|)",
    type=click.Choice([",", ";", "$'\\t'", '" "'], case_sensitive=False),
    default=",",
    show_default=True,
    required=False,
)
@click.version_option(__version__, "-v", "--version", is_flag=True)
def main(
    assembler,
    graph,
    contigs,
    paths,
    binned,
    output,
    prefix,
    max_iteration,
    diff_threshold,
    delimiter,
):
    """
    GraphBin: Refined Binning of Metagenomic Contigs using Assembly Graphs
    """

    # Setup logger
    # ---------------------------------------------------

    logger = logging.getLogger("GraphBin %s" % __version__)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    fileHandler = logging.FileHandler(f"{output}{prefix}graphbin.log")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    # Validate options
    # ---------------------------------------------------

    # Check if paths files is provided when the assembler type is SPAdes
    if assembler.lower() == "spades" and paths is None:
        logger.error("Please make sure to provide the path to the contigs.paths file.")
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    # Check if paths files is provided when the assembler type is Flye
    if assembler.lower() == "flye" and paths is None:
        logger.error("Please make sure to provide the path to the contigs.paths file.")
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    # Validate prefix
    if prefix != None:
        if not prefix.endswith("_"):
            prefix = prefix + "_"
    else:
        prefix = ""

    # Validate max_iteration
    if max_iteration <= 0:
        logger.error("Please enter a valid number for max_iteration")
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    # Validate diff_threshold
    if diff_threshold < 0:
        logger.error("Please enter a valid number for diff_threshold")
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    # # Remove previous files if they exist
    # if os.path.exists(f"{output}{prefix}graphbin.log"):
    #     os.remove(f"{output}{prefix}graphbin.log")
    # if os.path.exists(f"{output}{prefix}graphbin_output.csv"):
    #     os.remove(f"{output}{prefix}graphbin_output.csv")
    # if os.path.exists(f"{output}{prefix}graphbin_unbinned.csv"):
    #     os.remove(f"{output}{prefix}graphbin_unbinned.csv")

    # Make args object
    args = ArgsObj(
        assembler,
        graph,
        contigs,
        paths,
        binned,
        output,
        prefix,
        max_iteration,
        diff_threshold,
        delimiter,
    )

    # Run GraphBin
    # ---------------------------------------------------
    if assembler.lower() == "canu":
        graphbin_Canu.main(args)
    if assembler.lower() == "flye":
        graphbin_Flye.main(args)
    if assembler.lower() == "megahit":
        graphbin_MEGAHIT.main(args)
    if assembler.lower() == "miniasm":
        graphbin_Miniasm.main(args)
    if assembler.lower() == "sga":
        graphbin_SGA.main(args)
    if assembler.lower() == "spades":
        graphbin_SPAdes.main(args)

    # Exit program
    # --------------

    logger.info("Thank you for using GraphBin! Bye...!")

    logger.removeHandler(fileHandler)
    logger.removeHandler(consoleHeader)


if __name__ == "__main__":
    main()
