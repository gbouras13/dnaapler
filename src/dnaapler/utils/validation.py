import re
import shutil
import subprocess as sp
import sys
from pathlib import Path

import click
from Bio import SeqIO
from loguru import logger


def instantiate_dirs(output_dir: str, force: bool) -> None:
    """Checks the output directory
    :param out_dir: output directory path
    :param force: force flag
    :param logger: logger
    :return: out_dir: final output directory
    """

    # Checks the output directory
    # remove outdir on force
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Checking the output directory {output_dir}")
    if force is True:
        if Path(output_dir).exists():
            shutil.rmtree(output_dir)
        else:
            logger.info(
                "--force was specified even though the output directory does not already exist. Continuing."
            )
    else:
        if Path(output_dir).exists():
            logger.error(
                "Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory."
            )

    # instantiate outdir
    if Path(output_dir).exists() is False:
        Path(output_dir).mkdir(parents=True, exist_ok=True)


def validate_fasta(input_fasta: Path) -> None:
    """
    Validates  FASTA input - that the input is a FASTA with 1 sequence
    """
    logger.info(
        f"Checking that the input file {input_fasta} is in FASTA format and has only 1 entry."
    )
    # to get extension
    with open(input_fasta, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta):
            logger.info(f"{input_fasta} file checked.")
        else:
            logger.error(
                f"Error: {input_fasta} file is not in the FASTA format. Please check your input file"
            )

    with open(input_fasta, "r") as handle:
        # Check the number of records
        if len(list(SeqIO.parse(handle, "fasta"))) == 1:
            logger.info(f"{input_fasta} has only one entry.")
        else:
            logger.error(
                f"{input_fasta} has more than one entry. Please check your input FASTA file!"
            )


def validate_fasta_bulk(input_fasta: Path) -> None:
    """
    Validates  FASTA input - that the input is a FASTA with at least 1 sequence
    """
    logger.info(
        f"Checking that the input file {input_fasta} is in FASTA format and has at least 1 entry."
    )
    # to get extension
    with open(input_fasta, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta):
            logger.info(f"{input_fasta} file checked.")
        else:
            logger.error(
                f"Error: {input_fasta} file is not in the FASTA format. Please check your input file"
            )

    # with open(input_fasta, "r") as handle:
    #     # Check the number of records
    #     if len(list(SeqIO.parse(handle, "fasta"))) == 1:
    #         logger.error(
    #             f"{input_fasta} has only one entry, but more than one was expected. Please check your input FASTA file!"
    #         )
    #     else:
    #         logger.info(f"{input_fasta} has more than one entry.")


def validate_fasta_all(input_fasta: Path) -> None:
    """
    Validates  FASTA input - that the input is a FASTA with >= 1 sequence
    """
    logger.info(
        f"Checking that the input file {input_fasta} is in FASTA format and has at least 1 entry."
    )
    # to get extension
    with open(input_fasta, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta):
            logger.info(f"{input_fasta} file checked.")
        else:
            logger.error(
                f"Error: {input_fasta} file is not in the FASTA format. Please check your input file"
            )

    with open(input_fasta, "r") as handle:
        # Check the number of records
        if len(list(SeqIO.parse(handle, "fasta"))) == 1:
            logger.info(f"{input_fasta} has only one entry.")
        else:
            logger.info(f"{input_fasta} has more than one entry.")


def validate_custom_db_fasta(custom_fasta: Path) -> None:
    """
    Validates custom db FASTA input - ensures it is FASTA file with amino acids (.faa)
    """
    logger.info(
        f"Checking that the custom database file {custom_fasta} is in FASTA format with amino acid sequences."
    )

    with open(custom_fasta, "r") as handle:
        # checks if fasta
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta) is False:
            logger.error(
                f"Error: {custom_fasta} file is not in the FASTA format. Please check your input file"
            )

        # checks amino acids
    with open(custom_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if bool(is_protein_sequence(str(record.seq))) is False:
                logger.error(
                    f"Error: {custom_fasta} file does not seem to have amino acid sequences. Please check your input file"
                )

    logger.info(f"{custom_fasta} file checked.")


def is_protein_sequence(string):
    protein_letters = "acdefghiklmnpqrstvwy"
    nucleotide_letters = "acgnt"

    # Check if the string contains only nucleotide letters
    if all(letter.lower() in nucleotide_letters for letter in string):
        return False

    # Check if the string contains any protein letters
    return any(letter.lower() in protein_letters for letter in string)


def is_scientific_notation(evalue):
    """
    checks if evalue is scientific notation
    """
    # Define the regular expression pattern for scientific notation
    scientific_pattern = r"^[+\-]?(\d+(\.\d*)?|\.\d+)([eE][+\-]?\d+)?$"

    # Check if the number matches the scientific notation pattern
    return bool(re.match(scientific_pattern, evalue))


def is_numeric(evalue):
    """
    checks if evalue is numeric
    """
    try:
        float(evalue)  # Attempt to convert the value to a float
        return True
    except ValueError:
        return False


def check_evalue(evalue):
    """
    checks if the evalue is scientific notation or numeric
    """

    logger.info(f"You have specified an evalue of {evalue}.")

    if is_numeric(evalue) is False and is_scientific_notation(evalue) is False:
        logger.error(
            f"Error: evalue {evalue} is neither numeric nor in scientific notation. Please check your evalue."
        )


def validate_choice_autocomplete(ctx, param, value):
    """
    checks the click.Choice option for the autocomplete flag
    """
    choices = ["mystery", "nearest", "largest", "none"]
    if value not in choices:
        raise click.BadParameter(
            f"Invalid choice. Please choose from {', '.join(choices)}"
        )
    return value


def validate_choice_mode(ctx, param, value):
    """
    checks the click.Choice option for the mode flag in bulk subcommand
    """
    choices = ["chromosome", "phage", "plasmid", "custom", "archaea"]
    if value not in choices:
        raise click.BadParameter(f"Invalid choice. Choose from {', '.join(choices)}")
    return value


def validate_choice_db(ctx, param, value):
    """
    checks the click.Choice option for the mode flag in bulk subcommand
    """
    choices = [
        "all",
        "dnaa",
        "repa",
        "terl",
        "dnaa,repa",
        "dnaa,terl",
        "repa,terl",
        "cog1474",
        "dnaA,cog1474",
        "cog1474,terl",
        "repa,cog1474",
        "dnaa,repa,cog1474",
        "dnaa,terl,cog1474",
        "repa,terl,cog1474",
    ]
    if value not in choices:
        raise click.BadParameter(f"Invalid choice. Choose from {', '.join(choices)}")
    return value


def validate_ignore_file(ignore_file_path):
    try:
        # Open the file in read mode
        with open(ignore_file_path, "r") as file:
            # Read the first character
            first_char = file.read(1)
            # If the first character is not empty, will be true
            return bool(first_char)
    except FileNotFoundError:
        logger.error(f"{ignore_file_path} not found")
