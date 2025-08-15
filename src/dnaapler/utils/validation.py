import re
import shutil
import sys
from pathlib import Path

import click
from Bio import SeqIO
from loguru import logger

from dnaapler.utils.gfa import gfa_sequence_count, is_gfa


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


def is_fasta(input_file):
    """
    Check if the file is in FASTA format.
    """
    with open(input_file, "r") as handle:
        first_char = handle.read(1)
        if first_char != ">":
            return False
    with open(input_file, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)


def check_file_format(input_file):
    """
    Check if the file is in FASTA or GFA format.
    """
    if is_fasta(input_file):
        logger.info(f"{input_file} is in FASTA format.")
    elif is_gfa(input_file):
        logger.info(f"{input_file} is in GFA format.")
    else:
        logger.error(
            f"Error: {input_file} is not in FASTA or GFA format. Please check your input file."
        )


def number_of_sequences(input_file):
    """
    Return the number of sequences in a FASTA or GFA file.
    """
    if is_fasta(input_file):
        with open(input_file, "r") as handle:
            return sum(1 for _ in SeqIO.parse(handle, "fasta"))
    if is_gfa(input_file):
        return gfa_sequence_count(input_file)
    return 0


def validate_input(input_file: Path) -> None:
    """
    Validates input - that the input is a FASTA or GFA with 1 sequence
    """
    logger.info(
        f"Checking that the input file {input_file} is in FASTA or GFA format and has only 1 entry."
    )
    check_file_format(input_file)
    if number_of_sequences(input_file) == 1:
        logger.info(f"{input_file} has only one entry.")
    else:
        logger.error(
            f"{input_file} has more than one entry. Please check your input FASTA file!"
        )


def validate_input_bulk(input_file: Path) -> None:
    """
    Validates input - that the input is a FASTA or GFA with at least 1 sequence
    """
    logger.info(
        f"Checking that the input file {input_file} is in FASTA or GFA format and has at least 1 entry."
    )
    check_file_format(input_file)


def validate_input_all(input_file: Path) -> None:
    """
    Validates input - that the input is a FASTA or GFA with >= 1 sequence
    """
    logger.info(
        f"Checking that the input file {input_file} is in FASTA or GFA format and has at least 1 entry."
    )
    check_file_format(input_file)
    if number_of_sequences(input_file) == 1:
        logger.info(f"{input_file} has only one entry.")
    else:
        logger.info(f"{input_file} has more than one entry.")


def validate_custom_db_fasta(custom_fasta: Path) -> None:
    """
    Validates custom db FASTA input - ensures it is FASTA file with amino acids (.faa)
    """
    logger.info(
        f"Checking that the custom database file {custom_fasta} is in FASTA format with amino acid sequences."
    )
    if not is_fasta(custom_fasta):
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
        "cog1474,repa",
        "dnaa,cog1474,repa",
        "dnaa,cog1474,terl",
        "cog1474,repa,terl",
    ]
    if value not in choices:
        raise click.BadParameter(f"Invalid choice. Choose from {', '.join(choices)}")
    return value


def process_ignore_input(ignore_input):
    """
    Process ignore input which can be either:
    1. A comma separated list of chromosome names
    2. A path to a file containing chromosome names (one per line)
    3. "-" to read from stdin

    Returns a list of chromosome names to ignore
    """

    if not ignore_input or ignore_input == "":
        return []

    if ignore_input == "-":
        try:
            return [x.rstrip().split()[0] for x in sys.stdin if x.strip()]
        except Exception as e:
            logger.error(f"Error reading from stdin: {e}")
            sys.exit(1)

    # Check if it looks like a file path (contains path separators or has an extension)
    path_obj = Path(ignore_input)
    if "/" in ignore_input or "\\" in ignore_input or "." in path_obj.name:
        # It looks like a file path, so check if it exists
        if path_obj.is_file():
            with open(ignore_input) as f:
                return [x.rstrip().split()[0] for x in f if x.strip()]
        else:
            logger.error(f"{ignore_input} not found")
            sys.exit(1)
    else:
        # It's a comma separated list
        chromosomes = []
        parts = ignore_input.split(",")
        for part in parts:
            part = part.strip()
            if part:  # Only add non-empty parts
                chromosomes.append(part)
        return chromosomes
