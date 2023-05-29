import shutil
import sys
from pathlib import Path

from Bio import SeqIO
from loguru import logger


def instantiate_dirs(output_dir: str, force: bool):
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
    if force == True:
        if Path(output_dir).exists():
            shutil.rmtree(output_dir)
        else:
            logger.info(
                f"\n--force was specified even though the output directory does not already exist. Continuing. \n"
            )
    else:
        if Path(output_dir).exists():
            logger.error(
                f"Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory."
            )

    # instantiate outdir
    if Path(output_dir).exists() == False:
        Path(output_dir).mkdir(parents=True, exist_ok=True)


def validate_fasta(input_fasta: Path):
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


def validate_custom_db_fasta(custom_fasta: Path):
    """
    Validates custom db FASTA input - ensures it is FASTA file with amino acids (.faa)
    """
    logger.info(
        f"Checking that the custom database file {custom_fasta} is in FASTA format with amino acid sequences."
    )

    with open(custom_fasta, "r") as handle:
        # checks if fasta
        fasta = SeqIO.parse(handle, "fasta")
        if any(fasta) == False:
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
