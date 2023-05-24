
import os
import sys
from Bio import SeqIO
import shutil
from loguru import logger
import click
from pathlib import Path

def instantiate_dirs(output_dir: str, force: bool, ctx:click.Context):
	"""	Checks the output directory
	:param out_dir: output directory path 
	:param force: force flag
	:param logger: logger
	:return: out_dir: final output directory
	"""

	# Checks the output directory
	# remove outdir on force
	logger.info(f"Checking the output directory {output_dir}")
	if force == True:
		if Path(output_dir).exists():
			shutil.rmtree(output_dir)
		else:
			logger.info(f"\n--force was specified even though the output directory does not already exist. Continuing. \n")
	else:
		if Path(output_dir).exists():
			logger.error(f"Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory.")
			ctx.exit(2) 
	# instantiate outdir
	if Path(output_dir).exists() == False:
		Path(output_dir).mkdir(parents=True, exist_ok=True)


def validate_fasta(input_fasta: str, ctx:click.Context ):
	logger.info(f"Checking that the input file {input_fasta} is in FASTA format and has only 1 entry.")
	# to get extension
	with open(input_fasta, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")
		if any(fasta):
				logger.info(f"{input_fasta} file checked.")
		else:
			logger.error(f"Error: {input_fasta} file is not in the FASTA format. Please check your input file")
			ctx.exit(2) 

	with open(input_fasta, "r") as handle:
		# Check the number of records
		if len(list(SeqIO.parse(handle, "fasta"))) == 1:
			logger.info(f"{input_fasta} has only one entry.")
		else:
			logger.error(f"{input_fasta} has more than one entry. Please check your input FASTA file!")
			ctx.exit(2)




