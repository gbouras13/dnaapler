import os
import sys
from Bio import SeqIO
import shutil


def instantiate_dirs(output_dir, force, logger):
	"""
	Checks the output directory
    :param out_dir: output directory path 
    :param force: force flag
    :param logger: logger
    :return: out_dir: final output directory
	"""
	# remove outdir on force
	if force == True:
		if os.path.isdir(output_dir) == True:
			shutil.rmtree(output_dir)
		else:
			print("\n--force was specified even though the output directory does not already exist. Continuing \n")
	else:
		if os.path.isdir(output_dir) == True:
			sys.exit("\nOutput directory exists and force was not specified. Please specify -f or --force to overwrite the output directory. \n")  
	# instantiate outdir
	if os.path.isdir(output_dir) == False:
		os.mkdir(output_dir)
	return output_dir



def validate_fasta(file):
	# to get extension
	with open(file, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")
		if any(fasta):
				print("FASTA " +file + " checked")
		else:
			sys.exit("Error: Input file is not in the FASTA format.\n") 



