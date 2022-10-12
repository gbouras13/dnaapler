import argparse
import os
import sys
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import shutil
from version import __version__

v = __version__

### GLOBAL VARIABLES

def get_input():
	parser = argparse.ArgumentParser(description='dnaapler: Orients complete bacterial whole genome chromosome assemblies to begin with dnaa gene.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-c', '--chromosome', action="store", help='Bacterial chromosome assembly file in FASTA format.',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='Directory to write the output to.', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-f', '--force', help="Overwrites the output directory.", action="store_true" )
	parser.add_argument('-p', '--prefix', action="store", help='Prefix for output files. This is not required',  default='Default')
	parser.add_argument('-t', '--threads', action="store", help='Threads for BLAST. Defaults to 4.',  default='4')
	parser.add_argument('-V', '--version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir, force):
	# remove outdir on force
	if force == True:
		if os.path.isdir(output_dir) == True:
			shutil.rmtree(output_dir)
		else:
			print("\n--force was specified even though the outdir does not already exist. Continuing \n")
	else:
		if os.path.isdir(output_dir) == True:
			sys.exit("\nOutput directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory. \n")  
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



