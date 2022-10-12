#!/usr/bin/env python3
import input_commands
import processes
import os
import logging
import time
import datetime



if __name__ == "__main__":

    # get start time
    start_time = time.time()

    # getting time for log file 

    time_for_log = datetime.datetime.now().strftime("%m%d%Y_%H%M%S")

    args = input_commands.get_input()
        # set the prefix
    if args.prefix == "Default":
        prefix = "dnaapler"
    else:
        prefix = args.prefix
    
    out_dir = input_commands.instantiate_dirs(args.outdir, args.force) # incase there is already an outdir

    LOG_FILE = os.path.join(args.outdir, prefix + "_" + str(time_for_log) + ".log")
    logger = logging.getLogger()
    logging.basicConfig(level=logging.INFO,filename=LOG_FILE,format='%(asctime)s - %(levelname)s - %(message)s')
    print("Starting dnaapler.")
    logger.info("Starting dnaapler")
    print("Checking input fasta.")
    logger.info("Checking input fasta.")

    # instantiation/checking fastq 
    input_commands.validate_fasta(args.chromosome)

    print("Running Blast.")
    logger.info("Running Blast.")

    

    # run blast
    absolute_path = os.path.dirname(__file__)
    processes.blast(args.chromosome, args.outdir, prefix, absolute_path, args.threads, logger)

    # process output
    hits = processes.process_output(args.chromosome, args.outdir, prefix, logger) 
    


    # # extract phage if 2 hits, otherwise create empty file (for snakemake etc)
    # if hits == 2:
    #     print("Extracting hlb disrupting sequence.")
    #     logger.info("Extracting hlb disrupting sequence.")
    #     processes.extract_prophage(args.chromosome, args.outdir, prefix, logger)
    # else:
    #     processes.touch_output_files(args.outdir, prefix)

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("dnaapler has finished")
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")

    print("dnaapler has finished")
    print("Elapsed time: "+str(elapsed_time)+" seconds")

    




