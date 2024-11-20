import os
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from loguru import logger

from dnaapler.utils.constants import DNAAPLER_DB
from dnaapler.utils.external_tools import ExternalTool
from dnaapler.utils.processing import reorient_single_record_bulk
from dnaapler.utils.validation import validate_custom_db_fasta


def run_bulk_MMseqs2(
    ctx,
    input: Path,
    output: Path,
    prefix: str,
    gene: str,
    evalue: float,
    threads: int,
    custom_db: Path,
) -> bool:
    """
    Run the MMseqs2 part for dnaapler bulk (not processing).

    Args:
        ctx (): The context manager for managing the execution context.
        input (Path): Path to the input data for MMseqs2.
        output (Path): Directory where MMseqs2 results and logs will be saved.
        prefix (str): Prefix for output files.
        gene (str): Gene name to specify the MMseqs2 database ("dnaA", "repA", "terL", "all", "cog1474", or "custom").
        evalue (float): The E-value threshold for MMseqs2.
        threads (int): The number of threads to use for MMseqs2.
        custom_db (Path): Path to a custom database file (amino acid Fasta format) if 'gene' is set to 'custom'.

    Yields:
        bool: True if the MMseqs2 process succeeded, False otherwise.
    """

    # sets DB directory based of the gene name

    # defines db name
    db_name = "dnaA_db"
    if gene == "dnaA":
        db_name = "dnaA_db"
    elif gene == "repA":
        db_name = "repA_db"
    elif gene == "terL":
        db_name = "terL_db"
    elif gene == "all":
        db_name = "all_db"
    elif gene == "dnaA,repA":
        db_name = "dnaA_repA_db"
    elif gene == "dnaA,terL":
        db_name = "dnaA_terL_db"
    elif gene == "repA,terL":
        db_name = "repA_terL_db"
    elif gene == "cog1474":
        db_name = "cog1474_db"

    logdir = Path(f"{output}/logs")

    # for chromosome, plasmid or phage or all
    # runs MMseqs2
    logdir = Path(f"{output}/logs")
    MMseqs2_output_tmpdir = Path(f"{output}/tmp_MMseqs2_output")
    MMseqs2_output_file = Path(f"{output}/{prefix}_MMseqs2_output.txt")
    # matches the MMseqs2 ones to make subbing MMseqs2 for BLAST as easy as possible
    MMseqs2_columns = "query,qlen,target,tlen,alnlen,qstart,qend,tstart,tend,fident,nident,gapopen,mismatch,evalue,bits,qaln,taln"

    if gene != "custom":
        db = os.path.join(DNAAPLER_DB, db_name)

    elif gene == "custom":
        # validates custom fasta input for database
        validate_custom_db_fasta(Path(custom_db))

        # make db
        db_dir = os.path.join(output, "custom_db")
        Path(db_dir).mkdir(parents=True, exist_ok=True)
        custom_database = os.path.join(db_dir, "custom_db")

        makeMMseqs2db = ExternalTool(
            tool="mmseqs createdb",
            input=f" {custom_db}",
            output=f" {custom_database}",
            params=f"",
            logdir=logdir,
        )

        ExternalTool.run_tool(makeMMseqs2db, ctx)

        db = custom_database

    MMseqs2 = ExternalTool(
        tool="mmseqs easy-search",
        input=f"{input} {db}",
        output=f"{MMseqs2_output_file}",
        params=f"{MMseqs2_output_tmpdir} --search-type 2  --threads {threads} -e {evalue} --format-output {MMseqs2_columns}",
        logdir=logdir,
    )

    ExternalTool.run_tool(MMseqs2, ctx)
    from dnaapler.utils.util import remove_directory

    remove_directory(MMseqs2_output_tmpdir)


def bulk_process_MMseqs2_output_and_reorient(
    input: Path, MMseqs2_file: Path, output: Path, prefix: str
) -> None:
    """
    Processes the MMseqs2 output and saves reoriented and failed-to-reorient contigs.

    For dnaapler bulk

    Reoriented contigs are saved to 'output/prefix_reoriented.fasta'.
    Failed-to-reorient contigs are saved to 'output/prefix_failed_to_reorient.fasta'.

    Args:
        input (Path): Input file containing nucleotide sequences.
        MMseqs2_file (Path): MMseqs2 output file.
        output (Path): Output directory for saving files.
        prefix (str): Prefix for output files.

    Returns:
        None
    """

    # define colnames
    col_list = [
        "qseqid",
        "qlen",
        "sseqid",
        "slen",
        "length",
        "qstart",
        "qend",
        "sstart",
        "send",
        "pident",
        "nident",
        "gaps",
        "mismatch",
        "evalue",
        "bitscore",
        "qseq",
        "sseq",
    ]

    # read in the dataframe from MMseqs2
    try:
        MMseqs2_df = pd.read_csv(
            MMseqs2_file, delimiter="\t", index_col=False, names=col_list
        )

    except Exception:
        logger.error("There was an issue with parsing the MMseqs2 output file.")

    if isinstance(MMseqs2_df, pd.DataFrame) and MMseqs2_df.empty:
        logger.error(
            "There were 0 MMseqs2 hits. Please check your input file or try dnaapler custom. If you have assembled an understudied species, this may also be the cause."
        )

    # Initialize the list to store the IDs
    contigs = []
    starts = []
    strands = []
    top_hits = []
    top_hit_lengths = []
    covered_lens = []
    coverages = []
    idents = []
    identitys = []

    reoriented_output_file = os.path.join(output, f"{prefix}_reoriented.fasta")
    fail_reoriented_output_file = os.path.join(
        output, f"{prefix}_failed_to_reorient.fasta"
    )

    # Read the FASTA file and extract the IDs
    for record in SeqIO.parse(input, "fasta"):
        contig = record.description
        # need this to match for MMseqs2
        short_contig = record.id

        contigs.append(contig)

        # Filter the DataFrame where 'qseqid' matches 'contig'
        filtered_df = MMseqs2_df[MMseqs2_df["qseqid"] == short_contig]

        length_of_df = len(filtered_df)

        # no hits at all then write to the not-reoriented file
        if length_of_df == 0:
            with open(fail_reoriented_output_file, "a") as out_fa:
                SeqIO.write(record, out_fa, "fasta")

            # no hit save for the output DF
            message = "No_MMseqs2_hits"
            starts.append(message)
            strands.append(message)
            top_hits.append(message)
            top_hit_lengths.append(message)
            covered_lens.append(message)
            coverages.append(message)
            idents.append(message)
            identitys.append(message)

        else:
            # reset the index if you want to re-index the filtered DataFrame
            filtered_df.reset_index(drop=True, inplace=True)

            # top hit has a start of 1 ########
            # take the top hit - MMseqs2 sorts by bitscore
            # if the start is 1 of the top hit
            if filtered_df["qstart"][0] == 1:
                # writes to file
                with open(reoriented_output_file, "a") as out_fa:
                    SeqIO.write(record, out_fa, "fasta")

                # no hit save for the output DF
                message = "Contig_already_reoriented"

                starts.append(message)
                strands.append(message)
                top_hits.append(message)
                top_hit_lengths.append(message)
                covered_lens.append(message)
                coverages.append(message)
                idents.append(message)
                identitys.append(message)

            # top hit
            # prokaryotes can use AUG M, GUG V or UUG L as start codons - e.g. for Pseudomonas aeruginosa PA01  dnaA actually starts with V
            # Therefore, I will require the start codon to be the 1st in the searched sequence match - unlikely to not be but nonetheless
            # Sometimes, the top hit might not be the actual gene of interest (e.g. in phages if the terL is disrupted - like SAOMS1)
            # so in these cases, go down the list and check there is a hit with a legitimate start codon
            # I anticipate this will be rare usually, the first will be it :)
            else:
                if filtered_df["qseq"][0][0] in ["M", "V", "L"] and (
                    filtered_df["sstart"][0] == 1
                ):
                    # starts with valid start codon but needs orientation
                    (
                        start,
                        strand,
                        top_hit,
                        top_hit_length,
                        covered_len,
                        coverage,
                        ident,
                        identity,
                    ) = reorient_single_record_bulk(
                        filtered_df,
                        reoriented_output_file,
                        record,
                        overlapping_orf=False,
                    )

                else:  # top hit does not start with a valid start codon
                    (
                        start,
                        strand,
                        top_hit,
                        top_hit_length,
                        covered_len,
                        coverage,
                        ident,
                        identity,
                    ) = reorient_single_record_bulk(
                        filtered_df,
                        reoriented_output_file,
                        record,
                        overlapping_orf=True,
                    )

                    # because pyrodigal
                    if strand == -1:
                        strand = "reverse"
                    else:
                        strand = "forward"

                # save all the stats
                starts.append(start)
                strands.append(strand)
                top_hits.append(top_hit)
                top_hit_lengths.append(top_hit_length)
                covered_lens.append(covered_len)
                coverages.append(coverage)
                idents.append(ident)
                identitys.append(identity)

    # write the example info to file
    #
    # Create a dictionary from the lists
    data = {
        "Contig": contigs,
        "Start": starts,
        "Strand": strands,
        "Top_Hit": top_hits,
        "Top_Hit_Length": top_hit_lengths,
        "Covered_Length": covered_lens,
        "Coverage": coverages,
        "Identical_AAs": idents,
        "Identity_Percentage": identitys,
    }

    # Convert the dictionary to a DataFrame and save output
    bulk_summary_df = pd.DataFrame(data)
    bulk_summary_df.to_csv(
        os.path.join(output, prefix + "_bulk_reorientation_summary.tsv"),
        sep="\t",
        index=False,
    )
