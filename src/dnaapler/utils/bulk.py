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


def run_bulk_blast(
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
    Run the BLAST part for dnaapler bulk (not processing).

    Args:
        ctx (): The context manager for managing the execution context.
        input (Path): Path to the input data for BLAST.
        output (Path): Directory where BLAST results and logs will be saved.
        prefix (str): Prefix for output files.
        gene (str): Gene name to specify the BLAST database ("dnaA", "repA", "terL", "all", or "custom").
        evalue (float): The E-value threshold for BLAST.
        threads (int): The number of threads to use for BLAST.
        custom_db (Path): Path to a custom database file (amino acid Fasta format) if 'gene' is set to 'custom'.

    Yields:
        bool: True if the BLAST process succeeded, False otherwise.
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

    # for chromosome, plasmid or phage or all
    # runs blast
    if gene != "custom":
        # blast
        logdir = Path(f"{output}/logs")
        blast_output = os.path.join(output, f"{prefix}_blast_output.txt")

        db = os.path.join(DNAAPLER_DB, db_name)
        blast = ExternalTool(
            tool="blastx",
            input=f"-query {input}",
            output=f"-out {blast_output}",
            params=f'-db {db} -evalue  {evalue} -num_threads {threads} -outfmt " 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq "',
            logdir=logdir,
        )

        ExternalTool.run_tool(blast, ctx)

    # if gene == custom
    # maikes db then runs blast
    elif gene == "custom":
        # validates custom fasta input for database
        validate_custom_db_fasta(Path(custom_db))

        # make db
        db_dir = os.path.join(output, "custom_db")
        Path(db_dir).mkdir(parents=True, exist_ok=True)
        custom_db_fasta = os.path.join(db_dir, "custom_db.faa")
        shutil.copy2(custom_db, custom_db_fasta)

        logdir = Path(f"{output}/logs")
        # custom db
        # make custom db
        custom_database = os.path.join(db_dir, "custom_db")
        makeblastdb = ExternalTool(
            tool="makeblastdb",
            input=f"-in {custom_db_fasta}",
            output=f"-out {custom_database}",
            params="-dbtype prot ",
            logdir=logdir,
        )

        # chromosome path
        # blast
        logdir = Path(f"{output}/logs")
        blast_output = os.path.join(output, f"{prefix}_blast_output.txt")
        # dnaA da
        db = os.path.join(db_dir, "custom_db")
        blast = ExternalTool(
            tool="blastx",
            input=f"-query {input}",
            output=f"-out {blast_output}",
            params=f'-db {db} -evalue  {evalue} -num_threads {threads} -outfmt " 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq "',
            logdir=logdir,
        )

        tools_to_run = (makeblastdb, blast)
        ExternalTool.run_tools(tools_to_run, ctx)


def bulk_process_blast_output_and_reorient(
    input: Path, blast_file: Path, output: Path, prefix: str
) -> None:
    """
    Processes the BLAST output and saves reoriented and failed-to-reorient contigs.

    For dnaapler bulk

    Reoriented contigs are saved to 'output/prefix_reoriented.fasta'.
    Failed-to-reorient contigs are saved to 'output/prefix_failed_to_reorient.fasta'.

    Args:
        input (Path): Input file containing nucleotide sequences.
        blast_file (Path): BLAST output file.
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

    # read in the dataframe from BLAST
    try:
        blast_df = pd.read_csv(
            blast_file, delimiter="\t", index_col=False, names=col_list
        )

    except Exception:
        logger.error("There was an issue with parsing the BLAST output file.")

    if isinstance(blast_df, pd.DataFrame) and blast_df.empty:
        logger.error(
            "There were 0 BLAST hits. Please check your input file or try dnaapler custom. If you have assembled an understudied species, this may also be the cause."
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
        # need this to match for BLAST
        short_contig = record.id

        contigs.append(contig)

        # Filter the DataFrame where 'qseqid' matches 'contig'
        filtered_df = blast_df[blast_df["qseqid"] == short_contig]

        length_of_df = len(filtered_df)

        # no hits at all then write to the not-reoriented file
        if length_of_df == 0:
            with open(fail_reoriented_output_file, "a") as out_fa:
                SeqIO.write(record, out_fa, "fasta")

            # no hit save for the output DF
            message = "No_BLAST_hits"
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
            # take the top hit - blast sorts by bitscore
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
