import os
from pathlib import Path

import pandas as pd
import pyrodigal
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger


def process_MMseqs2_output_and_reorient(
    input: Path, MMseqs2_file: Path, out_file: Path, gene: str
) -> bool:
    """
    Processes MMseqs2 output and reorients the input sequence if necessary based on the top MMseqs2 hit.

    Args:
        input (Path): Input file containing a nucleotide sequence.
        MMseqs2_file (Path): MMseqs2 output file.
        out_file (Path): Output file to save the reoriented sequence.
        gene (str): The type of gene (e.g., "dnaA").

    Returns:
        bool: Whether a MMseqs2 hit with a valid start codon was identified.
    """

    #  instantiate the returned bool
    MMseqs2_success = True

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
        logger.info(
            "There were 0 MMseqs2 hits. Please check your input file or try dnaapler custom. If you have assembled an understudied species, this may also be the cause."
        )
        MMseqs2_success = False

    # top hit has a start of 1 ########
    # take the top hit - MMseqs2 sorts by bitscore
    # if the start is 1 of the top hit
    elif MMseqs2_df["qstart"][0] == 1:
        logger.info(
            f"Based on the MMseqs2 output top hit, your input is already oriented to begin with {gene}.\n"
            "The input file will be copied to the output so that your pipelines don't break :)"
        )
        # writes to file
        record = SeqIO.read(input, "fasta")
        with open(out_file, "w") as out_fa:
            SeqIO.write(record, out_fa, "fasta")

        MMseqs2_success == True

    # top hit
    # prokaryotes can use AUG M, GUG V or UUG L as start codons - e.g. for Pseudomonas aeruginosa PA01  dnaA actually starts with V
    # Therefore, I will require the start codon to be the 1st in the searched sequence match - unlikely to not be but nonetheless
    # Sometimes, the top hit might not be the actual gene of interest (e.g. in phages if the terL is disrupted - like SAOMS1)
    # so in these cases, go down the list and check there is a hit with a legitimate start codon
    # I anticipate this will be rare usually, the first will be it :)

    # After this issue (and consider pyrodigal is already in this program)
    # I have decided to implement reorienting by the CDS closest to the top hit if it doesn't start with a valid start codon
    # https://github.com/gbouras13/dnaapler/issues/44

    else:
        if MMseqs2_df["qseq"][0][0] in ["M", "V", "L"] and (
            MMseqs2_df["sstart"][0] == 1
        ):
            reorient_sequence(MMseqs2_df, input, out_file, gene, overlapping_orf=False)
        else:  # this will reorient the sequence with the orf that overlaps the tophit by the most
            # warn if the top hit doesnt begin with a valid start codon
            logger.warning(
                f"The top {gene} MMseqs2 hit did not have a MMseqs2 alignment beginning with a valid start codon M, V or L."
            )
            logger.warning(
                "The CDS with the most overlap with this MMseqs2 hit according to Pyrodigal will instead be used to reorient the genome."
            )
            reorient_sequence(MMseqs2_df, input, out_file, gene, overlapping_orf=True)

        MMseqs2_success = True

    return MMseqs2_success


def reorient_sequence(
    MMseqs2_df: pd.DataFrame,
    input: Path,
    out_file: Path,
    gene: str,
    overlapping_orf: bool,
) -> None:
    """
    Reorients the input sequence based on MMseqs2 results and a gene of interest.

    Args:
        MMseqs2_df (pd.DataFrame): DataFrame containing MMseqs2 results.
        input (Path): Input file containing a nucleotide sequence.
        out_file (Path): Output file to save the reoriented sequence.
        gene (str): The type of gene (e.g., "dnaA").
        overlapping_orf (bool): Indicates whether the top MMseqs2 hit overlaps with an existing open reading frame.

    Returns:
        None
    """
    # get the start

    ident = MMseqs2_df["nident"][0]
    top_hit = MMseqs2_df["sseqid"][0]
    top_hit_length = MMseqs2_df["slen"][0]
    covered_len = MMseqs2_df["length"][0]
    identity = round(ident / covered_len * 100, 2)
    coverage = round(covered_len / top_hit_length * 100, 2)

    # top MMseqs2 hit starts at 1
    if overlapping_orf is False:
        # basic stats
        dnaa_start_on_chromosome = MMseqs2_df["qstart"][0]
        strand = "forward"
        if MMseqs2_df["qstart"][0] > MMseqs2_df["qend"][0]:
            strand = "reverse"

        logger.info(
            f"{gene} gene identified. It starts at coordinate {dnaa_start_on_chromosome} on the {strand} strand in your input file."
        )
        logger.info(
            f"The best hit with a valid start codon in the database is {top_hit}, which has length of {top_hit_length} AAs."
        )
        logger.info(
            f"{covered_len} AAs were covered by the best hit, with an overall coverage of {coverage}%."
        )
        logger.info(
            f"{ident} AAs were identical, with an overall identity of {identity}%."
        )
        logger.info("Re-orienting.")

        ####################
        # reorientation
        ####################
        record = SeqIO.read(input, "fasta")
        # length of chromosome
        length = len(record.seq)

        # reorient to start at the terminase
        if strand == "forward":
            left = record.seq[(int(dnaa_start_on_chromosome) - 1) : length]
            right = record.seq[0 : int(dnaa_start_on_chromosome) - 1]
            total_seq = left + right

        # revese compliment if the strand is negative
        if strand == "reverse":
            record.seq = record.seq.reverse_complement()
            left = record.seq[(length - int(dnaa_start_on_chromosome)) : length]
            right = record.seq[0 : (length - int(dnaa_start_on_chromosome))]
            total_seq = left + right

        # updates the sequence
        record.seq = total_seq

        # writes to file
        with open(out_file, "w") as out_fa:
            SeqIO.write(record, out_fa, "fasta")

    # top MMseqs2 hit does not start at 1
    else:
        # get start and end of tophit
        start_tophit = MMseqs2_df["qstart"][0]
        end_tophit = MMseqs2_df["qend"][0]

        tophit_min = min(start_tophit, end_tophit)
        tophit_max = max(start_tophit, end_tophit)

        # there will only be 1 record
        for i, record in enumerate(SeqIO.parse(input, "fasta")):
            logger.info(
                f"The top MMseqs2 hit for {gene} is {top_hit}, which has length of {top_hit_length} AAs.  The alignment did not begin with a valid start codon."
            )
            logger.info(
                f"{covered_len} AAs were covered in the best hit, with an overall coverage of {coverage}%."
            )
            logger.info(
                f"{ident} AAs were identical, with an overall identity of {identity}%."
            )

            orf_finder = pyrodigal.GeneFinder(meta=True)
            genes = orf_finder.find_genes(str(record.seq))

            # dictionary to store distances to the start for each CDS
            overlap_dict = {}

            gene_index = 0

            for gene in genes:
                overlap_start = max(tophit_min, min(gene.begin, gene.end))
                overlap_end = min(tophit_max, max(gene.begin, gene.end))

                overlap_dict[gene_index] = max(
                    0, overlap_end - overlap_start
                )  # Ensure non-negative length
                gene_index += 1

            # Find the gene with the max overlap
            closest_gene_index = max(overlap_dict, key=lambda key: overlap_dict[key])

            # get strand
            strand = genes[closest_gene_index].strand

            # susie error 30-01-24 - misorienting on the negative strand
            # 'begin' just gives the lowest value, not the start, so was putting the terL at the end i.e. reorienting from the end of terL
            # therefore need to take end
            if strand == 1:
                start = genes[closest_gene_index].begin
            elif strand == -1:
                start = genes[closest_gene_index].end

            start = genes[closest_gene_index].begin

            if strand == 1:
                strand_eng = "forward"
            else:
                strand_eng = "negative"

            logger.info(
                f"The CDS most overlapping the tophit has a start coordinate of {start}."
            )
            logger.info(
                f"The CDS most overlapping the tophit is on the {strand_eng} strand."
            )
            logger.info("Re-orienting.")

            reorient_sequence_random(input, out_file, start, strand)


def reorient_sequence_random(
    input: Path, out_file: Path, start: int, strand: int
) -> None:
    """
    Reorients the input sequence based on the provided start coordinate and strand information.

    For mystery nearest and largest

    Args:
        input (Path): Input file containing a nucleotide sequence.
        out_file (Path): Output file to save the reoriented sequence.
        start (int): The start coordinate for reorientation.
        strand (int): The strand information for reorientation (1 for forward, -1 for reverse).

    Returns:
        None
    """
    ####################
    # reorientation
    ####################
    record = SeqIO.read(input, "fasta")
    # length of chromosome
    length = len(record.seq)

    # reorient to start at the terminase
    # forward is 1 from pyrodigal
    if strand == 1:
        left = record.seq[(int(start) - 1) : length]
        right = record.seq[0 : int(start) - 1]
        total_seq = left + right

    # revese compliment if the strand is -1 from pyrodigal
    if strand == -1:
        record.seq = record.seq.reverse_complement()
        left = record.seq[(length - int(start)) : length]
        right = record.seq[0 : (length - int(start))]
        total_seq = left + right

    # updates the sequence
    record.seq = total_seq

    # writes to file
    with open(out_file, "w") as out_fa:
        SeqIO.write(record, out_fa, "fasta")


def reorient_single_record_bulk(
    MMseqs2_df: pd.DataFrame,
    out_file: Path,
    record: SeqIO.SeqRecord,
    overlapping_orf: bool,
) -> tuple:
    """
    Reorients a single DNA sequence record in dnaapler bulk or dnaapler all based on MMseqs2 results.

    Args:
        MMseqs2_df (pd.DataFrame): DataFrame containing MMseqs2 results for the record.
        out_file (Path): Output file to save the reoriented sequence.
        record (SeqIO.SeqRecord): Sequence record to be reoriented.
        overlapping_orf (bool): Indicates whether the top MMseqs2 hit starts with a valid start codon (False) or not (True).

    Returns:
        Tuple of reorientation details:
        - start (int): Start coordinate after reorientation.
        - strand (str): Reoriented strand ('forward' or 'reverse').
        - top_hit (str): Identifier of the top MMseqs2 hit.
        - top_hit_length (int): Length of the top MMseqs2 hit.
        - covered_len (int): Length of the sequence covered by the top MMseqs2 hit.
        - coverage (float): Percentage coverage of the sequence.
        - ident (int): Number of identical amino acids.
        - identity (float): Percentage identity.
    """

    top_hit = MMseqs2_df["sseqid"][0]
    top_hit_length = MMseqs2_df["slen"][0]
    covered_len = MMseqs2_df["length"][0]
    coverage = round(covered_len / top_hit_length * 100, 2)
    ident = MMseqs2_df["nident"][0]
    identity = round(ident / covered_len * 100, 2)

    # default - when tophit starts with a valid start codon
    if overlapping_orf is False:
        # get the start
        dnaa_start_on_chromosome = MMseqs2_df["qstart"][0]
        strand = "forward"
        if MMseqs2_df["qstart"][0] > MMseqs2_df["qend"][0]:
            strand = "reverse"

        ####################
        # reorientation
        ####################

        # length of chromosome
        length = len(record.seq)

        # reorient to start at the terminase
        if strand == "forward":
            left = record.seq[(int(dnaa_start_on_chromosome) - 1) : length]
            right = record.seq[0 : int(dnaa_start_on_chromosome) - 1]
            total_seq = left + right

        # revese compliment if the strand is negative
        if strand == "reverse":
            record.seq = record.seq.reverse_complement()
            left = record.seq[(length - int(dnaa_start_on_chromosome)) : length]
            right = record.seq[0 : (length - int(dnaa_start_on_chromosome))]
            total_seq = left + right

        # updates the sequence
        record.seq = total_seq

        # writes to file
        with open(out_file, "a") as out_fa:
            SeqIO.write(record, out_fa, "fasta")

        return (
            dnaa_start_on_chromosome,
            strand,
            top_hit,
            top_hit_length,
            covered_len,
            coverage,
            ident,
            identity,
        )
    # tophit does not start with 1
    else:
        # get start and end of tophit
        start_tophit = MMseqs2_df["qstart"][0]
        end_tophit = MMseqs2_df["qend"][0]

        tophit_min = min(start_tophit, end_tophit)
        tophit_max = max(start_tophit, end_tophit)

        # find the genes
        logger.warning(
            f"The top MMseqs2 hit for the contig {record.id} did not begin with a valid start codon."
        )
        logger.warning(
            "Searching with pyrodigal for the CDS overlapping the most with the top MMseqs2 hit to reorient with."
        )

        orf_finder = pyrodigal.GeneFinder(meta=True)
        genes = orf_finder.find_genes(str(record.seq))

        # dictionary to store distances to the start for each CDS
        overlap_dict = {}

        gene_index = 0

        for gene in genes:
            overlap_start = max(tophit_min, min(gene.begin, gene.end))
            overlap_end = min(tophit_max, max(gene.begin, gene.end))

            overlap_dict[gene_index] = max(
                0, overlap_end - overlap_start
            )  # Ensure non-negative length
            gene_index += 1

        # Find the gene with the max overlap
        # issue #77 - where there are no overlaps
        if len(overlap_dict) == 0:
            logger.warning(
                f"{record.id} has no CDS. Therefore, skipping reorientation (you should probably filter this contig out!)"
            )
            start = 1
            strand = "forward"
            top_hit = "MMseqs2_hit_but_no_overlapping_CDS"
            top_hit_length = "MMseqs2_hit_but_no_overlapping_CDS"
            covered_len = "MMseqs2_hit_but_no_overlapping_CDS"
            coverage = "MMseqs2_hit_but_no_overlapping_CDS"
            ident = "MMseqs2_hit_but_no_overlapping_CDS"
            identity = "MMseqs2_hit_but_no_overlapping_CDS"

        # 99.99% of scenarios
        else:
            closest_gene_index = max(overlap_dict, key=lambda key: overlap_dict[key])

            strand = genes[closest_gene_index].strand

            # susie error 30-01-24 - misorienting on the negative strand
            # 'begin' just gives the lowest value, not the start, so was putting the terL at the end i.e. reorienting from the end of terL
            # therefore need to take end
            if strand == 1:
                start = genes[closest_gene_index].begin
            elif strand == -1:
                start = genes[closest_gene_index].end

            ####################
            # reorientation
            ####################

            # length of chromosome
            length = len(record.seq)

            # reorient to start at the terminase
            if strand == 1:
                left = record.seq[(int(start) - 1) : length]
                right = record.seq[0 : int(start) - 1]
                total_seq = left + right

            # revese compliment if the strand is negative
            if strand == -1:
                record.seq = record.seq.reverse_complement()
                left = record.seq[(length - int(start)) : length]
                right = record.seq[0 : (length - int(start))]
                total_seq = left + right

            # updates the sequence
            record.seq = total_seq

        # writes to file
        with open(out_file, "a") as out_fa:
            SeqIO.write(record, out_fa, "fasta")

        # update strand in the all output if fallback

        if strand == 1:
            strand = "forward"
        else:
            strand = "reverse"

        return (
            start,
            strand,
            top_hit,
            top_hit_length,
            covered_len,
            coverage,
            ident,
            identity,
        )


def reorient_sequence_and_append(
    record: SeqIO.SeqRecord, out_file: Path, start: int, strand: int
) -> None:
    """
    Reorients a DNA sequence record and appends it to an output file.

    Args:
        record (SeqIO.SeqRecord): Sequence record to be reoriented.
        out_file (Path): Output file to save the reoriented sequence.
        start (int): Start coordinate after reorientation.
        strand (int): Strand orientation (1 for forward, -1 for reverse).

    Returns:
        None
    """

    ####################
    # reorientation
    ####################

    # length of record
    length = len(record.seq)

    # reorient to start at the terminase
    # forward is 1 from pyrodigal
    if strand == 1:
        left = record.seq[(int(start) - 1) : length]
        right = record.seq[0 : int(start) - 1]
        total_seq = left + right

    # revese compliment if the strand is -1 from pyrodigal
    if strand == -1:
        record.seq = record.seq.reverse_complement()
        left = record.seq[(length - int(start)) : length]
        right = record.seq[0 : (length - int(start))]
        total_seq = left + right

    # updates the sequence
    record.seq = total_seq
    # update that it has been rotated
    record.description = record.description + " rotated=True"

    # writes to file
    with open(out_file, "a") as out_fa:
        SeqIO.write(record, out_fa, "fasta")


# function to touch create a file
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)


# to create empty fasta file files
def touch_output_files(out_dir, prefix):
    touch_file(os.path.join(out_dir, prefix + "_reoriented.fasta"))
