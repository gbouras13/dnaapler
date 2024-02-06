import os
import random
from pathlib import Path
from typing import Tuple

import pandas as pd
import pyrodigal
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger

from dnaapler.utils.processing import (
    reorient_sequence_and_append,
    reorient_single_record_bulk,
)


def all_process_blast_output_and_reorient(
    input: Path,
    blast_file: Path,
    output: Path,
    prefix: str,
    ignore_list: Path,
    autocomplete: str,
    seed_value: int,
    custom_db: str,
) -> None:
    """Processes the blast output,reorients and saves all contigs into os.path.join(output, f"{prefix}_reoriented.fasta")

    :param input: input file
    :param blast_file: blast output file
    :param output: output directory
    :param prefix: prefix
    :param ignore_list: list containing contigs (short_contig) to ignore
    :return:
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
    genes = []  # for the gene it reoriented it by
    starts = []
    strands = []
    top_hits = []
    top_hit_lengths = []
    covered_lens = []
    coverages = []
    idents = []
    identitys = []

    reoriented_output_file = os.path.join(output, f"{prefix}_reoriented.fasta")

    # Read the FASTA file and extract the IDs
    for record in SeqIO.parse(input, "fasta"):
        contig = record.description
        # need this to match for BLAST
        short_contig = record.id
        contigs.append(contig)

        # Filter the DataFrame where 'qseqid' matches 'contig'
        filtered_df = blast_df[blast_df["qseqid"] == short_contig]

        length_of_df = len(filtered_df)

        if short_contig in ignore_list:
            # write contig anyway
            with open(reoriented_output_file, "a") as out_fa:
                SeqIO.write(record, out_fa, "fasta")

            # no hit save for the output DF
            message = "Contig_ignored"
            genes.append(message)
            starts.append(message)
            strands.append(message)
            top_hits.append(message)
            top_hit_lengths.append(message)
            covered_lens.append(message)
            coverages.append(message)
            idents.append(message)
            identitys.append(message)

        # no hits at all then autcomplete
        elif length_of_df == 0:
            if autocomplete == "none":
                # write contig anyway
                with open(reoriented_output_file, "a") as out_fa:
                    SeqIO.write(record, out_fa, "fasta")

                # no hit save for the output DF
                message = "No_BLAST_hits"
                genes.append(message)
                starts.append(message)
                strands.append(message)
                top_hits.append(message)
                top_hit_lengths.append(message)
                covered_lens.append(message)
                coverages.append(message)
                idents.append(message)
                identitys.append(message)

            else:
                message = f"autocomplete_method_{autocomplete}"
                (start, strand) = run_autocomplete_record(
                    record, autocomplete, reoriented_output_file, seed_value
                )

                if strand == 1:
                    strand_eng = "forward"
                else:
                    strand_eng = "negative"

                genes.append(message)
                starts.append(start)
                strands.append(strand_eng)
                top_hits.append(message)
                top_hit_lengths.append(message)
                covered_lens.append(message)
                coverages.append(message)
                idents.append(message)
                identitys.append(message)

        else:  # there is at least one BLAST hit
            # determine the numbers of repA, dnaA and terL

            # UniRef90 is in string for repA
            # phrog is in string for terL
            # DNAA is in string for dnaA

            # counts
            counts = {
                "dnaA": filtered_df[
                    filtered_df["sseqid"].str.contains("DNAA", case=False)
                ].shape[0],
                "terL": filtered_df[
                    filtered_df["sseqid"].str.contains("phrog", case=False)
                ].shape[0],
                "repA": filtered_df[
                    filtered_df["sseqid"].str.contains("UniRef90", case=False)
                ].shape[0],
            }

            # if there are hits to more than 1 of dnaA, terL, repA, implement logic
            # to prefer dnaA, repA then terL (in that order)
            if (counts["dnaA"] > 0) + (counts["terL"] > 0) + (counts["repA"] > 0) >= 2:
                # prefer dnaA if it is greater than zero
                if counts["dnaA"] > 0:
                    # keep only the hits where dnaA is found
                    filtered_df = filtered_df[
                        filtered_df["sseqid"].str.contains("DNAA")
                    ]

                else:  # where there is repA and terL, keep repA
                    filtered_df = filtered_df[
                        filtered_df["sseqid"].str.contains("UniRef90")
                    ]

            # reset the index if you want to re-index the filtered DataFrame
            filtered_df.reset_index(drop=True, inplace=True)

            # top hit has a start of 1 ########
            # take the top hit - blast sorts by bitscore
            # if the start is 1 of the top hit
            if filtered_df["qstart"][0] == 1:
                # update the record description to contain 'rotated=True' akin to how unicycler does it
                record.description = record.description + " rotated=True"

                # writes to file
                with open(reoriented_output_file, "a") as out_fa:
                    SeqIO.write(record, out_fa, "fasta")

                # no hit save for the output DF
                message = "Contig_already_reoriented"

                genes.append(message)
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
                # get gene
                # set as dnaA by default
                gene = "dnaA"

                # for plasmids
                if "UniRef90" in filtered_df["sseqid"][0]:
                    gene = "repA"
                # for phages
                if "phrog" in filtered_df["sseqid"][0]:
                    gene = "terL"
                # custom
                if custom_db is not None:
                    gene = "custom"

                # update the record description to contain 'rotated=True' akin to how unicycler does it
                record.description = (
                    f"{record.description} rotated=True rotated_gene={gene}"
                )

                if filtered_df["qseq"][0][0] in ["M", "V", "L"] and (
                    filtered_df["sstart"][0] == 1
                ):
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

                else:  # top hit doesn't have a valid start codon - get most overlapping CDS
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

                # save all the stats
                genes.append(gene)
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
        "Gene_Reoriented": genes,
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
        os.path.join(output, prefix + "_all_reorientation_summary.tsv"),
        sep="\t",
        index=False,
    )


def run_autocomplete_record(
    record: SeqRecord, autocomplete: str, reoriented_output_file: Path, seed_value: int
) -> Tuple[int, int]:
    """
    Perform sequence reorientation and annotation based on specified autocomplete method.

    Args:
        record (SeqRecord): A Bio.SeqRecord object representing a nucleotide sequence.
        autocomplete (str): The autocomplete method to determine the gene to annotate ('mystery', 'nearest', 'largest').
        reoriented_output_file (Path): The output file path where the reoriented sequence will be saved.
        seed_value (int): Seed value for random number generation.

    Returns:
        Tuple[int, int]: A tuple containing the start coordinate and strand information of the chosen gene.
    """
    logger.warning(f"There was no blastx hit for contig {record.id}.")
    logger.warning(f"Running {autocomplete} on contig {record.id}.")

    # get number of records of input
    orf_finder = pyrodigal.GeneFinder(meta=True)

    # set seed
    random.seed(int(seed_value))

    # there will only be 1 record
    genes = orf_finder.find_genes(str(record.seq))
    # get number of genes
    gene_count = len(genes)

    if autocomplete == "mystery":
        # ensure has > 3 genes
        if gene_count < 4:
            logger.error(
                f"{record.id} has less than 4 CDS. You shouldn't be using -a mystery "
            )

        # ensure not first or last gene
        reorient_gene_number = random.randint(2, gene_count - 1)

        logger.info(
            f"Gene number {reorient_gene_number} was selected for contig {record.id}."
        )
        start = genes[reorient_gene_number].begin
        strand = genes[reorient_gene_number].strand

        if strand == 1:
            strand_eng = "forward"
        else:
            strand_eng = "negative"

        logger.info(
            f"Your random CDS has a start coordinate of {start} for contig {record.id}."
        )
        logger.info(
            f"Your random CDS is on the {strand_eng} strand for contig {record.id}."
        )

    elif autocomplete == "nearest":
        # ensure has > 1 genes
        if gene_count < 2:
            logger.error(
                f"{record.id} has less than 2 CDS. You shouldn't be using -a nearest"
            )

        reorient_gene_number = 1

        start = genes[reorient_gene_number].begin
        strand = genes[reorient_gene_number].strand

        if strand == 1:
            strand_eng = "forward"
        else:
            strand_eng = "negative"

        logger.info(
            f"Your nearest CDS has a start coordinate of {start} for contig {record.id}."
        )
        logger.info(
            f"Your nearest CDS is on the {strand_eng} strand for contig {record.id}."
        )

    elif autocomplete == "largest":
        if gene_count < 4:
            logger.error(
                f"{record.id} has less than 4 CDS. You probably shouldn't be using -a largest"
            )

        # iterate over dict

        size_dict = {}
        gene_index = 0

        for gene in genes:
            size = abs(gene.end - gene.begin)
            size_dict[gene_index] = size
            gene_index += 1

        # Find the gene with the max overlap
        largest_gene_index = max(size_dict, key=lambda key: size_dict[key])

        start = genes[largest_gene_index].begin
        strand = genes[largest_gene_index].strand
        max_size = round((size_dict[largest_gene_index] / 3), 2)

        if strand == 1:
            strand_eng = "forward"
        else:
            strand_eng = "negative"

        logger.info(
            f"Your largest CDS has a start coordinate of {start} for contig {record.id} with max size {max_size}."
        )
        logger.info(
            f"Your largest CDS is on the {strand_eng} strand for contig {record.id}."
        )

    reorient_sequence_and_append(record, reoriented_output_file, start, strand)

    return start, strand
