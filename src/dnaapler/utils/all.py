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


def all_process_blast_output_and_reorient(
    input, blast_file, output, prefix, ignore_list
):
    """Processes the blast output, reorients and saves all contigs into os.path.join(output, f"{prefix}_all_reoriented.fasta")

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

    reoriented_output_file = os.path.join(output, f"{prefix}_all_reoriented.fasta")

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

        # no hits at all then just write to the file
        elif length_of_df == 0:
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
                gene_found = False
                for i in range(0, len(filtered_df.qseq)):
                    if filtered_df["qseq"][i][0] in ["M", "V", "L"] and (
                        filtered_df["sstart"][i] == 1
                    ):
                        # update the record description to contain 'rotated=True' akin to how unicycler does it
                        record.description = record.description + " rotated=True"

                        # get gene
                        # set as dnaA by default
                        gene = "dnaA"

                        # for plasmids
                        if "UniRef90" in filtered_df["sseqid"][i]:
                            gene = "repA"
                        # for phages
                        if "phrog" in filtered_df["sseqid"][i]:
                            gene = "terL"

                        # if already reoriented
                        if filtered_df["qstart"][i] == 1:
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

                        else:
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
                                filtered_df, reoriented_output_file, record, i
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

                        gene_found = True

                        break

                if (
                    gene_found is False
                ):  # where there is no reorientiation no blast hits at a;;
                    with open(reoriented_output_file, "a") as out_fa:
                        SeqIO.write(record, out_fa, "fasta")

                    # no hit save for the output DF
                    message = "No_BLAST_hits"

                    genes.append(gene)
                    starts.append(message)
                    strands.append(message)
                    top_hits.append(message)
                    top_hit_lengths.append(message)
                    covered_lens.append(message)
                    coverages.append(message)
                    idents.append(message)
                    identitys.append(message)

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
