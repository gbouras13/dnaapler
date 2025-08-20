#!/usr/bin/env python3
"""dnaapler"""

import os
from pathlib import Path

import click
from loguru import logger

from dnaapler.utils.all import all_process_MMseqs2_output_and_reorient
from dnaapler.utils.bulk import (
    bulk_process_MMseqs2_output_and_reorient,
    run_bulk_MMseqs2,
)
from dnaapler.utils.cds_methods import (
    run_largest,
    run_MMseqs2_based_method,
    run_mystery,
    run_nearest,
)
from dnaapler.utils.external_tools import ExternalTool
from dnaapler.utils.gfa import finalise_gfa, prep_gfa
from dnaapler.utils.processing import rotate_input
from dnaapler.utils.util import (
    begin_dnaapler,
    check_duplicate_headers,
    end_dnaapler,
    get_version,
    print_citation,
    remove_file,
    run_autocomplete,
)
from dnaapler.utils.validation import (
    check_evalue,
    instantiate_dirs,
    process_ignore_input,
    validate_choice_autocomplete,
    validate_choice_db,
    validate_choice_mode,
    validate_custom_db_fasta,
    validate_input,
    validate_input_all,
    validate_input_bulk,
)

"""
some code adapted from tbpore https://github.com/mbhall88/tbpore
"""

log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-i",
            "--input",
            help="Path to input file in FASTA or GFA format",
            type=click.Path(),
            required=True,
        ),
        click.option(
            "-o",
            "--output",
            default="output.dnaapler",
            show_default=True,
            type=click.Path(),
            help="Output directory ",
        ),
        click.option(
            "-t",
            "--threads",
            help="Number of threads to use with MMseqs2",
            default=1,
            show_default=True,
        ),
        click.option(
            "-p",
            "--prefix",
            default="dnaapler",
            help="Prefix for output files",
            show_default=True,
        ),
        click.option(
            "-f",
            "--force",
            is_flag=True,
            help="Force overwrites the output directory",
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


def autocomplete_options(func):
    """
    autocomplete options
    """
    options = [
        click.option(
            "-a",
            "--autocomplete",
            type=click.STRING,
            callback=validate_choice_autocomplete,
            default="none",
            help="Choose an option to autocomplete reorientation if MMseqs2 based approach fails.\nMust be one of: none, mystery, largest, or nearest [default: none]",
        ),
        click.option(
            "--seed_value",
            help="Random seed to ensure reproducibility.",
            type=int,
            default=13,
            show_default=True,
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
def main_cli():
    1 + 1


"""
Chromosome command
"""


@main_cli.command(
    short_help="Reorients your genome to begin with the dnaA chromosomal replication initiation gene"
)
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="E-value for MMseqs2",
    show_default=True,
)
@autocomplete_options
def chromosome(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    autocomplete,
    seed_value,
    **kwargs,
):
    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--autocomplete": autocomplete,
        "--seed_value": seed_value,
    }

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "dnaA"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)
    logger.info(
        f"You have chosen {autocomplete} method to reorient your sequence if the MMseqs2 based method fails."
    )

    # validates fasta or gfa
    validate_input(input)
    input_is_gfa, input, original_gfa = prep_gfa(input, output)

    # validate E-value
    check_evalue(evalue)

    # run MMseqs2
    MMseqs2_success = run_MMseqs2_based_method(
        ctx, input, output, prefix, gene, evalue, threads
    )

    # run autocomplete if MMseqs2 reorientation failed
    run_autocomplete(
        MMseqs2_success, autocomplete, ctx, input, seed_value, output, prefix
    )

    if input_is_gfa:
        finalise_gfa(
            input, original_gfa, os.path.join(output, f"{prefix}_reoriented.fasta")
        )

    # end dnaapler
    end_dnaapler(start_time)


"""
archaea command
"""


@main_cli.command(
    short_help="Reorients your genome to begin with the archaeal COG1474 Orc1/cdc6 origin recognition complex gene"
)
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="E-value for MMseqs2",
    show_default=True,
)
@autocomplete_options
def archaea(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    autocomplete,
    seed_value,
    **kwargs,
):
    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--autocomplete": autocomplete,
        "--seed_value": seed_value,
    }

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "cog1474"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)
    logger.info(
        f"You have chosen {autocomplete} method to reorient your sequence if the MMseqs2 based method fails."
    )

    # validates fasta or gfa
    validate_input(input)
    input_is_gfa, input, original_gfa = prep_gfa(input, output)

    # validate E-value
    check_evalue(evalue)

    # run MMseqs2
    MMseqs2_success = run_MMseqs2_based_method(
        ctx, input, output, prefix, gene, evalue, threads
    )

    # run autocomplete if MMseqs2 reorientation failed
    run_autocomplete(
        MMseqs2_success, autocomplete, ctx, input, seed_value, output, prefix
    )

    if input_is_gfa:
        finalise_gfa(
            input, original_gfa, os.path.join(output, f"{prefix}_reoriented.fasta")
        )

    # end dnaapler
    end_dnaapler(start_time)


"""
Plasmid command
"""


@main_cli.command(
    short_help="Reorients your genome to begin with the repA replication initiation gene"
)
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="E-value for MMseqs2",
    show_default=True,
)
@autocomplete_options
def plasmid(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    autocomplete,
    seed_value,
    **kwargs,
):
    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--autocomplete": autocomplete,
        "--seed_value": seed_value,
    }

    # defines gene
    gene = "repA"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)
    logger.info(
        f"You have chosen {autocomplete} method to reorient your sequence if the MMseqs2 based method fails."
    )

    # validates fasta or gfa
    validate_input(input)
    input_is_gfa, input, original_gfa = prep_gfa(input, output)

    # validate E-value
    check_evalue(evalue)

    # run MMseqs2
    MMseqs2_success = run_MMseqs2_based_method(
        ctx, input, output, prefix, gene, evalue, threads
    )

    # run autocomplete if MMseqs2 reorientation failed
    run_autocomplete(
        MMseqs2_success, autocomplete, ctx, input, seed_value, output, prefix
    )

    if input_is_gfa:
        finalise_gfa(
            input, original_gfa, os.path.join(output, f"{prefix}_reoriented.fasta")
        )

    # end dnaapler
    end_dnaapler(start_time)


"""
Phage command
"""


@main_cli.command(
    short_help="Reorients your genome to begin with the terL large terminase subunit"
)
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="E-value for MMseqs2",
    show_default=True,
)
@autocomplete_options
def phage(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    autocomplete,
    seed_value,
    **kwargs,
):
    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--autocomplete": autocomplete,
        "--seed_value": seed_value,
    }

    # defines gene
    gene = "terL"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)
    logger.info(
        f"You have chosen {autocomplete} method to reorient your sequence if the MMseqs2 based method fails."
    )

    # validates fasta or gfa
    validate_input(input)
    input_is_gfa, input, original_gfa = prep_gfa(input, output)

    # validate E-value
    check_evalue(evalue)

    # runs and processes MMseqs2
    MMseqs2_success = run_MMseqs2_based_method(
        ctx, input, output, prefix, gene, evalue, threads
    )

    # run autocomplete if MMseqs2 reorientation failed
    run_autocomplete(
        MMseqs2_success, autocomplete, ctx, input, seed_value, output, prefix
    )

    if input_is_gfa:
        finalise_gfa(
            input, original_gfa, os.path.join(output, f"{prefix}_reoriented.fasta")
        )

    # end dnaapler
    end_dnaapler(start_time)


"""
custom command
"""


@main_cli.command(short_help="Reorients your genome with a custom database")
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="E-value for MMseqs2",
    show_default=True,
)
@click.option(
    "-c",
    "--custom_db",
    help="FASTA file with amino acids that will be used as a custom MMseqs2 database to reorient your sequence however you want.",
    type=click.Path(),
    required=True,
)
@autocomplete_options
def custom(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    custom_db,
    autocomplete,
    seed_value,
    **kwargs,
):
    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--autocomplete": autocomplete,
        "--seed_value": seed_value,
        "--custom_db": custom_db,
    }

    # defines gene
    gene = "custom"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)
    logger.info(
        f"You have chosen {autocomplete} method to reorient your sequence if the MMseqs2 based method fails."
    )

    # validates fasta or gfa
    validate_input(input)
    input_is_gfa, input, original_gfa = prep_gfa(input, output)

    # validate E-value
    check_evalue(evalue)

    # validates custom fasta input for database
    validate_custom_db_fasta(Path(custom_db))

    # make db
    db_dir = os.path.join(output, "custom_db")
    Path(db_dir).mkdir(parents=True, exist_ok=True)

    logdir = Path(f"{output}/logs")

    # custom db
    # make custom db
    custom_database = os.path.join(db_dir, "custom_db")

    makeMMseqs2db = ExternalTool(
        tool="mmseqs",
        input=f"createdb {custom_db}",
        output=f" {custom_database}",
        params="",
        logdir=logdir,
    )

    ExternalTool.run_tool(makeMMseqs2db, ctx)

    # runs and processes MMseqs2
    MMseqs2_success = run_MMseqs2_based_method(
        ctx, input, output, prefix, gene, evalue, threads
    )

    # run autocomplete if MMseqs2 reorientation failed
    run_autocomplete(
        MMseqs2_success, autocomplete, ctx, input, seed_value, output, prefix
    )

    if input_is_gfa:
        finalise_gfa(
            input, original_gfa, os.path.join(output, f"{prefix}_reoriented.fasta")
        )

    # end dnaapler
    end_dnaapler(start_time)


"""
mystery command
"""


@main_cli.command(short_help="Reorients your genome with a random CDS")
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "--seed_value",
    help="Random seed to ensure reproducibility.",
    type=int,
    default=13,
    show_default=True,
)
def mystery(ctx, input, output, threads, prefix, seed_value, force, **kwargs):
    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--seed_value": seed_value,
    }

    # defines gene
    gene = "mystery"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)

    # validates fasta or gfa
    validate_input(input)

    # run the mystery workflow
    run_mystery(ctx, input, seed_value, output, prefix)

    # finish dnaapler
    end_dnaapler(start_time)


"""
nearest command
"""


@main_cli.command(
    short_help="Reorients your genome the begin with the first CDS as called by pyrodigal"
)
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
def nearest(ctx, input, output, threads, prefix, force, **kwargs):
    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
    }

    # defines gene
    gene = "nearest"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)

    # validates fasta or gfa
    validate_input(input)
    input_is_gfa, input, original_gfa = prep_gfa(input, output)

    # run the nearest workflow
    run_nearest(ctx, input, output, prefix)

    if input_is_gfa:
        finalise_gfa(
            input, original_gfa, os.path.join(output, f"{prefix}_reoriented.fasta")
        )

    # finish dnaapler
    end_dnaapler(start_time)


"""
largest command
"""


@main_cli.command(
    short_help="Reorients your genome the begin with the largest CDS as called by pyrodigal"
)
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
def largest(ctx, input, output, threads, prefix, force, **kwargs):
    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
    }

    # defines gene
    gene = "largest"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)

    # validates fasta or gfa
    validate_input(input)
    input_is_gfa, input, original_gfa = prep_gfa(input, output)

    # run the nearest workflow
    run_largest(ctx, input, output, prefix)

    if input_is_gfa:
        finalise_gfa(
            input, original_gfa, os.path.join(output, f"{prefix}_reoriented.fasta")
        )

    # finish dnaapler
    end_dnaapler(start_time)


"""
bulk subcommand
"""


@main_cli.command(short_help="Reorients multiple genomes to begin with the same gene")
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="E-value for MMseqs2",
    show_default=True,
)
@click.option(
    "-m",
    "--mode",
    type=click.STRING,
    callback=validate_choice_mode,
    default="none",
    help="Choose an mode to reorient in bulk.\nMust be one of: chromosome, plasmid, phage or custom [default: chromosome]",
)
@click.option(
    "-c",
    "--custom_db",
    help="FASTA file with amino acids that will be used as a custom MMseqs2 database to reorient your sequence however you want. Must be specified if -m custom is specified.",
    type=click.Path(),
)
def bulk(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    mode,
    force,
    custom_db,
    **kwargs,
):
    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "dnaA"
    if mode == "chromosome":
        gene = "dnaA"
    elif mode == "plasmid":
        gene = "repA"
    elif mode == "phage":
        gene = "terL"
    elif mode == "archaea":
        gene = "cog1474"
    elif mode == "custom":
        gene = "custom"

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--mode": mode,
        "--custom_db": custom_db,
    }

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)

    # validates fasta or gfa
    validate_input_bulk(input)
    check_duplicate_headers(input)
    input_is_gfa, input, original_gfa = prep_gfa(input, output)

    # validate E-value
    check_evalue(evalue)

    # check custom db
    logdir = Path(f"{output}/logs")
    if mode == "custom":
        if custom_db is None:
            logger.error(
                "You have specified dnaapler bulk -m custom without specifying a custom database using -c. Please try again using -c to input a custom database."
            )
        else:  # makes the database
            # validates custom fasta input for database
            validate_custom_db_fasta(Path(custom_db))

            # make db
            db_dir = os.path.join(output, "custom_db")
            Path(db_dir).mkdir(parents=True, exist_ok=True)

            # custom db
            # make custom db
            custom_database = os.path.join(db_dir, "custom_db")

            makeMMseqs2db = ExternalTool(
                tool="mmseqs",
                input=f"createdb {custom_db}",
                output=f" {custom_database}",
                params="",
                logdir=logdir,
            )

            ExternalTool.run_tool(makeMMseqs2db, ctx)

    else:
        if custom_db is not None:
            logger.info(
                f"You have specified a custom database using -c but you have specified -m {mode}  not -m custom. Ignoring the custom database and continuing."
            )

    # runs  MMseqs2
    run_bulk_MMseqs2(ctx, input, output, prefix, gene, evalue, threads, custom_db)

    # rerorients MMseqs2
    MMseqs2_file = os.path.join(output, f"{prefix}_MMseqs2_output.txt")
    bulk_process_MMseqs2_output_and_reorient(input, MMseqs2_file, output, prefix)

    if input_is_gfa:
        finalise_gfa(
            input, original_gfa, os.path.join(output, f"{prefix}_reoriented.fasta")
        )

    # end dnaapler
    end_dnaapler(start_time)


"""
all subcommand
"""


@main_cli.command(
    short_help="Reorients contigs to begin with any of dnaA, repA, terL or archaeal COG1474 Orc1/cdc6"
)
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="E-value for MMseqs2",
    show_default=True,
)
@click.option(
    "--ignore",
    default="",
    help="Text file listing contigs (one per row) that are to be ignored OR comma separated list of contig names to ignore OR '-' to read from stdin",
    type=str,
    show_default=False,
)
@click.option(
    "--db",
    "--database",
    default="all",
    type=click.STRING,
    callback=validate_choice_db,
    help="Lets you choose a subset of databases rather than all 4. Must be one of: 'all', 'dnaa', 'repa', terl', 'cog1474', 'dnaa,repa', 'dnaa,terl', 'repa,terl',  'dnaA,cog1474', 'cog1474,terl', 'cog1474,repa', 'dnaa,cog1474,repa', 'dnaa,cog1474,terl' or 'cog1474,repa,terl'",
    show_default=True,
)
@click.option(
    "-c",
    "--custom_db",
    default="",
    help="FASTA file with amino acids that will be used as a custom MMseqs2 database to reorient your sequence however you want.",
    type=click.Path(),
)
@autocomplete_options
def all(
    ctx,
    input,
    output,
    threads,
    prefix,
    evalue,
    force,
    autocomplete,
    seed_value,
    ignore,
    db,
    custom_db,
    **kwargs,
):
    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    params = {
        "--input": input,
        "--output": output,
        "--threads": threads,
        "--force": force,
        "--prefix": prefix,
        "--evalue": evalue,
        "--autocomplete": autocomplete,
        "--seed_value": seed_value,
        "--ignore": ignore,
        "--custom_db": custom_db,
        "--db": db,
    }

    # defines gene
    gene = "all"

    # other options

    if db == "dnaa":
        gene = "dnaA"
    elif db == "repa":
        gene = "repA"
    elif db == "terl":
        gene = "terL"
    elif db == "dnaa,repa":
        gene = "dnaA,repA"
    elif db == "dnaa,terl":
        gene = "dnaA,terL"
    elif db == "repa,terl":
        gene = "repA,terL"
    elif db == "cog1474":
        gene = "cog1474"
    elif db == "dnaA,cog1474":
        gene = "dnaA,cog1474"
    elif db == "cog1474,terl":
        gene = "cog1474,terL"
    elif db == "repa,cog1474":
        gene = "repA,cog1474"
    elif db == "dnaa,repa,cog1474":
        gene = "dnaa,repa,cog1474"
    elif db == "dnaa,terl,cog1474":
        gene = "dnaa,terl,cog1474"
    elif db == "repa,terl,cog1474":
        gene = "repa,terl,cog1474"

    # custom
    if custom_db != "":
        gene = "custom"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene, params)

    # validates fasta or gfa
    validate_input_all(input)
    check_duplicate_headers(input)
    input_is_gfa, input, original_gfa = prep_gfa(input, output)

    # validate E-value
    check_evalue(evalue)

    # process ignore input (file path or comma/space separated list)
    ignore_list = process_ignore_input(ignore)
    if ignore_list:
        logger.info(f"You have specified contigs to ignore: {', '.join(ignore_list)}")
    else:
        logger.info("No contigs specified to ignore.")

    if gene == "custom":
        # validates custom fasta input for database
        validate_custom_db_fasta(Path(custom_db))

        # make db
        db_dir = os.path.join(output, "custom_db")
        Path(db_dir).mkdir(parents=True, exist_ok=True)
        logdir = Path(f"{output}/logs")

        # custom db
        # make custom db
        custom_database = os.path.join(db_dir, "custom_db")

        makeMMseqs2db = ExternalTool(
            tool="mmseqs",
            input=f"createdb {custom_db}",
            output=f" {custom_database}",
            params="",
            logdir=logdir,
        )

        ExternalTool.run_tool(makeMMseqs2db, ctx)
    else:
        custom_db = None

    # rotate all replicons by half the length of the contig
    # the rotated input for MMSeqs2 will have the original contigs + the rotated ones with suffix "rotated_"
    rotated_input = os.path.join(output, "rotated_input.fasta")
    rotate_input(input, rotated_input)

    # runs bulk MMseqs2
    run_bulk_MMseqs2(
        ctx, rotated_input, output, prefix, gene, evalue, threads, custom_db=custom_db
    )

    # rerorients MMseqs2
    MMseqs2_file = os.path.join(output, f"{prefix}_MMseqs2_output.txt")

    all_process_MMseqs2_output_and_reorient(
        input,
        MMseqs2_file,
        output,
        prefix,
        ignore_list,
        autocomplete,
        seed_value,
        custom_db=custom_db,
        gene=gene,
    )

    # remove the rotated input
    remove_file(Path(rotated_input))

    if input_is_gfa:
        finalise_gfa(
            input, original_gfa, os.path.join(output, f"{prefix}_reoriented.fasta")
        )

    # end dnaapler
    end_dnaapler(start_time)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


# main_cli.add_command(run)
main_cli.add_command(citation)


def main():
    main_cli()


if __name__ == "__main__":
    main()
