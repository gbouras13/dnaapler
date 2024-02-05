#!/usr/bin/env python3
"""dnaapler"""

import os
import shutil
from pathlib import Path

import click
from loguru import logger

from dnaapler.utils.all import all_process_blast_output_and_reorient
from dnaapler.utils.bulk import bulk_process_blast_output_and_reorient, run_bulk_blast
from dnaapler.utils.cds_methods import (
    run_blast_based_method,
    run_largest,
    run_mystery,
    run_nearest,
)
from dnaapler.utils.constants import DNAAPLER_DB
from dnaapler.utils.external_tools import ExternalTool
from dnaapler.utils.util import (
    begin_dnaapler,
    check_duplicate_headers,
    end_dnaapler,
    get_version,
    print_citation,
    run_autocomplete,
)
from dnaapler.utils.validation import (
    check_evalue,
    instantiate_dirs,
    validate_choice_autocomplete,
    validate_choice_db,
    validate_choice_mode,
    validate_custom_db_fasta,
    validate_fasta,
    validate_fasta_all,
    validate_fasta_bulk,
    validate_ignore_file,
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
            help="Path to input file in FASTA format",
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
            help="Number of threads to use with BLAST",
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
            help="Choose an option to autocomplete reorientation if BLAST based approach fails.\nMust be one of: none, mystery, largest, or nearest [default: none]",
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


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="e value for blastx",
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
    """Reorients your genome to begin with the dnaA chromosomal replication initiation gene"""

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "dnaA"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)
    logger.info(
        f"You have chosen {autocomplete} method to reorient your sequence if the BLAST based method fails."
    )

    # validates fasta
    validate_fasta(input)

    # validate e value
    check_evalue(evalue)

    # run BLAST
    blast_success = run_blast_based_method(
        ctx, input, output, prefix, gene, evalue, threads
    )

    # run autocomplete if BLAST reorientation failed
    run_autocomplete(
        blast_success, autocomplete, ctx, input, seed_value, output, prefix
    )

    # end dnaapler
    end_dnaapler(start_time)


"""
Plasmid command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="e value for blastx",
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
    """Reorients your genome to begin with the repA replication initiation gene"""

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "repA"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)
    logger.info(
        f"You have chosen {autocomplete} method to reorient your sequence if the BLAST based method fails."
    )

    # validates fasta
    validate_fasta(input)

    # validate e value
    check_evalue(evalue)

    # run BLAST
    blast_success = run_blast_based_method(
        ctx, input, output, prefix, gene, evalue, threads
    )

    # run autocomplete if BLAST reorientation failed
    run_autocomplete(
        blast_success, autocomplete, ctx, input, seed_value, output, prefix
    )

    # end dnaapler
    end_dnaapler(start_time)


"""
Phage command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="e value for blastx",
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
    """Reorients your genome to begin with the terL large terminase subunit"""

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "terL"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)
    logger.info(
        f"You have chosen {autocomplete} method to reorient your sequence if the BLAST based method fails."
    )

    # validates fasta
    validate_fasta(input)

    # validate e value
    check_evalue(evalue)

    # runs and processes BLAST
    blast_success = run_blast_based_method(
        ctx, input, output, prefix, gene, evalue, threads
    )

    # run autocomplete if BLAST reorientation failed
    run_autocomplete(
        blast_success, autocomplete, ctx, input, seed_value, output, prefix
    )

    # end dnaapler
    end_dnaapler(start_time)


"""
custom command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="e value for blastx",
    show_default=True,
)
@click.option(
    "-c",
    "--custom_db",
    help="FASTA file with amino acids that will be used as a custom blast database to reorient your sequence however you want.",
    type=click.Path(),
    required=True,
)
@click.option(
    "--ignore",
    default="",
    help="Text file listing contigs (one per row) that are to be ignored",
    type=click.Path(),
    show_default=False,
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
    ignore,
    **kwargs,
):
    """Reorients your genome with a custom database"""

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "custom"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)
    logger.info(
        f"You have chosen {autocomplete} method to reorient your sequence if the BLAST based method fails."
    )

    # validates fasta
    validate_fasta(input)

    # validate e value
    check_evalue(evalue)

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

    ExternalTool.run_tool(makeblastdb, ctx)

    # runs and processes BLAST
    blast_success = run_blast_based_method(
        ctx, input, output, prefix, gene, evalue, threads
    )

    # run autocomplete if BLAST reorientation failed
    run_autocomplete(
        blast_success, autocomplete, ctx, input, seed_value, output, prefix
    )

    # end dnaapler
    end_dnaapler(start_time)


"""
mystery command
"""


@main_cli.command()
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
    """Reorients your genome with a random CDS"""

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "mystery"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta
    validate_fasta(input)

    # run the mystery workflow
    run_mystery(ctx, input, seed_value, output, prefix)

    # finish dnaapler
    end_dnaapler(start_time)


"""
nearest command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
def nearest(ctx, input, output, threads, prefix, force, **kwargs):
    """Reorients your genome the begin with the first CDS as called by pyrodigal"""

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "nearest"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta
    validate_fasta(input)

    # run the nearest workflow
    run_nearest(ctx, input, output, prefix)

    # finish dnaapler
    end_dnaapler(start_time)


"""
nearest command
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
def largest(ctx, input, output, threads, prefix, force, **kwargs):
    """Reorients your genome the begin with the largest CDS as called by pyrodigal"""

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

    # defines gene
    gene = "nearest"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta
    validate_fasta(input)

    # run the nearest workflow
    run_largest(ctx, input, output, prefix)

    # finish dnaapler
    end_dnaapler(start_time)


"""
bulk subcommand
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="e value for blastx",
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
    help="FASTA file with amino acids that will be used as a custom blast database to reorient your sequence however you want. Must be specified if -m custom is specified.",
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
    """Reorients multiple genomes to begin with the same gene"""

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
    elif mode == "custom":
        gene == "custom"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta
    validate_fasta_bulk(input)
    check_duplicate_headers(input)

    # validate e value
    check_evalue(evalue)

    # check custom db
    if mode == "custom":
        if custom_db == None:
            logger.error(
                "You have specified dnaapler bulk -m custom without specifying a custom database using -c. Please try again using -c to input a custom database."
            )
        else:  # makes the database
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

            ExternalTool.run_tool(makeblastdb, ctx)
    else:
        if custom_db != None:
            logger.info(
                f"You have specified a custom database using -c but you have specified -m {mode}  not -m custom. Ignoring the custom database and continuing."
            )

    # runs  BLAST
    run_bulk_blast(ctx, input, output, prefix, gene, evalue, threads, custom_db)

    # rerorients blast
    blast_file = os.path.join(output, f"{prefix}_blast_output.txt")
    bulk_process_blast_output_and_reorient(input, blast_file, output, prefix)

    # end dnaapler
    end_dnaapler(start_time)


"""
all subcommand
"""


@main_cli.command()
@click.help_option("--help", "-h")
@click.version_option(get_version(), "--version", "-V")
@click.pass_context
@common_options
@click.option(
    "-e",
    "--evalue",
    default="1e-10",
    help="e value for blastx",
    show_default=True,
)
@click.option(
    "--ignore",
    default="",
    help="Text file listing contigs (one per row) that are to be ignored",
    type=click.Path(),
    show_default=False,
)
@click.option(
    "--db",
    "--database",
    default="all",
    type=click.STRING,
    callback=validate_choice_db,
    help="Lets you choose a subset of databases rather than all 3. Must be one of: 'all', 'dnaa', 'repa', terl', 'dnaa,repa', 'dnaa,terl' or 'repa,terl' ",
    show_default=True,
)
@click.option(
    "-c",
    "--custom_db",
    default="",
    help="FASTA file with amino acids that will be used as a custom blast database to reorient your sequence however you want.",
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
    """Reorients contigs to begin with any of dnaA, repA or terL"""

    # validates the directory  (need to before I start dnaapler or else no log file is written)
    instantiate_dirs(output, force)

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

    # custom
    if custom_db != "":
        gene = "custom"

    # initial logging etc
    start_time = begin_dnaapler(input, output, threads, gene)

    # validates fasta
    validate_fasta_all(input)
    check_duplicate_headers(input)

    # validate e value
    check_evalue(evalue)

    # create flag for ignore
    if ignore == "":
        ignore_flag = False
    else:
        ignore_flag = True
    # checks if the ignore file exists and contains text
    if ignore_flag == True:
        logger.info(f"You have specified contigs to ignore in {ignore}.")
        exists_contains_txt = validate_ignore_file(ignore)

    if gene == "custom":
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

        ExternalTool.run_tool(makeblastdb, ctx)
    else:
        custom_db = None

    # runs bulk BLAST
    run_bulk_blast(
        ctx, input, output, prefix, gene, evalue, threads, custom_db=custom_db
    )

    # rerorients blast
    blast_file = os.path.join(output, f"{prefix}_blast_output.txt")

    ### ignore
    # list is empty
    ignore_list = []
    if ignore_flag == True:
        if exists_contains_txt is False:
            logger.warning(f"{ignore} contains no text. No contigs will be ignored")
        else:
            # gets all contigs in the ignore
            # will split by space so short_contig only (to match BLAST)
            with open(ignore) as f:
                ignore_dict = {x.rstrip().split()[0] for x in f}
            ignore_list = list(ignore_dict)

    all_process_blast_output_and_reorient(
        input,
        blast_file,
        output,
        prefix,
        ignore_list,
        autocomplete,
        seed_value,
        custom_db=custom_db,
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
