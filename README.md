[![CI](https://github.com/gbouras13/dnaapler/actions/workflows/ci.yaml/badge.svg)](https://github.com/gbouras13/dnaapler/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/gbouras13/dnaapler/branch/refactor/graph/badge.svg?token=4B1T2PGM9V)](https://codecov.io/gh/gbouras13/dnaapler)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/dnaapler/badges/version.svg)](https://anaconda.org/bioconda/dnaapler)
[![Bioconda Downloads](https://img.shields.io/conda/dn/bioconda/dnaapler)](https://img.shields.io/conda/dn/bioconda/dnaapler)
[![PyPI version](https://badge.fury.io/py/dnaapler.svg)](https://badge.fury.io/py/dnaapler)
[![Downloads](https://static.pepy.tech/badge/dnaapler)](https://pepy.tech/project/dnaapler)


# dnaapler

## Quick Start

```
# creates conda environment with dnaapler
conda create -n dnaapler_env dnaapler

# activates conda environment
conda activate dnaapler_env

# runs dnaapler chromosome
dnaapler chromosome -i input.fasta -o output_directory_path -p my_bacteria_name -t 8
```


## Description


`dnaapler` is a simple python program that takes a single nucleotide input sequence (in FASTA format), finds the desired start gene using `blastx` against an amino acid sequence database, checks that the start codon of this gene is found, and if so, then reorients the chromosome to begin with this gene on the forward strand. 

It was originally designed to replicate the reorientation functionality of [Unicycler](https://github.com/rrwick/Unicycler/blob/main/unicycler/gene_data/repA.fasta) with dnaA, but for for long-read first assembled chromosomes. I have extended it to work with plasmids (`dnaapler plasmid`) and phages (`dnaapler phage`), or for any input FASTA desired with `dnaapler custom`, `dnaapler mystery` or `dnaapler nearest`.

For bacterial chromosomes, `dnaapler chromosome` should ensure the chromosome breakpoint never interrupts genes or mobile genetic elements like prophages. It is intended to be used with good-quality completed bacterial genomes, generated with methods such as [Trycycler](https://github.com/rrwick/Trycycler/wiki), [Dragonflye](https://github.com/rpetit3/dragonflye) or my own pipleine [hybracter](https://github.com/gbouras13/hybracter).

## Commands

* `dnaapler chromosome`: Reorients your sequence to begin with the dnaA chromosomal replication initiator gene
* `dnaapler plasmid`: Reorients your sequence to begin with the repA plasmid replication initiation gene
* `dnaapler phage`: Reorients your sequence to begin with the terL large terminase subunit gene
* `dnaapler custom`: Reorients your sequence to begin with a custom amino acid FASTA format gene that you specify
* `dnaapler mystery`: Reorients your sequence to begin with a random CDS
* `dnaapler nearest`: Reorients your sequence to begin with the first CDS (nearest to the start). Designed for fixing sequences where a CDS spans the breakpoint.

## Installation

`dnaapler` requires only BLAST as an external dependency. 

Installation from conda is recommended as this will install BLAST automatically.

### Conda

`dnaapler` is available on bioconda.

```
conda install -c bioconda dnaapler
```

### Pip

You can also install `dnaapler` with pip.

```
pip install dnaapler
```

You will need to install BLAST separately.

e.g.
`conda install -c bioconda blast`

## Usage

```
Usage: dnaapler [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  chromosome  Reorients your sequence to begin with the dnaA chromosomal...
  citation    Print the citation(s) for this tool
  custom      Reorients your sequence with a custom database
  mystery     Reorients your sequence with a random gene
  nearest     Reorients your sequence the begin with the first CDS as..
  phage       Reorients your sequence to begin with the terL large...
  plasmid     Reorients your sequence to begin with the repA replication...
  ```

  ```
  Usage: dnaapler chromosome [OPTIONS]

  Reorients your sequence to begin with the dnaA chromosomal replication
  initiation gene

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  -i, --input PATH       Path to input file in FASTA format  [required]
  -o, --output PATH      Output directory   [default: output.dnaapler]
  -t, --threads INTEGER  Number of threads to use with BLAST.  [default: 1]
  -p, --prefix TEXT      Prefix for output files.  [default :dnaapler]
  -f, --force            Force overwrites the output directory
  ```

The reoriented output FASTA will be `{prefix}_reoriented.fasta` in the specified output directory.

## Example Usage

```
dnaapler chromosome -i input.fasta -o output_directory_path -p my_bacteria_name -t 8
```

## Databases

`dnaapler chromosome` uses 584 proteins downloaded from Swissprot with the query "Chromosomal replication initiator protein DnaA" on 24 May 2023 as its database for dnaA. All hits from the query were also filtered to ensure "GN=dnaA" was included in the header of the FASTA entry.

`dnaapler plasmid` uses the repA database curated by Ryan Wick in [Unicycler](https://github.com/rrwick/Unicycler/blob/main/unicycler/gene_data/repA.fasta).

`dnaapler phage` uses a terL database curated using [PHROGs](https://phrogs.lmge.uca.fr). I downloaded all the AA sequences of the 55 phrogs annotated as 'large terminase subunit', combined them depduplicated them using [seqkit](https://github.com/shenwei356/seqkit) `seqkit rmdup -s -o terL.faa phrog_terL.faa`.

`dnaapler custom` uses a custom amino acid FASTA format file that you specify using `-c`. 

The matching is strict - it requires a strong BLAST match (default e-value 1E-10), and the first amino acid of a BLAST hit gene to be identified as Methionine, Valine or Leucine, the 3 most used start codons in bacteria/phages. 

For the most commonly studied microbes (ESKAPE pathogens, etc), the dnaA database should suffice.

If you try `dnaapler` on a more novel or under-studied microbe with a dnaA gene that has little sequence similarity to the database, you may need to provide your own dnaA gene(s) in amino acid FASTA format using `dnaapler custom`.

After this [issue](https://github.com/gbouras13/dnaapler/issues/1), `dnaapler mystery` was added. It predicts all ORFs in the input using [pyrodigal](https://github.com/althonos/pyrodigal), then picks a random gene to re-orient your sequence with

## Motivation

1. I couldn't get [Circlator](https://sanger-pathogens.github.io/circlator/) to work and it is no longer supported.
2. [berokka](https://github.com/tseemann/berokka) doesn't orient chromosomes to begin with dnaa.
3. After reading Ryan Wick's masterful bacterial genome assembly [tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/wiki), I realised that it is probably optimal to run 2 polishing steps, once before then once after rotating the chromosome, to ensure the breakpoint is polished. Further, for some "complete" long read bacterial assemblies that didn't circularise properly, I figured that as long as you have a complete assembly (even if not "circular" as marked as in Flye), polishing after a re-orientation would be likely to circularise the chromosome. A bit like Ryan's [rotate_circular_gfa.py](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/blob/main/scripts/rotate_circular_gfa.py) script, without the requirement of strict circularity.
4. While researching MGEs in _S. aureus_ whole genome sequences, I repeatedly found instances where MGEs were interrupted by the chromosome breakpoint. So I thought I'd add a tool to automate it in my pipeline. 
5. It's probably good to have all your sequences start at the same location for synteny analyses.

## Acknowledgements

Thanks to Torsten Seemann, Ryan Wick and the Circlator team for their existing work in the space. Also to [Michael Hall](https://github.com/mbhall88), whose repository [tbpore](https://github.com/mbhall88/tbpore) I took and adapted a lot of scaffolding code from because he writes really nice code, [Rob Edwards](https://github.com/linsalrob), because everything always comes back to phages, and especially [Vijini Mallawaarachchi](https://github.com/Vini2) who taught me how to actually do something resembling legitimate software development.


