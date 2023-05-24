[![CI](https://github.com/gbouras13/dnaapler/actions/workflows/dnaapler_test.yml/badge.svg)](https://github.com/gbouras13/dnaapler/actions/workflows/dnaapler_test.yml)

[![codecov](https://codecov.io/gh/gbouras13/dnaapler/branch/refactor/graph/badge.svg?token=4B1T2PGM9V)](https://codecov.io/gh/gbouras13/dnaapler)


# dnaapler

Refactoring in progress :)

```
git clone "https://github.com/gbouras13/dnaapler"
cd dnaapler

pip install -e .
dnaapler --help
```


`dnaapler` is a simple python program that takes a single nucleotide sequence (in FASTA format) as input, finds the desired start gene using BLAST, checks that the start of the gene is found, and if so, then reorients the chromosome to begin with this genes on the forward strand. 

It was designed to replicate the reorientation functionality of [Unicycler](https://github.com/rrwick/Unicycler/blob/main/unicycler/gene_data/repA.fasta), but for FASTA input, and primarily for long-read first assembled chromosomes. But I have extended it to work with plasmids and phages, and any input FASTA desired.



For bacterial chromosome, `dnaapler chromosome` should ensure the chromosome breakpoint never interrupts genes or mobile genetic elements like prophages. It is intended to be used with good-quality completed bacterial genomes, generated with methods such as [Trycycler](https://github.com/rrwick/Trycycler/wiki) or [Dragonflye](https://github.com/rpetit3/dragonflye).

Databases
=============

`dnaapler chromosome` uses 733 proteins downloaded from Uniprot with the query "Chromosomal replication initiator protein DnaA" on 24 May 2023 as its database for dnaA. 

`dnaapler plasmid` uses the repA database curated by Ryan Wick in [Unicycler](https://github.com/rrwick/Unicycler/blob/main/unicycler/gene_data/repA.fasta).

`dnaapler phage` uses a terL database I curated using [PHROGs](https://phrogs.lmge.uca.fr). I downloaded all the AA sequences of the 55 phrogs annotated as 'large terminase subunit', combined them depduplicated them using [seqkit](https://github.com/shenwei356/seqkit).

```
seqkit rmdup -s -o terL.faa phrog_terL.faa
```

It is strict - it requires a strong BLAST match, and the first amino acid of the putative dnaA gene BLAST sequence to be identified as Methionine, Valine or Leucine, the 3 most used start codons in bacteria. 

For the most commonly studied microbes (ESKAPE pathogens, etc), the dnaA database should suffice.

If you try dnaapler on a more novel or under-studied microbe with a dnaA gene that has little sequence similarity to the database, you may need to provide your own dnaA gene in amino acid blast format using `dnaapler custom`.

Thanks to Torsten Seemann, Ryan Wick and the Circlator team for their work in the space and general ideas. Also to Michael Hall, whose repository [tbpore](https://github.com/mbhall88/tbpore) I took and adapted a lot of scaffolding code from, and [Rob Edwards](https://github.com/linsalrob), because everything always comes back to phages.

Motivation
------------

1. I couldn't get [Circlator](https://sanger-pathogens.github.io/circlator/) to work and it is no longer supported.
2. [berokka](https://github.com/tseemann/berokka) doesn't orient chromosomes to begin with dnaa.
3. After reading Ryan Wick's masterful bacterial genome assembly [tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/wiki), I realised that it is probably optimal to run 2 polishing steps, once before then once after rotating the chromosome, to ensure the breakpoint is polished. Further, for some "complete" assemblies that didn't circularise properly, I figured that as long as you have a complete assembly (even if not "circular" as marked as in Flye), polishing after a re-orientation would be likely to circularise the chromosome. A bit like Ryan's [rotate_circular_gfa.py](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/blob/main/scripts/rotate_circular_gfa.py) script, without the requirement of strict circularity.
4. While researching MGEs in some of my s aureus whole genome sequences, I repeatedly found instances where MGEs were interrupted by the chromosome breakpoint. So I thought I'd add a tool to automate it in my pipeline. 

Polishing Afterwards
-----------

I recommend that you undertake 2 rounds of polishing. The first prior to running dnaapler, and then again after. I'd highly recommend a conservative polisher like [Polypolish](https://github.com/rrwick/Polypolish) if you have short reads, otherwise 2 rounds of medaka.

Installation
----------

dnaapler requires only blast and biopython.

The easiest way to install is via conda either manually

```
git clone https://github.com/gbouras13/dnaapler.git
cd dnaapler
conda env create -f environment.yml
conda activate dnaapler_env
dnaapler.py -h
```

or via my conda channel.

```
conda install -c gbouras13 dnaapler
```
