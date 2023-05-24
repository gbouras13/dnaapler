# dnaapler

```
git clone "https://github.com/gbouras13/dnaapler"
cd dnaapler

pip install -e .
dnaapler --help


```


dnaapler is a simple python program that takes a single complete whole genome chromosome as input, finds the dnaa gene using BLAST, checks that the start of the gene is found, and if so, then reorients the chromosome to begin with this dnaa on the forward strand. This will ensure the chromosome breakpoint never interrupts genes or mobile genetic elements like prophages.

dnaapler is intended to be used with good-quality completed genomes, generated with methods such as [Trycycler](https://github.com/rrwick/Trycycler/wiki) or with Flye and Illumina polishing.

dnaapler uses 650 dnaA proteins downloaded from Uniprot with the query "Chromosomal replication initiator protein DnaA" on 12-10-22 as its database. 

It is strict - it requires a strong BLAST match, and the first amino acid of the putative dnaA gene BLAST sequence to be identified as Methionine, Valine or Leucine, the 3 most used start codons in bacteria. 

For the most commonly studied microbes (ESKAPE pathogens, etc), the database should suffice.

If you try dnaapler on a more novel or under-studied microbe, you may need to provide your own dnaA gene in amino acid blast format. I will add this functionality later if required.

Thanks to Torsten Seemann, Ryan Wick and the Circlator team for their work in the space and general ideas.

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

Usage
----------

```
usage: dnaapler.py [-h] -c CHROMOSOME [-o OUTDIR] [-f] [-p PREFIX] [-t THREADS] [-V]

dnaapler: Orients complete bacterial whole genome chromosome assemblies to begin with dnaa gene.

optional arguments:
  -h, --help            show this help message and exit
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Bacterial chromosome assembly file in FASTA format.
  -o OUTDIR, --outdir OUTDIR
                        Directory to write the output to.
  -f, --force           Overwrites the output directory.
  -p PREFIX, --prefix PREFIX
                        Prefix for output files. This is not required
  -t THREADS, --threads THREADS
                        Threads for BLAST. Defaults to 4.
  -V, --version         show program's version number and exit
```
