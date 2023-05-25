[![CI](https://github.com/gbouras13/dnaapler/actions/workflows/dnaapler_test.yml/badge.svg)](https://github.com/gbouras13/dnaapler/actions/workflows/dnaapler_test.yml)

[![codecov](https://codecov.io/gh/gbouras13/dnaapler/branch/refactor/graph/badge.svg?token=4B1T2PGM9V)](https://codecov.io/gh/gbouras13/dnaapler)


# dnaapler

Refactoring is in progress :) - install from source only for now. Will put on pypi and conda soon.

If you have any ideas or suggestions, please make an issue.

Installation
----------

dnaapler requires only BLAST as an external dependency.

For now, please install form source

```
git clone "https://github.com/gbouras13/dnaapler"
cd dnaapler

pip install -e .
dnaapler --help
```

and make sure you have blast installed and available.

```
conda install blast
```


Usage
----------

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
  phage       Reorients your sequence to begin with the terL large...
  plasmid     Reorients your sequence to begin with the repA replication...
  ```


Description
----------

`dnaapler` is a simple python program that takes a single nucleotide input sequence (in FASTA format) as input, finds the desired start gene using blastx, checks that the start of a gene is found, and if so, then reorients the chromosome to begin with this genes on the forward strand. 

It was designed to replicate the reorientation functionality of [Unicycler](https://github.com/rrwick/Unicycler/blob/main/unicycler/gene_data/repA.fasta) with dnaA, but for FASTA input, and primarily for long-read first assembled chromosomes. And I have extended it to work with plasmids and phages, or for any input FASTA desired.

For bacterial chromosome, `dnaapler chromosome` should ensure the chromosome breakpoint never interrupts genes or mobile genetic elements like prophages. It is intended to be used with good-quality completed bacterial genomes, generated with methods such as [Trycycler](https://github.com/rrwick/Trycycler/wiki) or [Dragonflye](https://github.com/rpetit3/dragonflye) or my own pipleine [hybracter](https://github.com/gbouras13/hybracter).

Databases
=============

`dnaapler chromosome` uses 733 proteins downloaded from Uniprot with the query "Chromosomal replication initiator protein DnaA" on 24 May 2023 as its database for dnaA. 

`dnaapler plasmid` uses the repA database curated by Ryan Wick in [Unicycler](https://github.com/rrwick/Unicycler/blob/main/unicycler/gene_data/repA.fasta).

`dnaapler phage` uses a terL database I curated using [PHROGs](https://phrogs.lmge.uca.fr). I downloaded all the AA sequences of the 55 phrogs annotated as 'large terminase subunit', combined them depduplicated them using [seqkit](https://github.com/shenwei356/seqkit) `seqkit rmdup -s -o terL.faa phrog_terL.faa`.

`dnaapler custom` uses a custom amino acid FASTA format gene(s) that you specify using `-c`. 

The matching is strict - it requires a strong BLAST match, and the first amino acid of a top BLAST hit gene to be identified as Methionine, Valine or Leucine, the 3 most used start codons in bacteria/phages. 

For the most commonly studied microbes (ESKAPE pathogens, etc), the dnaA database should suffice.

If you try dnaapler on a more novel or under-studied microbe with a dnaA gene that has little sequence similarity to the database, you may need to provide your own dnaA gene(s) in amino acid FASTA format using `dnaapler custom`.

After Erin Young's [issue](https://github.com/gbouras13/dnaapler/issues/1), I've also added:

`dnaapler mystery` which predicts all ORFs in the input, then picks a random sequence to re-orient your sequence with.s


Motivation
------------

1. I couldn't get [Circlator](https://sanger-pathogens.github.io/circlator/) to work and it is no longer supported.
2. [berokka](https://github.com/tseemann/berokka) doesn't orient chromosomes to begin with dnaa.
3. After reading Ryan Wick's masterful bacterial genome assembly [tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/wiki), I realised that it is probably optimal to run 2 polishing steps, once before then once after rotating the chromosome, to ensure the breakpoint is polished. Further, for some "complete" assemblies that didn't circularise properly, I figured that as long as you have a complete assembly (even if not "circular" as marked as in Flye), polishing after a re-orientation would be likely to circularise the chromosome. A bit like Ryan's [rotate_circular_gfa.py](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/blob/main/scripts/rotate_circular_gfa.py) script, without the requirement of strict circularity.
4. While researching MGEs in some of my s aureus whole genome sequences, I repeatedly found instances where MGEs were interrupted by the chromosome breakpoint. So I thought I'd add a tool to automate it in my pipeline. 
5. It's probably good to have all your sequences start at the same location for synteny analyses.

Polishing Afterwards
-----------

I recommend that you undertake 2 rounds of polishing. The first prior to running dnaapler, and then again after. I'd highly recommend a conservative polisher like [Polypolish](https://github.com/rrwick/Polypolish) if you have short reads, otherwise 2 rounds of medaka.

Acknowledgements
=============

Thanks to Torsten Seemann, Ryan Wick and the Circlator team for their existing work in the space. Also to [Michael Hall](https://github.com/mbhall88), whose repository [tbpore](https://github.com/mbhall88/tbpore) I took and adapted a lot of scaffolding code from because he writes really nice code, [Rob Edwards](https://github.com/linsalrob), because everything always comes back to phages, and [Vijini Mallawaarachchi]https://github.com/Vini2] who taught me how to actually do something resembling legitimate software development.


