---
title: 'Dnaapler: A tool to reorient circular microbial genomes'
tags:  
  - Python
  - microbiology
  - genomics
authors:
  - name: George Bouras
    orcid: 0000-0002-5885-4186
    affiliation: "1, 2"
  - name: Michael Roach
    orcid: 0000-0003-1488-5148
    affiliation: 3
affiliations:
  - name: Adelaide Medical School, Faculty of Health and Medical Sciences, The University of Adelaide, Adelaide, South Australia 5005, Australia
    index: 1
  - name: The Department of Surgery – Otolaryngology Head and Neck Surgery, Central Adelaide Local Health Network, Adelaide, South Australia 5000, Australia 
    index: 2
  - name: Flinders Accelerator for Microbiome Exploration, College of Science and Engineering, Flinders University, Bedford Park, Adelaide, South Australia 5042, Australia
    index: 3
date: 04 August 2023  
bibliography: paper.bib
---

# Summary

Due to algorithms employed in assembling genomes from sequencing reads, microbes with circular genomes such as bacteria, plasmids and bacteriophages often begin at a random position. Such arbitrary startpoints can affect downstream genome annotation and analysis; they may occur within coding sequences (CDS), can disrupt the prediction potential mobile genetic elements like prophages, and make pangenome analyses based on gene order difficult. Therefore, microbial sequences are often required to be reoriented to begin by convention with certain genes: the dnaA chromosomal replication initiator gene for bacterial chromosomes, the repA plasmid replication initiation gene for plasmids and the terL large terminase subunit gene for bacteriophages. Here we present Dnaapler, a flexible microbial sequence reorientation tool that allows for rapid and consistent orientation of circular microbial genomes such as bacteria, plasmids and bacteriophages. 

# Statement of need

Circlator [@Hunt:2015] is the most commonly used dedicated tool for reorienting bacterial genomes. However, Circlator was designed for bacterial chromosomes and plasmids only, is no longer supported by its developers, has a number of external dependencies and requires corrected reads in FASTA or FASTQ format along with the FASTA genome assembly as input. Otherwise, genome reorientation is usually done manually or with custom scripts on a genome-by-genome and project-by-project basis, making integration into assembly workflows difficult. 

We propose Dnaapler, a simple command-line tool written in Python 3 that can easily be integrated into assembly workflows. Dnaapler takes only a single FASTA formatted genome file as input. It uses the Basic Local Alignment Search Tool (BLAST) [@altschul:1990; @mount:2007] or Pyrodigal [@Larralde:2022] [@Hyatt:2010] depending on the chosen subcommand for reorientation. Specifically, Dnaapler 'chromosome', 'phage' and 'plasmid' subcommands use blastx (protein databases are searched using a translated nucleotide query) to search for the dnaA, terL or repA gene respectively in the input genomes, using built-in amino acid databases for each gene. Dnaapler will then check that the first amino acid of the identified start site begins with either Methionine, Valine or Leucine, the 3 most used gene start codons in bacteria and bacteriophages, and will then reorient the genome to begin with this gene forward stranded. If the 'custom' subcommand is selected, the same process will be conducted but with a user specified amino acid FASTA formatted input database. If the 'mystery' or 'nearest' subcommands are selected, Pyrodigal will be used to predict all coding sequences, and the genome will be reoriented to begin with either a random (mystery) or the first (nearest) CDS, respectively. Dnaapler returns an output directory containing a log file and a the genome reoriented as a FASTA formatted file.

Dnaapler has already been integrated into the United States of America StaPH-B (State Public Health Lab Bioinformatics) consortium [Docker image collection](https://github.com/StaPH-B/docker-builds).

# Availability

Dnaapler is distributed on PyPI. A [Conda](https://conda.io/) package is
also available in the Bioconda channel [@Bioconda:2018]. The source code is available on [GitHub](https://github.com/gbouras13/dnaapler),
and features tests via a Continuous Integration workflow implemented with Github actions.

# Acknowledgements
We would like to thank Michael B. Hall for providing some code snippets particularly the external tool class from his great tool tbpore, Ryan Wick for curating a repA database from Unicycler and Sarah Vreugde for her supervision.

# References