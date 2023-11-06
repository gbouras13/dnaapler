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
  - name: Susanna R. Grigson
    orcid: 0000-0003-4738-3451
    affiliation: 3
  - name: Bhavya Papudeshi
    orcid: 0000-0001-5359-3100
    affiliation: 3
  - name: Vijini Mallawaarachchi
    orcid: 0000-0002-2651-8719
    affiliation: 3
  - name: Michael J. Roach
    orcid: 0000-0003-1488-5148
    affiliation: 3
affiliations:
  - name: Adelaide Medical School, Faculty of Health and Medical Sciences, The University of Adelaide, Adelaide, South Australia 5005, Australia
    index: 1
  - name: The Department of Surgery – Otolaryngology Head and Neck Surgery, Central Adelaide Local Health Network, Adelaide, South Australia 5000, Australia 
    index: 2
  - name: Flinders Accelerator for Microbiome Exploration, College of Science and Engineering, Flinders University, Bedford Park, Adelaide, South Australia 5042, Australia
    index: 3
date: 08 August 2023  
bibliography: paper.bib
---

# Summary

Microorganisms found in natural environments are fundamental components of ecosystems and play vital roles in various ecological processes. Studying their genomes can provide valuable insights into the diversity, functionality, and evolution of microbial life, as well as their impacts on human health. Once the genetic material is extracted from environmental samples, it undergoes sequencing using technologies like whole genome sequencing (WGS). The raw sequence data is then analysed, and computational methods are applied to assemble the fragmented sequences and reconstruct the complete microbial genomes [@Wick:2021] [@Mallawaarachchi:2023] [@Bouras1:2023]. 

Many biological entities including Bacteria, Archaea, plasmids, bacteriophages and other viruses can have circular genomes. Once assembled, a circular genome sequence  is represented as a linear character string and labelled in some way to indicate that it should be circular. The point at which the linear sequence begins is random, due to the nature of the algorithms employed in assembling genomes from sequencing reads. Such arbitrary startpoints can affect downstream genome annotation and analysis; they may occur within coding sequences (CDS), can disrupt the prediction potential of mobile genetic elements like prophages, and make pangenome analyses based on gene order difficult. Therefore, microbial sequences are often required to be reoriented to begin by convention with certain genes: the dnaA chromosomal replication initiator gene for bacterial chromosomes, the repA plasmid replication initiation gene for plasmids and the terL large terminase subunit gene for bacteriophages as shown in \autoref{fig:workflow}. Here we present Dnaapler, a flexible microbial sequence reorientation tool that allows for rapid and consistent orientation of circular microbial genomes such as Bacteria, plasmids and bacteriophages. Dnaapler is hosted on GitHub at [github.com/gbouras13/dnaapler](https://github.com/gbouras13/dnaapler).

![Example microbial genome assembly workflow.\label{fig:workflow}](Dnaapler_figure.png){width=100%}

# Statement of need

Circlator [@Hunt:2015] is the most commonly used dedicated tool for reorienting bacterial genomes. However, Circlator was designed for bacterial chromosomes and plasmids only, is no longer supported by its developers, has several burdensome external dependencies, and requires the corrected reads in FASTA or FASTQ format along with the FASTA genome assembly as input. Alternatively, genome reorientation is often performed manually or with custom scripts on a genome-by-genome and project-by-project basis, making integration into assembly workflows difficult, and creating inconsistencies between different projects and researchers. We propose Dnaapler, a light-weight command-line tool written in Python that can be easily integrated into assembly workflows. Dnaapler takes only a FASTA formatted genome file as input. It uses BLAST [@Altschul:1990] [@Mount:2007] — its only external dependency — or Pyrodigal [@Larralde:2022] [@Hyatt:2010] depending on the chosen subcommand for reorientation. A list of the subcommands provided in Dnaapler are as follows:

| Subcommand       | Database used                                                                 | Gene used to reorient                       |
|------------|-------------------------------------------------------------------------------|---------------------------------------------|
| chromosome       | Custom database downloaded from Swissprot                                     | dnaA chromosomal replication initiator gene |
| plasmid          | repA database curated from Unicycler [@Wick:2017]                             | repA plasmid replication initiation gene    |
| phage            | Prokaryotic Virus Remote Homologous Groups database (PHROGs) [@Terzian:2021]  | terL large terminase subunit gene           |
| all              | Chromosome, plasmid and phage databases combined                              | dnaA, repA and terL                         |
| custom           | User specified                                                                | Custom gene                                 |
| mystery          | Pyrodigal predicted coding sequences                                          | Random CDS                                  |
| nearest          | Pyrodigal predicted coding sequences                                          | First CDS (nearest to the start)            |
| largest          | Pyrodigal predicted coding sequences                                          | Largest CDS                                 |
| bulk             | Either chromosome, plasmid, phage or custom. Requires multiple input contigs. | dnaA, repA, terL or a custom gene           |


Specifically, Dnaapler 'chromosome', 'phage' and 'plasmid' subcommands use blastx (protein databases are searched using a translated nucleotide query) to search for the dnaA, terL or repA gene respectively in the input genomes, using built-in amino acid databases for each gene. Dnaapler 'all' will run a blastx search against all three databases simultaneously, prioritising dnaA hits then repA and finally terL if multiple genes have hits. Taking the top blastx hit, Dnaapler will check that the first amino acid of the BLAST alignment begins with Methionine, Valine, or Leucine (the 3 most used gene start codons in bacteria and bacteriophages). If it does, then it will then reorient the genome to begin at that position in the forward direction. If it does not, then Pyrodigal will be used to predict all coding sequences. Dnaapler will calculate the CDS with the most overlap to the top blastx hit, and reorient the genome to begin with the start codon of that CDS in the forward direction. 

If the 'custom' subcommand is selected, the same process will be conducted but with a user specified amino acid FASTA formatted input database. If the 'mystery', 'nearest' or 'largest' subcommands are selected, Pyrodigal will be used to predict all coding sequences, and the genome will be reoriented to begin with either a random (mystery), the first (nearest) CDS, or the largest CDS respectively. Dnaapler returns an output directory containing a log file and the genome reoriented as a FASTA formatted file. Finally, the 'bulk' subcommand can be used to reorient multiple input contigs (in a mulitFASTA format file) using either the chromosome, plasmid, phage or custom reorientation strategies.

Examples of Dnaapler's functionality on the C333 _Staphylococcus aureus_ chromosome and the C333 Sa3int prophage (GenBank accession GCA_030288915.1, Sample Number SAMN32360890 from BioProject PRJNA914892 from [@Houtak:2023]) are shown below using the plotting functionalities of Bakta v1.8.2 [@Schwengers:2021] and Pharokka v1.5.1 [@Bouras2:2023].

![Example Dnaapler phage reorientation of the c333 Sa3int prophage as a circular genomic map from Pharokka beginning at the top of the circle. Each coloured arrow represents a CDS. The large terminase subunit gene is labelled as terL. Dnaapler identified the terL gene as beginning with coordinate 19146 on the forward strand. \label{fig:prophage}](C333_phage_combined.png){width=100%}

![Example Dnaapler chromosome reorientation of the C333 chromosome as a circular genomic map from Bakta beginning at the top of the circle. Each thin line represents a CDS, with the positive stranded CDSs denoted in the outer ring and the negatived stranded CDSs in the inner ring. The position of the chromosomal replication initiator gene is labelled as dnaA. The red and green ring denotes the GC content while the blue and yellow ring denotes the GC skew. Dnaapler identified the dnaA gene as beginning with coordinate 466140 on the reverse strand. \label{fig:prophage}](C333_chromosome_combined.png){width=100%}


# Availability

Dnaapler is distributed on PyPI. A [Conda](https://conda.io/) package is
also available in the Bioconda channel [@Bioconda:2018]. The source code is available on [GitHub](https://github.com/gbouras13/dnaapler),
and features Continuous Integration tests and test coverage, and Continuous Deployment using GitHub actions. Dnaapler has already been integrated into the United States of America StaPH-B (State Public Health Lab Bioinformatics) consortium [Docker image collection](https://github.com/StaPH-B/docker-builds).

# Acknowledgements
We would like to thank Michael B. Hall for providing some code snippets particularly the external tool class from his tool [tbpore](https://github.com/mbhall88/tbpore), Ryan Wick for curating a repA database from Unicycler and Sarah Vreugde and Robert A. Edwards for their supervision.

# References