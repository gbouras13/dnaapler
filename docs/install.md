# Installing `dnaapler`

`dnaapler` requires only BLAST v2.9 or higher as an external dependency. 

The easiest way to install dnaapler is via conda. This method is highly recommended, because this will install BLAST automatically.

If you need instructions on how to install conda, please see the end of this page.

Conda
-----

`dnaapler` is available on bioconda.

```
conda install -c bioconda dnaapler
```

Pip
----

You can also install `dnaapler` with pip.

```
pip install dnaapler
```

If you install `dnaapler` with pip, then you will then need to install BLAST v 2.9 or higher separately. It will need to be available in the `$PATH` or else `dnaapler` will not work. 

Beginner Conda Installation
----

If you are new to using the command-line, please install conda using the following instructions.

First install [Anaconda](https://www.anaconda.com/products/distribution). There are lots of options but the two best in our opinion are:

   * [miniforge](https://github.com/conda-forge/miniforge).
   * [miniconda](https://docs.conda.io/en/latest/miniconda.html).
  
Please follow the instructions at the links to install based on your computer architecture. 

* Note: We would recommend install miniforge as it will automatically install mamba, which is much faster than base conda. 

After your installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

We would recommend installing dnaapler into a fresh environment. Assuming you installed miniforge, to create a environment called dnaaplerENV with dnaapler installed:

* To create a conda environment called `dnaapler_env`

```
conda create -n dnaapler_env
```

* To activate the environment

```
conda activate dnaapler_env
```

* To install dnaapler

```
mamba install -c bioconda dnaapler
```

* Once that has finished downloading and installing, you can check installation worked using:

```
dnaapler -h
```

* You should see:

```
Usage: dnaapler [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help     Show this message and exit.
  -V, --version  Show the version and exit.

Commands:
  all         Reorients contigs to begin with any of dnaA, repA or terL
  bulk        Reorients multiple genomes to begin with the same gene
  chromosome  Reorients your genome to begin with the dnaA chromosomal...
  citation    Print the citation(s) for this tool
  custom      Reorients your genome with a custom database
  largest     Reorients your genome the begin with the largest CDS as...
  mystery     Reorients your genome with a random CDS
  nearest     Reorients your genome the begin with the first CDS as...
  phage       Reorients your genome to begin with the terL large...
  plasmid     Reorients your genome to begin with the repA replication...
```