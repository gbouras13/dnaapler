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

```
mamba create -n dnaaplerENV dnaapler
conda activate dnaaplerENV
dnaapler --help
```
