# Installing `dnaapler`

`dnaapler` requires only BLAST v2.9 or higher as an external dependency. 

The easiest way to install dnaapler is via conda. For inexperienced command line users, this method is highly recommended, because this will install BLAST automatically.

If you need instructions on how to install conda, please see the end of this page.

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

You will then need to install BLAST v 2.9 or higher separately.

e.g.

```
conda install -c bioconda blast>2.8
```

### Beginner Conda Installation

If you are new to using the command-line, please install conda using the following instructions.

1. Install [Anaconda](https://www.anaconda.com/products/distribution). I would recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html).
2. Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate one on the [miniconda](https://docs.conda.io/en/latest/miniconda.html) website).

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

* For Mac (Intel, will also work with M1):

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

3. Install miniconda and follow the prompts.

`sh Miniconda3-latest-Linux-x86_64.sh`

4. After installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

5. After this, conda should be installed (you may need to restart your terminal). It is recommended that mamba is also installed, as it will solve the enviroment quicker than conda:

`conda install mamba`

6. Finally, I would recommend installing dnaapler into a fresh environment. For example to create an environment called dnaaplerENV with dnaapler installed:

```
mamba create -n dnaaplerENV dnaapler
conda activate dnaaplerENV
dnaapler --help
```
