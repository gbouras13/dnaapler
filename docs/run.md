
# Running `dnaapler`

For all subcommands, `dnaapler` requires an input FASTA file using the `-i` or `--input` parameters. 

It is also highly recommended to specify an output directory using the `-o` or `--output` parameters, otherwise `dnaapler` will write the output to a directory named `output.dnaapler` by default.

You can modify the prefix for the output files from `dnaapler` to whatever you please with the `-p` or `--prefix` parameters.

You can use BLAST with multiple threads using the `-t` or `--threads` parameters and modify the BLAST evalue with the `-e` or `--evalue` parameter.

`dnaapler` will not overwrite an output directory if it already exists by default. To force overwrite, please use `-f` or `--force`.

Finally, for the BLAST based subcommands (`chromosome`, `phage`, `plasmid`, `custom` or `all`), if no BLAST hit is found, by default `dnaapler` will error and exit. 

However, you can decide to autocomplete `dnaapler` using the `-a` or `--autocomplete` parameters along with `mystery` or `nearest`, which will then run those subcommands to reorient your sequence.

Also, a seed value using `--seed_value` can be specified with `dnaapler` to ensure that `dnaapler mystery` (or when austocomplete is used with `-a mystery`) to ensure `dnaapler` is reproducible in workflows.


### all

`dnaapler all` is designed to simultaneously orient multiple contigs that can be a mix of chromosomes, plasmids and phages. It will also work on just 1 contig.

If a contig has BLAST hits for both dnaA and terL or repA, dnaA will be chosen for reorientation.

If a contig has BLAST hits for both terL and repA (but not dnaA), repA will be chosen for reorientation.

You can also specify a text file with `--ignore` that lists all contigs (based on their header) to be ignored during reorientation.

e.g. the file (`ignored_contigs.txt`) needs to be formatted as follows:

```
contig_1
contig_2
```

Example usage to reorient a number of contigs in `input.fasta`, ignoring all contigs with headers denoted in  `ignored_contigs.txt`

```
dnaapler all -i input.fasta -o output_directory_path -t 8  --ignore ignored_contigs.txt
```

```
Usage: dnaapler all [OPTIONS]

  Reorients contigs to begin with any of dnaA, repA or terL

Options:
  -h, --help               Show this message and exit.
  -V, --version            Show the version and exit.
  -i, --input PATH         Path to input file in FASTA format  [required]
  -o, --output PATH        Output directory   [default: output.dnaapler]
  -t, --threads INTEGER    Number of threads to use with BLAST  [default: 1]
  -p, --prefix TEXT        Prefix for output files  [default: dnaapler]
  -f, --force              Force overwrites the output directory
  -e, --evalue TEXT        e value for blastx  [default: 1e-10]
  --ignore PATH            Text file listing contigs (one per row) that are to
                           be ignored
  -a, --autocomplete TEXT  Choose an option to autocomplete reorientation if
                           BLAST based approach fails. Must be one of: none,
                           mystery, largest, or nearest [default: none]
  --seed_value INTEGER     Rand
```


### chromosome

Example usage with `mystery` as the autocomplete command and a random seed of 245 for reproducibility and with 8 threads for BLAST:

```
dnaapler chromosome -i input.fasta -o output_directory_path -p my_bacteria_name -t 8 -a mystery --seed_value 245
```

```
Usage: dnaapler chromosome [OPTIONS]

  Reorients your genome to begin with the dnaA chromosomal replication
  initiation gene

Options:
  -h, --help               Show this message and exit.
  -V, --version            Show the version and exit.
  -i, --input PATH         Path to input file in FASTA format  [required]
  -o, --output PATH        Output directory   [default: output.dnaapler]
  -t, --threads INTEGER    Number of threads to use with BLAST  [default: 1]
  -p, --prefix TEXT        Prefix for output files  [default: dnaapler]
  -f, --force              Force overwrites the output directory
  -e, --evalue TEXT        e value for blastx  [default: 1e-10]
  -a, --autocomplete TEXT  Choose an option to autocomplete reorientation if
                           BLAST based approach fails. Must be one of: none,
                           mystery or nearest [default: none]
  --seed_value INTEGER     Random seed to ensure reproducibility.  [default:
                           13]
```

### phage

Example usage with no autocomplete command:

```
dnaapler phage -i input.fasta -o output_directory_path -p my_phage_name -t 8 
```

```
Usage: dnaapler phage [OPTIONS]

  Reorients your genome to begin with the terL large terminase subunit

Options:
  -h, --help               Show this message and exit.
  -V, --version            Show the version and exit.
  -i, --input PATH         Path to input file in FASTA format  [required]
  -o, --output PATH        Output directory   [default: output.dnaapler]
  -t, --threads INTEGER    Number of threads to use with BLAST  [default: 1]
  -p, --prefix TEXT        Prefix for output files  [default: dnaapler]
  -f, --force              Force overwrites the output directory
  -e, --evalue TEXT        e value for blastx  [default: 1e-10]
  -a, --autocomplete TEXT  Choose an option to autocomplete reorientation if
                           BLAST based approach fails. Must be one of: none,
                           mystery or nearest [default: none]
  --seed_value INTEGER     Random seed to ensure reproducibility.  [default:
                           13]
```

### plasmid

Example usage with no autocomplete command:

```
dnaapler plasmid -i input.fasta -o output_directory_path -p my_plasmid_name -t 8 
```

```
Usage: dnaapler plasmid [OPTIONS]

  Reorients your genome to begin with the repA replication initiation gene

Options:
  -h, --help               Show this message and exit.
  -V, --version            Show the version and exit.
  -i, --input PATH         Path to input file in FASTA format  [required]
  -o, --output PATH        Output directory   [default: output.dnaapler]
  -t, --threads INTEGER    Number of threads to use with BLAST  [default: 1]
  -p, --prefix TEXT        Prefix for output files  [default: dnaapler]
  -f, --force              Force overwrites the output directory
  -e, --evalue TEXT        e value for blastx  [default: 1e-10]
  -a, --autocomplete TEXT  Choose an option to autocomplete reorientation if
                           BLAST based approach fails. be one of: none,
                           mystery or nearest [default: none]
  --seed_value INTEGER     Random seed to ensure reproducibility.  [default:
                           13]
```

### custom

To run `dnaapler custom`, you need to prefix an Amino Acid FASTA file containing the desired custom database gene using `-c` or `--custom_db`. 

Example usage:

```
dnaapler custom -i input.fasta -o output_directory_path -p my_plasmid_name -t 8 -c custom_db.faa
```

```
Usage: dnaapler custom [OPTIONS]

  Reorients your genome with a custom database

Options:
  -h, --help               Show this message and exit.
  -V, --version            Show the version and exit.
  -i, --input PATH         Path to input file in FASTA format  [required]
  -o, --output PATH        Output directory   [default: output.dnaapler]
  -t, --threads INTEGER    Number of threads to use with BLAST  [default: 1]
  -p, --prefix TEXT        Prefix for output files  [default: dnaapler]
  -f, --force              Force overwrites the output directory
  -e, --evalue TEXT        e value for blastx  [default: 1e-10]
  -c, --custom_db PATH     FASTA file with amino acids that will be used as a
                           custom blast database to reorient your sequence
                           however you want.  [required]
  -a, --autocomplete TEXT  Choose an option to autocomplete reorientation if
                           BLAST based approach fails. Must be one of: none,
                           mystery or nearest [default: none]
  --seed_value INTEGER     Random seed to ensure reproducibility.  [default:
                           13]
```


### mystery

`dnaapler mystery` will reorient your genome to begin with a random coding sequence (CDS) (as predicted by [Pyrodigal](https://github.com/althonos/pyrodigal)).

Example usage:

```
dnaapler mystery -i input.fasta -o output_directory_path -t 8 
```

```
Usage: dnaapler mystery [OPTIONS]

  Reorients your genome with a random CDS

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  -i, --input PATH       Path to input file in FASTA format  [required]
  -o, --output PATH      Output directory   [default: output.dnaapler]
  -t, --threads INTEGER  Number of threads to use with BLAST  [default: 1]
  -p, --prefix TEXT      Prefix for output files  [default: dnaapler]
  -f, --force            Force overwrites the output directory
  --seed_value INTEGER   Random seed to ensure reproducibility.  [default: 13]
```

### nearest

`dnaapler nearest` will reorient your genome to begin the first coding sequence (CDS) as predicted by [Pyrodigal](https://github.com/althonos/pyrodigal).

Example usage:

```
dnaapler nearest -i input.fasta -o output_directory_path -t 8 
```

```
Usage: dnaapler nearest [OPTIONS]

  Reorients your genome the begin with the first CDS as called by pyrodigal

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  -i, --input PATH       Path to input file in FASTA format  [required]
  -o, --output PATH      Output directory   [default: output.dnaapler]
  -t, --threads INTEGER  Number of threads to use with BLAST  [default: 1]
  -p, --prefix TEXT      Prefix for output files  [default: dnaapler]
  -f, --force            Force overwrites the output directory
```

### bulk


`dnaapler bulk` is designed to simultaneously orient multiple genomes.
You must also specify `-m` or `--mode` with either `chromosome`, `phage`, `plasmid` or `custom` to tell `dnaapler` what mode to run. It will default to `-m chromosome`. Additionally, if you choose `-m custom`, then you must also specify a custom database amino acid file using `-c` or `--custom_db`.

Your input FASTA must also have at least 2 contigs.

Example usage to reorient a number of bacterial chromosomes in `input.fasta` to begin with the dnaA gene:

```
dnaapler bulk -i input.fasta -o output_directory_path -t 8  -m chromosome
```

```
Usage: dnaapler bulk [OPTIONS]

  Reorients multiple genomes to begin with the same gene

Options:
  -h, --help             Show this message and exit.
  -V, --version          Show the version and exit.
  -i, --input PATH       Path to input file in FASTA format  [required]
  -o, --output PATH      Output directory   [default: output.dnaapler]
  -t, --threads INTEGER  Number of threads to use with BLAST  [default: 1]
  -p, --prefix TEXT      Prefix for output files  [default: dnaapler]
  -f, --force            Force overwrites the output directory
  -e, --evalue TEXT      e value for blastx  [default: 1e-10]
  -m, --mode TEXT        Choose an mode to reorient in bulk. Must be one of:
                         chromosome, plasmid, phage or custom [default:
                         chromosome]
  -c, --custom_db PATH   FASTA file with amino acids that will be used as a
                         custom blast database to reorient your sequence
                         however you want. Must be specified if -m custom is
                         specified.
```

