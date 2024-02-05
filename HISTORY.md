# History

# 0.7.0 (2024-02-05)

* Adds `-c/--custom_db` with `dnaapler all` to allow specifying custom databases with `dnaapler all`.

# 0.6.0 (2024-01-31)

* Fixes bug where if the starting gene (dnaA/terL/repA) was on the reverse strand and the top BLAST hit did not find the start codon, it would reorient the replicon to begin at the end of the starting gene, not the start. Thanks @susiegriggo

# 0.5.2 (2024-01-24)

* Bumps version to include updated citation

# 0.5.1 (2024-01-09)

* With `dnaapler all`, adds the reoriented gene to the header (thanks @ammaraziz https://github.com/gbouras13/dnaapler/issues/67)
* Adds `--db` parameter to `dnaapler all`  allowing specifying a subset of genes to make up the database. In particular, if you have bacteria and plasmids, `--db dnaa,repa` should speed up Dnaapler's runtime quite a bit (thanks @oschwengers https://github.com/gbouras13/dnaapler/issues/63)

# 0.5.0 (2023-12-03)

* JOSS release with minor typos and bug fixes from v0.4

# 0.4.0 (2023-10-25)

* Implemented a modification to the logic for all cases where the top blastx hit alignment does not begin with a start codon. In this case, dnaapler will find the CDS according the pyrodigal that has the most overlap with the top hit alignment. Thanks @simone-pignotti for this suggestion [here](https://github.com/gbouras13/dnaapler/issues/44).
* Changes `dnaapler all` output FASTA to `_reoriented.fasta` instead of `_all_reoriented.fasta` for consistency with all other commands (except `dnaapler bulk`).
* Adds `-a` or `--autocomplete` option with `dnaapler all`.
* Adds `dnaapler largest` and `-a largest` as an option to orient your sequence beginning with the largest 

# 0.3.2 (2023-09-20)

* Changes `Orffinder` to `Genefinder`  to support `pyrodigal` v3.
* Updates dependency to `pyrodigal >=v3`.

# 0.3.1 (2023-09-01)

* Minor release to fix an error with dnaapler all #38 thanks @samnooij

# 0.3.0 (2023-08-18)

* `dnaapler all` subcommand added thanks @alexweisberg
* `dnaapler all` implements `--ignore` to ignore some contigs


# 0.2.0 (2023-08-08)

* `dnaapler nearest` subcommand added
* `dnaapler bulk` subcommand added
* dnaA database filtered to keep only bona-file dnaA genes (i.e. GN=dnaA)
* Adds `-e` parameter to vary BLAST evalue if desired
* Adds `-a` autocomplete parameter if user wants to reorient sequences with mystery or nearest methods in case the BLAST based method fails

# 0.1.0 (June 2023)

* Completely overhauled
* First stable released with pypi and conda 
* `dnaapler chromosome` added
* `dnaapler custom` added
* `dnaapler mystery` added 
* `dnaapler phage` added
* `dnaapler plasmid` added


# 0.0.1 (2022-10-12)

* First release (conda only `conda install -c gbouras dnaapler`)
