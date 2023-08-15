# Output

### files

dnaapler creates a number of output files. For all subcommands that are not `dnaapler bulk`, inside the output directory you will find:

* A `{prefix}_reoriented.fasta` containing the reoriented genome.

* A `{prefix}_blast_output.txt` if a BLAST based method is used. This file will contain the raw blastx results in BLAST output format 6.

* A `.log` file containing the dnaapler log output (as printed to terminal).

* A  `logs` directory containing BLAST logs (if a BLAST based method was selected)

* If `dnaapler custom` is run, then a `custom_db` directory will also be present, containing the custom BLAST directory used by `dnaapler`.

### bulk

If you run `dnaapler bulk`, the output will be different. There will still be log files and a `{prefix}_blast_output.txt` file. The difference are:

* There will be 2 output `.fasta` files. One will be `{prefix}_reoriented.fasta` containing the reoriented contigs, while any contigs that failed to reorient will be in `{prefix}_failed_to_reorient.fasta`.

* There will be a `{prefix}_bulk_reorientation_summary.tsv` summary file containing the reorientation information for each contig. 

This summary file will contain columns for the contig name (Contig), the start coordinate of the reorientation (Start), the strandedness of this coordintate (Strand), the top BLAST hit name in the database (Top_Hit), the top hit length in amino acids (Top_Hit_Length), the length of the top BLAST hit covered in amino acids (Covered_Length), the coverage percentage (Coverage), the number of identical amino acids (Identical_AAs) and the identify percentage (Identity_Percentage).

For example for an input file with 3 contigs, where the first had no BLAST hit (and so failed to reorient), the second was already reoriented and the third was successfully reoriented, the summary file will look like:

| **Contig** | **Start**                 | **Strand**                | **Top_Hit**               | **Top_Hit_Length**        | **Covered_Length**        | **Coverage**              | **Identical_AAs**         | **Identity_Percentage**       |
|------------|---------------------------|---------------------------|---------------------------|---------------------------|---------------------------|---------------------------|---------------------------|---------------------------|
| contig_1    | No_BLAST_hits             | No_BLAST_hits             | No_BLAST_hits             | No_BLAST_hits             | No_BLAST_hits             | No_BLAST_hits             | No_BLAST_hits             | No_BLAST_hits             |
| contig_2       | Contig_already_reoriented | Contig_already_reoriented | Contig_already_reoriented | Contig_already_reoriented | Contig_already_reoriented | Contig_already_reoriented | Contig_already_reoriented | Contig_already_reoriented |
| contig_3     | 466148                    | reverse                   | sp\|Q6GD89\|DNAA_STAAS    | 453                       | 453                       | 100                       | 453                       | 100                       |

### all

If you run `dnaapler all`, the output will be different again. There will still be log files and a `{prefix}_blast_output.txt` file. 

* There will be 1 output `.fasta` file. One will be `{prefix}_all_reoriented.fasta` containing all the contigs.
  
* In this FASTA file, all  contigs that were reoriented will be indicated in the contig FASTA header with `rotated=True`.

* There will be a `{prefix}_all_reorientation_summary.tsv` summary file containing the reorientation information for each contig. 

This summary file will be the same as for `bulk`, but with an extra column `Gene_Reoriented` that denotes which gene was detected in each contig (dnaA, repA or terL).
