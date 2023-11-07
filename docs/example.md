
# `dnaapler` examples

You can try running `dnaapler` yourself using test data found in the `tests/test_data/overall_inputs` as shown below. These examples assumed you have cloned the dnaapler repository from GitHub and have moved into the directory e.g. and have `dnaapler` [installed](install.md) :

```
git clone "https://github.com/gbouras13/dnaapler.git"
cd dnaapler
```

## Chromosome

This chromosome is from _Staphylococcus aureus_ isolate C333 taken from [Houtak et al](https://www.biorxiv.org/content/10.1101/2023.03.28.534496v1), GenBank accession GCA_030288915.1, Sample Number SAMN32360890 from BioProject PRJNA914892.

To run `dnaapler chromosome` to reorient the C333 chromosome to begin with the dnaA gene 

```
dnaapler chromosome -i tests/test_data/overall_inputs/chromosome.fasta -o C333_dnaapler -t 8 -p C333
```

The output should look like:

```
2023-11-07 12:08:51.243 | INFO     | dnaapler.utils.validation:instantiate_dirs:23 - Checking the output directory C333_dnaapler
2023-11-07 12:08:51.251 | INFO     | dnaapler.utils.util:begin_dnaapler:71 - You are using dnaapler version 0.4.0
2023-11-07 12:08:51.252 | INFO     | dnaapler.utils.util:begin_dnaapler:72 - Repository homepage is https://github.com/gbouras13/dnaapler
2023-11-07 12:08:51.252 | INFO     | dnaapler.utils.util:begin_dnaapler:73 - Written by George Bouras: george.bouras@adelaide.edu.au
2023-11-07 12:08:51.252 | INFO     | dnaapler.utils.util:begin_dnaapler:74 - Your input FASTA is tests/test_data/overall_inputs/chromosome.fasta
2023-11-07 12:08:51.252 | INFO     | dnaapler.utils.util:begin_dnaapler:75 - Your output directory  is C333_dnaapler
2023-11-07 12:08:51.252 | INFO     | dnaapler.utils.util:begin_dnaapler:76 - You have specified 8 threads to use with blastx
2023-11-07 12:08:51.252 | INFO     | dnaapler.utils.util:begin_dnaapler:77 - You have specified dnaA gene(s) to reorient your sequence
2023-11-07 12:08:51.252 | INFO     | dnaapler.utils.util:check_blast_version:115 - Checking BLAST installation.
2023-11-07 12:08:51.290 | INFO     | dnaapler.utils.util:check_blast_version:135 - BLAST version found is v2.14.1.
2023-11-07 12:08:51.290 | INFO     | dnaapler.utils.util:check_blast_version:145 - BLAST version is ok.
2023-11-07 12:08:51.290 | INFO     | dnaapler.utils.util:check_pyrodigal_version:90 - Checking pyrodigal installation.
2023-11-07 12:08:51.290 | INFO     | dnaapler.utils.util:check_pyrodigal_version:101 - Pyrodigal version is v3.1.1
2023-11-07 12:08:51.290 | INFO     | dnaapler.utils.util:check_pyrodigal_version:102 - Pyrodigal version is ok.
2023-11-07 12:08:51.290 | INFO     | dnaapler:chromosome:170 - You have chosen none method to reorient your sequence if the BLAST based method fails.
2023-11-07 12:08:51.290 | INFO     | dnaapler.utils.validation:validate_fasta:46 - Checking that the input file tests/test_data/overall_inputs/chromosome.fasta is in FASTA format and has only 1 entry.
2023-11-07 12:08:51.309 | INFO     | dnaapler.utils.validation:validate_fasta:53 - tests/test_data/overall_inputs/chromosome.fasta file checked.
2023-11-07 12:08:51.319 | INFO     | dnaapler.utils.validation:validate_fasta:62 - tests/test_data/overall_inputs/chromosome.fasta has only one entry.
2023-11-07 12:08:51.320 | INFO     | dnaapler.utils.validation:check_evalue:187 - You have specified an evalue of 1e-10.
2023-11-07 12:08:51.321 | INFO     | dnaapler.utils.external_tools:run:49 - Started running blastx -db ~/dnaapler/src/dnaapler/db/dnaA_db -evalue 1e-10 -num_threads 8 -outfmt ' 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq ' -out C333_dnaapler/C333_blast_output.txt -query tests/test_data/overall_inputs/chromosome.fasta ...
2023-11-07 12:09:01.769 | INFO     | dnaapler.utils.external_tools:run:51 - Done running blastx -db ~/dnaapler/src/dnaapler/db/dnaA_db -evalue 1e-10 -num_threads 8 -outfmt ' 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq ' -out C333_dnaapler/C333_blast_output.txt -query tests/test_data/overall_inputs/chromosome.fasta
2023-11-07 12:09:01.776 | INFO     | dnaapler.utils.processing:reorient_sequence:145 - dnaA gene identified. It starts at coordinate 466140 on the reverse strand in your input file.
2023-11-07 12:09:01.776 | INFO     | dnaapler.utils.processing:reorient_sequence:148 - The best hit with a valid start codon in the database is sp|Q6GKU4|DNAA_STAAR, which has length of 453 AAs.
2023-11-07 12:09:01.776 | INFO     | dnaapler.utils.processing:reorient_sequence:151 - 453 AAs were covered by the best hit, with an overall coverage of 100.0%.
2023-11-07 12:09:01.776 | INFO     | dnaapler.utils.processing:reorient_sequence:154 - 452 AAs were identical, with an overall identity of 99.78%.
2023-11-07 12:09:01.776 | INFO     | dnaapler.utils.processing:reorient_sequence:157 - Re-orienting.
2023-11-07 12:09:01.798 | INFO     | dnaapler.utils.util:end_dnaapler:158 - dnaapler has finished
2023-11-07 12:09:01.798 | INFO     | dnaapler.utils.util:end_dnaapler:159 - Elapsed time: 10.55 seconds
```

In the results in the output directory, you will see that the `C333_reorientation_summary.tsv` file shows that `dnaapler` has identified the C333 genome to begin with coordinate

A comparison of circular genomic maps of the C333 chromosome before and after `dnaapler` made with [Bakta v1.8.2](https://github.com/oschwengers/bakta) can be seen below. The genome begins at the top of the circle. Each thin line represents a CDS, with the positive stranded CDSs denoted in the outer ring and the negatived stranded CDSs in the inner ring. The position of the chromosomal replication initiator gene is labelled as dnaA. The red and green ring denotes the GC content while the blue and yellow ring denotes the GC skew. 

![Image](C333_chromosome_combined.png)

## Phage

This phage is the Sa3int prophage from the _Staphylococcus aureus_ isolate C333 described above.

To run `dnaapler phage` to reorient the C333 prophage to begin with the terL (terminase large subunit) gene: 

```
dnaapler phage -i tests/test_data/overall_inputs/C333_sa3int_phage.fasta -o C333_phage_dnaapler -t 8 -p C333_phage
```

The output should look like:

```
2023-11-07 12:24:14.227 | INFO     | dnaapler.utils.validation:instantiate_dirs:23 - Checking the output directory C333_phage_dnaapler
2023-11-07 12:24:14.234 | INFO     | dnaapler.utils.util:begin_dnaapler:71 - You are using dnaapler version 0.4.0
2023-11-07 12:24:14.234 | INFO     | dnaapler.utils.util:begin_dnaapler:72 - Repository homepage is https://github.com/gbouras13/dnaapler
2023-11-07 12:24:14.234 | INFO     | dnaapler.utils.util:begin_dnaapler:73 - Written by George Bouras: george.bouras@adelaide.edu.au
2023-11-07 12:24:14.234 | INFO     | dnaapler.utils.util:begin_dnaapler:74 - Your input FASTA is tests/test_data/overall_inputs/C333_sa3int_phage.fasta
2023-11-07 12:24:14.234 | INFO     | dnaapler.utils.util:begin_dnaapler:75 - Your output directory  is C333_phage_dnaapler
2023-11-07 12:24:14.234 | INFO     | dnaapler.utils.util:begin_dnaapler:76 - You have specified 8 threads to use with blastx
2023-11-07 12:24:14.234 | INFO     | dnaapler.utils.util:begin_dnaapler:77 - You have specified terL gene(s) to reorient your sequence
2023-11-07 12:24:14.234 | INFO     | dnaapler.utils.util:check_blast_version:115 - Checking BLAST installation.
2023-11-07 12:24:14.309 | INFO     | dnaapler.utils.util:check_blast_version:135 - BLAST version found is v2.14.1.
2023-11-07 12:24:14.309 | INFO     | dnaapler.utils.util:check_blast_version:145 - BLAST version is ok.
2023-11-07 12:24:14.309 | INFO     | dnaapler.utils.util:check_pyrodigal_version:90 - Checking pyrodigal installation.
2023-11-07 12:24:14.309 | INFO     | dnaapler.utils.util:check_pyrodigal_version:101 - Pyrodigal version is v3.1.1
2023-11-07 12:24:14.309 | INFO     | dnaapler.utils.util:check_pyrodigal_version:102 - Pyrodigal version is ok.
2023-11-07 12:24:14.309 | INFO     | dnaapler:phage:298 - You have chosen none method to reorient your sequence if the BLAST based method fails.
2023-11-07 12:24:14.309 | INFO     | dnaapler.utils.validation:validate_fasta:46 - Checking that the input file tests/test_data/overall_inputs/C333_sa3int_phage.fasta is in FASTA format and has only 1 entry.
2023-11-07 12:24:14.316 | INFO     | dnaapler.utils.validation:validate_fasta:53 - tests/test_data/overall_inputs/C333_sa3int_phage.fasta file checked.
2023-11-07 12:24:14.316 | INFO     | dnaapler.utils.validation:validate_fasta:62 - tests/test_data/overall_inputs/C333_sa3int_phage.fasta has only one entry.
2023-11-07 12:24:14.316 | INFO     | dnaapler.utils.validation:check_evalue:187 - You have specified an evalue of 1e-10.
2023-11-07 12:24:14.317 | INFO     | dnaapler.utils.external_tools:run:49 - Started running blastx -db ~/dnaapler/src/dnaapler/db/terL_db -evalue 1e-10 -num_threads 8 -outfmt ' 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq ' -out C333_phage_dnaapler/C333_phage_blast_output.txt -query tests/test_data/overall_inputs/C333_sa3int_phage.fasta ...
2023-11-07 12:24:15.180 | INFO     | dnaapler.utils.external_tools:run:51 - Done running blastx -db ~/dnaapler/src/dnaapler/db/terL_db -evalue 1e-10 -num_threads 8 -outfmt ' 6 qseqid qlen sseqid slen length qstart qend sstart send pident nident gaps mismatch evalue bitscore qseq sseq ' -out C333_phage_dnaapler/C333_phage_blast_output.txt -query tests/test_data/overall_inputs/C333_sa3int_phage.fasta
2023-11-07 12:24:15.188 | INFO     | dnaapler.utils.processing:reorient_sequence:145 - terL gene identified. It starts at coordinate 19146 on the forward strand in your input file.
2023-11-07 12:24:15.188 | INFO     | dnaapler.utils.processing:reorient_sequence:148 - The best hit with a valid start codon in the database is phrog_9_p344137, which has length of 553 AAs.
2023-11-07 12:24:15.188 | INFO     | dnaapler.utils.processing:reorient_sequence:151 - 553 AAs were covered by the best hit, with an overall coverage of 100.0%.
2023-11-07 12:24:15.188 | INFO     | dnaapler.utils.processing:reorient_sequence:154 - 552 AAs were identical, with an overall identity of 99.82%.
2023-11-07 12:24:15.188 | INFO     | dnaapler.utils.processing:reorient_sequence:157 - Re-orienting.
2023-11-07 12:24:15.189 | INFO     | dnaapler.utils.util:end_dnaapler:158 - dnaapler has finished
2023-11-07 12:24:15.189 | INFO     | dnaapler.utils.util:end_dnaapler:159 - Elapsed time: 0.96 seconds
```

In the results in the output directory, you will see that the `C333_phage_reorientation_summary.tsv` file shows that `dnaapler` has identified the C333 genome to begin with coordinate 19146 on the forward strand.

A comparison of circular genomic maps of the C333 Sa3int prophage before and after `dnaapler` made with [Pharokka v1.5.1](https://github.com/gbouras13/pharokka) can be seen below. The genome begins at the top of the circle. Each coloured arrow represents a CDS.

![Image](C333_phage_combined.png)
