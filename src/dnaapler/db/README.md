# DB construction

* v1 (MMSeqs2) differs from v0 (BLAST) slightly as MMSeqs2 headers can't pass the pipe "|" https://github.com/soedinglab/MMseqs2/issues/804
* Therefore, changing this to "~" - this is not present in any of the FASTA files

```
for file in *.faa; do
    sed -i '' 's/|/~/g' "$file"
done
``` 

* To create the MMSeqs2 compatible DBs 

```
for file in *.faa; do
    # Extract the base filename without the extension
    base_name="${file%.faa}"
    # Run mmseqs createdb for each file
    mmseqs createdb "$file" "${base_name}_db"
done
```
