`dnaapler` is a simple python program that takes a single nucleotide input sequence (in FASTA or GFA format), finds the desired start gene using `blastx` against an amino acid sequence database, checks that the start codon of this gene is found, and if so, then reorients the chromosome to begin with this gene on the forward strand. 

It was originally designed to replicate the reorientation functionality of [Unicycler](https://github.com/rrwick/Unicycler/blob/main/unicycler/gene_data/repA.fasta) with dnaA, but for for long-read first assembled chromosomes. I have extended it to work with plasmids (`dnaapler plasmid`), phages (`dnaapler phage`) and archaea (`dnaapler archaea`), or for any input FASTA or GFA desired with `dnaapler custom`,`dnaapler largest`, `dnaapler mystery` or `dnaapler nearest`.

If your input FASTA or GFA is mixed and you have 1 or more contigs (e.g. has chromosome and plasmids), you should use `dnaapler all`, with the option to ignore some contigs with the `--ignore` parameter. This is probably the most useful command for most users.

When provided with a [GFA](https://gfa-spec.github.io/GFA-spec/GFA1.html) file, `dnaapler` will process only circular sequences – those with a single circularising link and no additional links – while leaving all other sequences unchanged. The output format will match the input: FASTA input produces FASTA output, and GFA input produces GFA output.

Additionally, you can also reorient multiple bacterial chromosomes/plasmids/phages at once using the `dnaapler bulk` subcommand - it will give you more information about what contigs couldn't be rotated which may be useful.

For bacterial chromosomes, `dnaapler chromosome` should ensure the chromosome breakpoint never interrupts genes or mobile genetic elements like prophages. It is intended to be used with good-quality completed bacterial genomes, generated with methods such as [Autocycler](https://github.com/rrwick/Autocycler/wiki), [Dragonflye](https://github.com/rpetit3/dragonflye) or [hybracter](https://github.com/gbouras13/hybracter).

A clear schematic of what dnaapler is trying to achieve is presented below:

![Image](Dnaapler_figure.png)