# smStructure-seq

- We start with the raw pacbio subreads as the bam and corresponding pbi index files to generate the Single-Molecule Consensus Reads(HiFi Reads)  

- Generate consensus reads using CCS version 4.2.0(https://github.com/PacificBiosciences/ccs) with minimum 10 passes

- External tools that we need for the raw reads analysis namely, `pbccs`, `blasr`, `lima` and `bam2fasta` can be optionally installed inside the containers through the bioconda channel.

- The other dependencies, `bam2fastx` version(https://github.com/PacificBiosciences/bam2fastx version: 21-05-2019) and `removesmartbell` utility from BBTools bioinformatics tools (version BBMap_38.60)(http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0185056) are also installed inside the container

