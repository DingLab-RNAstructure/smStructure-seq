# smStructure-seq

- We start with the raw pacbio subreads as the bam and corresponding pbi index files to generate the Single-Molecule Consensus Reads(HiFi Reads)  

- Generate consensus reads using CCS version 4.2.0(https://github.com/PacificBiosciences/ccs) with 10 passes.

- External tools that we need for the raw reads analysis namely, `pbccs`, `blasr`, `lima` and `bam2fasta` can be optionally through the bioconda channel.

- The other dependencies, `bam2fastx` version(https://github.com/PacificBiosciences/bam2fastx version: 21-05-2019) and `removesmartbell` utility from BBTools bioinformatics tools (version BBMap_38.60)(http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0185056) were also used.

step 1: Raw subreads to consensus HiFi reads
--------------------------------------------

- We used the CCS version 4.2.0 to generate Highly Accurate Single-Molecule Consensus Reads (HiFi Reads) from (https://github.com/PacificBiosciences/ccs)
- Both the input and output reads should have the associated index files.

```
# CCS: Generate Highly Accurate Single-Molecule Consensus Reads (HiFi Reads) 
# ccs 4.2.0 (commit v4.2.0)
ccs  -j 64  --minPasses=10  raw_reads/m54270_200819_091243.subreads.bam  consensus_reads/consensus_reads/enriched_coolair.ccs2.bam 
```

step 2: De-multiplex the consensus HiFi reads
--------------------------------------------

- Demultiplex the reads if you have have multiplexed the reads before the sequencing
- We use `primers.fa` with appended the directional primers names as 5p and 3p for the forward and reverse direction respectively.
- Using the lima version (lima 1.11.0),  a PacBio Barcode Demultiplexer and Primer Remover tool(https://github.com/PacificBiosciences/barcoding) 

```
# De-multiplex 
lima  --ccs consensus_reads/enriched_coolair.ccs2.bam   primers.fa   output.bam  --same  --split-bam-named --bam-handles 11  --bam-handles-verbose
```

step 3: Convert De-multiplexed consensus HiFi reads to Fasta 
------------------------------------------------------------

- Convert the reads to FASTA format for the downstream proceeses
- The conversion allows us to manually inspect and to be an input to many other toools
- bam2fasta version 1.3.1 from bam2fastx package (https://github.com/PacificBiosciences/bam2fastx)

```
# bam (or .subreadset.xml)  to Fasta
bam2fasta -u  --output  R1_5p    output.R1_5p--R1_5p.subreadset.xml
bam2fasta -u  --output  R2_5p    output.R2_5p--R2_5p.subreadset.xml
```

Step 4: Check and remove SMRTbell Barcoded Adapters
---------------------------------------------------

- Check and remove SMRTbell Barcoded Adapters(if any)
- We check for any residual SMRT Bell adapters with removesmartbell utility from BBMap version (38.60) (https://jgi.doe.gov/data-and-tools/bbtools)
- Please note this step may be optional

```
removesmartbell.sh in=R1_5p.fasta   out=R1_5p.fa  split=f adapter=CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC,GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
removesmartbell.sh in=R2_5p.fasta   out=R2_5p.fa  split=f adapter=CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC,GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
```

Step 5: BLASR mapping 
---------------------

- We have used the PacBiolong read aligner BLASR tool(version 5.1) to map the clean reads to the transcriptome reference (https://github.com/pacificbiosciences/blasr)
- The specific parameters we found useful hitPolicy as leftmost
- We used m5 format as output(as this allows to manual inspect the mapping results)

```
blasr --hitPolicy leftmost  --nproc 8  R1_5p.fasta   cool6.fasta  --minMatch 10   -m 5  --out    R1_5p.m5 
blasr --hitPolicy leftmost  --nproc 8  R2_5p.fasta   cool6.fasta  --minMatch 10   -m 5  --out    R2_5p.m5 
```

Step 6: Generate bitvectors from m5
-----------------------------------

- Oberved mutations profile for a transcript were made as the binary vectors representatation using the script `m5_to_bitvectors.py`
- The bitvectors represents sites of mutation and marked as 1 or 0 otherwise for the wild type state


```
python3 m5_to_bitvectors.py  --input_file  R1_5p.m5  --transcript COOLAIR3  --reference_file  cool6.fasta  --output_file  R1_5p.bit
python3 m5_to_bitvectors.py  --input_file  R2_5p.m5  --transcript COOLAIR3  --reference_file  cool6.fasta  --output_file  R2_5p.bit
```


Step 7: Merge bit vectors 
-------------------------

- This step is optional if you have not pooled the samples
- It gives us an opportunity to compare the replicates or samples
- Here we merged the pooled samples using the script `merge_bitvectors.py` 


```
python3 merge_bitvectors.py   --bit_file R1_5p.bit   R2_5p.bit    --output_file   merged_R1_R2.bit

```


