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

Step 8: Context Free RNA folding
--------------------------------

- We use context free RNA folding with contrafold
- Using the ContraFold version 2.02 (http://contra.stanford.edu/contrafold/contrafold_v2_02.tar.gz)
- The folded files were further proceessed through the forgi tool to get the Forging encoding the corresponding secondary structure as denoted as Forgi vectors
- The bespoke script `fold-contrafold-uniq-bits-vectors.py` fold and generate Forgi secondary structure representation for all the transcripts
- Forgi tool and the associated utility rnaConvert was used(https://github.com/ViennaRNA/forgi)
- Corresponding script for the conversion of folded files to the dotracket formatted is `fold2dotbracketFasta.py` 

```
# run contrafold on merged_R1_R2.bit followed by forgi (https://github.com/ViennaRNA/forgi/blob/master/examples/rnaConvert.py)
python3 fold-contrafold-uniq-bits-vectors.py  --bit_file merged_R1_R2.bit --reference_file  cool6.fasta --transcript COOLAIR3 --size_file  sizer.tab
```

Step 9: PCA Projection
----------------------

- We use principal component analysis on the Forgi vectors using Scikit-learn: Machine Learning in Python(Pedregosa et al., JMLR 12, pp. 2825-2830, 2011)
- The projection helps us to view the grouping and separation the latent classes or configuration of structures
- Our bespoke script `run-pca-on-forgi-vectors.py` gives an output as a png and pdf format along with csv file for further manual insepction. 
 

```
python3  run-pca-on-forgi-vectors.py --input_file forgi-vect-ser.txt --tag COOLAIR3_R1_R2-forgi-pca  --csv_file COOLAIR3_R1_R2-forgi-pca.csv
```

Step 10: K-Means clustering
---------------------------

- Given a desired number of clusters we perform K-means clustering on the PCA projection using Scikit-learn
- User can specify the number of clusters,  we have set the default value to three 
- Clusters 1, 2 3 are coloured in the plot as green, orange and pink respectively.
- Final outputs will be a figure and a csv for this example `COOLAIR3_R1_R2-clusters.png`  and `COOLAIR3_R1_R2-clusters.csv`, repectively.


```
python3  draw-kmeans-clusters.py --input_file COOLAIR3_R1_R2-forgi-pca.csv --tag  COOLAIR3_R1_R2 --num_clusters 3

#inputs:
# input_file = COOLAIR3_R1_R2-forgi-pca.csv
# tag = COOLAIR3_R1_R2
# num_clusters = 3

#final outputs:
# - COOLAIR3_R1_R2-clusters.csv
# - COOLAIR3_R1_R2-clusters.png
# - COOLAIR3_R1_R2-clusters.pdf
```




