# NLR-Annotator version 2


**Under construction...please don't use yet!** 






----------------------

----------------------



- Code rewritten
- Integrated motif search tool. No dependency of MEME.

## Introduction
NLR-Annotator is a tool to annotate loci associated with NLRs in large sequences.
It is searching for amino acid motifs within all 6 frames of a nucleotide sequence.
An NLR locus is defined from first to last motif that can be associated with an NLR. 

It does NOT predict genes. A predicted NLR locus might be a pseudogene. 
If it is overlapping with a gene, the actual gene start or intron-exon boundaries are not given. 
It just points you to the loci that might be worth investigating, which we hope you will find useful.
Details are published in [Steuernagel et al.: The NLR-Annotator tool enables annotation of the intracellular immune receptor repertoire, Plant Physiology, 2020](https://www.ncbi.nlm.nih.gov/pubmed/32184345)



### JRE 1.6

Make sure you have the Java Runtime Environments 1.8 or higher. Download from http://java.com



## Installation
JRE is installed. Download jar file from [here](https://github.com/steuernb/NLR-Annotator/blob/NLR-Annotator-2/NLR-Annotator-v2.0.jar) 

You will need two config files, [mot.txt](https://github.com/steuernb/NLR-Annotator/blob/NLR-Annotator-2/src/mot.txt) and [store.txt](https://github.com/steuernb/NLR-Annotator/blob/NLR-Annotator-2/src/store.txt).


## Running NLR-Annotator

Run NLR-Annotator with `java -jar NLR-Annotator-v2.0.jar -i input.fasta -x mot.txt -y store.txt -o output.txt`

Replace "input.fasta" with the file (in fasta format) that contains the nucleotide sequences you want to annotate. 
If mot.txt and store.txt are not in your current directory, add the path.

If you use the `-t` parameter to run NLR-Annotator with multiple threads, please make sure the machine you are running it on has enough cores available.

In case you get an out-of-memory exception, add more to the java virtual machine, e.g. `java -Xmx8000M -jar NLR-Annotator.jar ...` will allow for 8000 MB.


### All parameters

 parameter  | argument | description
---|--- | ---
-i | input.fasta | Input file in fasta format. 
-o | output.txt | output file in tabular format
-g | output.gff | output file in gff format
-b  |output.bed  | output file in bed format
-m |output.motifs.bed | output file of the motifs in bed format
-a | output.nbarkMotifAlignment.fasta| output file of the nb-arc motifs as multiple alignment. This file can be used as input to generate a phylogenetic tree. 
-f | genome.fasta output.nlr.fasta flanking | Write fasta of nlr loci. This parameter requires 3 arguments. The first is the original (not chopped) input sequence. The second is the file that is being generated. The third is the length of flanking sequence around the loci.
-t | 1 | number of threads
-n | 1000 | number of sequences per thread (default 1000)
-x | mot.txt | motif file, contains internal config setting for motifs
-y | store.txt | store file, contains internal config settings for motifs

## Acknowledgements



