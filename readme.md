# NLR-Annotator version 2


**Under construction...please don't use yet!** 






----------------------

----------------------



- Code rewritten
- Integrated motif search tool. No dependency of MEME.

## Introduction
NLR-Annotator is a tool to annotate loci associated with NLRs in large sequences.
It is searching for amino acid motifs within all 6 frames of a nucleotide sequence.
An NLR locus is defined from first to last motif that can be associated with an NLR. It does NOT predict genes. A predicted NLR locus might be a pseudogene. if it is overlapping with a gene, the actual gene start or intron-exon boundaries are not given. It just points you to the loci that might be worth investigating, which we hope you will find useful.
Details are published in [Steuernagel et al.: The NLR-Annotator tool enables annotation of the intracellular immune receptor repertoire, Plant Physiology, 2020](https://www.ncbi.nlm.nih.gov/pubmed/32184345)



### JRE 1.6

Make sure you have the Java Runtime Environments 1.6 or higher. Download from http://java.com



## Installation
Make sure MEME suite and JRE are installed. Download three jar files (ChopSequence.jar, NLR-Parser.jar and NLR-Annotator.jar) as well as the motif definition file (meme.xml) from the [release](https://github.com/steuernb/NLR-Annotator/releases). Done.


## Running NLR-Annotator pipeline

