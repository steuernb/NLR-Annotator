# NLR-Annotator


This branch has an updated version of NLR-Parser. Now this should be compatible with up-to-date meme-suite versions. However, this is not fully tested yet. It should work but please let me know if you find any issues.

There is a sequence file for testing NLR-Parser: [Rgenes.CDS.fasta](https://github.com/steuernb/NLR-Annotator/blob/nlr_parser3/Rgenes.CDS.fasta). (Please note that these sequences are CDS of NLR genes. This is for testing NLR-Parser only!) 

## Introduction
NLR-Annotator is a tool to annotate loci associated with NLRs in large sequences.
It is searching for amino acid motifs within all 6 frames of a nucleotide sequence.
An NLR locus is defined from first to last motif that can be associated with an NLR. It does NOT predict genes. A predicted NLR locus might be a pseudogene. if it is overlapping with a gene, the actual gene start or intron-exon boundaries are not given. It just points you to the loci that might be worth investigating, which we hope you will find useful.


## Workflow
The NLR-Annotator pipeline consists of three steps: 

1. **Chopping the input sequence into overlapping sub-sequences.** This allows the usage of the Motif Alignment Search Tool (MAST) and also allows parallelization. The sub-sequences are overlapping to ensure no locus is missed because it was split by the chopping.
2. **Running NLR-Parser.** The chopped sub sequences are searched for NLR-associated motifs. The motifs are defined by [Jupe *et al*. 2012](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-75) (Table 1). The concept of NLR-Parser has been [published](https://academic.oup.com/bioinformatics/article/31/10/1665/177009) but for NLR-Annotator, the updated version bundled with this software is required.
3. **Running NLR-Annotator.** This program integrates the motifs annotated in the sub-sequences and searches for the actual NLR-Loci. Several output formats are provided.


## Requirements

### MEME suite

The MEME suite is available at http://meme.nbcr.net/meme/
Don't worry about setting up the Apache webserver. You just need MAST, so the quick install is sufficient.

### JRE 1.6

Make sure you have the Java Runtime Environments 1.6 or higher. Download from http://java.com



## Installation
Make sure MEME suite and JRE are installed. Download three jar files (ChopSequence.jar and NLR-Annotator.jar) as well as the motif definition file (meme.xml) from the [release](https://github.com/steuernb/NLR-Annotator/releases). 

For testing, please use the version of NLR-Parser in this branch: [NLR-Parser3.jar](https://github.com/steuernb/NLR-Annotator/blob/nlr_parser3/NLR-Parser3.jar)




## Running NLR-Annotator pipeline

### Chopping sequences. 
Input sequences are required to be nucleotide sequences and in FASTA format. Files may be gzip compressed. 


Usage 

```
java -jar ChopSequence.jar -i <inputsequence.fasta> -o <outputsequence.fasta> -l <sub-sequence length> -p <length of overlap>
```


Parameters
  
 parameter  | argument
---|---
-i | input file with sequence in fasta format
-o | output file with  sequence in fasta format
-l | (integer) length of sub-sequence (default 20000)
-p | (integer) length of overlap (default 5000)



### NLR-Parser

Usage

```
java -jar NLR-Parser.jar -t <number of threads> -y <path/to/meme/bin/mast> -x <path/to/meme.xml> -i <sub-seqeunces.fasta> -c <output.nlr.xml>
```

 parameter  | argument
---|---
-x | The path to the meme.xml (The motif definitions)
-y | The path to the mast installation (including the mast command. e.g. /programs/meme/bin/mast )
-i | The file with the chopped sequences. (Output from ChopSequence.jar)
-c | The output file that will be created. This is in xml format and will be the input for NLR-Annotator.


### NLR-Annotator

Usage

```
java -jar NLR-Annotator -i <nlr.xml> -o <output.nlr.txt>
```

 parameter  | argument | description
---|--- | ---
-i | input.xml | Input file in xml format. This is what comes out from NLR-Parser -c of the chopped inputSequence
-o | output.txt | output file in tabular format
-g | output.gff | output file in gff format
-b  |output.bed  | output file in bed format
-m |output.motifs.bed | output file of the motifs in bed format
-a | output.nbarkMotifAlignment.fasta| output file of the nb-arc motifs as multiple alignment. This file can be used as input to generate a phylogenetic tree. 
-f | genome.fasta output.nlr.fasta flanking | Write fasta of nlr loci. This parameter requires 3 arguments. The first is the original (not chopped) input sequence. The second is the file that is being generated. The third is the length of flanking sequence around the loci.
-distanceWithinMotifCombination | integer | (default:500)
-distanceForElongating | integer | (default:2500)
-distanceBetweenMotifCombinations |  integer | (default:10000)
			



