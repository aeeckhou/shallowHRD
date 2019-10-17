# shallowHRD

This method uses shallow Whole Genome Sequencing (sWGS ~ 0.5-2x) and the segmentation of a tumor genomic profile to infer the Homologous Recombination status of a tumor based on the number of Large-scale State Transtions (LSTs).

## Introduction

*shallowHRD* is a R script that can be launched from the command line. It relies on a ratio file characterizing the normalized read counts of a shallow Whole Genome Sequencing (0.5-2x) in sliding windows along the genome and its segmentation. It was developped on the [ControlFREEC](http://boevalab.inf.ethz.ch/FREEC/tutorial.html)'s output (Boeva,V. et al., 2011) but is adapted to other similar softwares (see sections "run shallowHRD" and "Nota Bene"). 

Softwares such as ControlFREEC count reads in sliding windows, correct the read count for GCcontent and low mappability region and then segment the genomic profile. *shallowHRD*, based on a inferred cut-off representing a one copy difference, will smooth the segmentation in a step wise manner, using first large segments, reintegrating small segments afterwards and then filtering small interstitial CNAs. The HR status is estimated based on the number of Large-scale State Transitions (LSTs) along the genome.

## Requirements

* R installed (tested with v.3.5.1)
* The following packages installed : 
  * ggpubr (tested with v.0.2.3)
  * gridExtra (tested with v.2.3)
  * DescTools (tested with v.0.99.28)

Tested on Linux and Mac.

## Prerequisities

First, FASTQ files should be aligned to the hg19 reference genome (using [BWA-MEM](https://github.com/lh3/bwa) for instance) and supplementary & duplicate reads removed from the BAM files, using [Samtools](http://www.htslib.org/doc/samtools.html) and [PicardTools' MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates), respectively.

Then, the BAM file should then be processed by a software such as ControlFREEC. The recommended options for controlFREEC are indicated in a config file example in the repository (*controlfreec_config_file_example*). The window size was fixed here to 20kb (coverage > 0.5x) and the parameters were set for a sensitive segmentation. The window size can however be increased up to ~60kb if necessary depending on the coverage, with a step size half its length.

Finally, the file *cytoBand_adapted_fr_hg19.csv* (available in the repository) has to be downloaded. 

The R packages needed can be installed with the script *install_packages.R* (available in repository) and the command line :

```
/path/to/Rscript /path/to/install.packages.R
```

## Run *shallowHRD*

To run *shallowHRD* only one ratio file is needed (formated in ControlFREEC's output).

The name of the file should be in this format : *SAMPLE_NAME.bam_ratio.txt*. <br/>

*shallowHRD* will rely on the *first four columns* of the input file (tabulated and with Chromosome in number !) : <br/>
Chromosome &nbsp; Start &nbsp; Ratio &nbsp; RatioMedian <br/>
1 &nbsp;&nbsp; 1 &nbsp;&nbsp; -1 &nbsp;&nbsp; -1 <br/>
1 &nbsp;&nbsp; 20001 &nbsp;&nbsp; -1 &nbsp;&nbsp; -1 <br/>
. &nbsp;&nbsp; . &nbsp;&nbsp; . &nbsp;&nbsp; . <br/>
. &nbsp;&nbsp; . &nbsp;&nbsp; . &nbsp;&nbsp; . <br/>

The command line to launch *shallowHRD* is (absolute or relative paths) :

```
/path/to/Rscript /path/to/shallowHRD.R /path/to/SAMPLE_NAME.bam_ratio.txt /path/to/output_directory /path/to/cytoBand_adapted_fr_hg19.csv
```

One file named *example.bam_ratio.txt* is downloadable in the repository to try *shallowHRD*.

## Outputs

All the figures and files created by the script will be available in the output directory.

The summary plot figure recapitulating all the information will look like this :

![alt text](https://github.com/aeeckhou/shallowHRD/blob/master/summary_plot_example.jpeg)

A : Genomic profile with LSTs in green (the entire processed segmentation is represented in red if there are no LST) <br/>
B : Density used to fix the difference between a copy level <br/>
C : Graphe representing the value of each final segment (small blue circles) - <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If the segmentation is good, the different copy number should appear clearly with disctinct steps <br/>
D : Table recapitulating different data, including the case quality and the final diagnostic for the HR status

## Nota Bene

1. The overall pipeline works also on WGS with a higher coverage

2. *shallowHRD* can be adapted to other segmentation outputs <br/> 
[QDNAseq](https://github.com/ccagc/QDNAseq) for instance : just comment the two script's lines for the log2 transformation of the data <br/>
(line 427 + 428 - sub-section "Fast gathering" of shallowHRD.R - 10/17/2019)

## Contact

Don't hesitate to contact us for any questions, problems or adaptation of the method !

alexandre.eeckhoutte@curie.fr <br/>
tatiana.popova@curie.fr <br/>
marc-henri.stern@curie.fr


## Publication

Submitted to __________.

"shallowHRD : a cheap and fast automatized method to evaluate Homologous Recombination Deficiency using shallow Whole Genome Sequencing" A. Eeckhoutte et al., 2019
