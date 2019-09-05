# shallowHRD

This method uses shallow Whole Genome Sequencing (sWGS ~ 1x) and its segmentation (via ControlFREEC) to infer the Homologous Recombination status of a tumor based on the number of Large-scale State Transtions (LSTs).


## Requirements

* R installed
* Linux or MACS OS
* The following packages installed : 
  * ggpubr
  * gridExtra
  * DescTools


## Prerequisities

shallowHRD is a R script that can be launched from the command line. It relies on the [ControlFREEC](http://boevalab.inf.ethz.ch/FREEC/tutorial.html)'s output (V. Boeava et al., 2011) on sWGS (0.5-2x). ControlFREEC counts reads in overlapping windows and corrects the read count for GCcontent, removing low mappability windows.

First, FASTQ files should be aligned to the hg19 reference genome (using [BWA-MEM](https://github.com/lh3/bwa) for instance) and supplementary & duplicate reads removed from the BAM files, using [Samtools](http://www.htslib.org/doc/samtools.html) and [PicardTools' MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates), respectively.

BAM files should then be processed by ControlFREEC. The recommended options are indicated in a config file example in the repository (controlfreec_config_file_example). Briefly, depending on the coverage of the BAM files (0.5-2x), the window size should varies between 20kb and 60kb, with a step size half its length (the lower the coverage, the larger the window).

Finally, the file *cytoBand_adapted_fr_hg19.csv* (available in the repository) has to be downloaded. 


## Run shallowHRD

To run shallowHRD, two files of the ControlFREEC's output are needed : the ratio and the info file. The names of the files should be in this format : *SAMPLE_NAME.bam_ratio.txt* and *SAMPLE_NAME.bam_info.txt*. Create a directory named after your sample (SAMPLE_NAME in the format exemple) and put the ratio and info file inside. Finally, put the *cytoBand_adapted_fr_hg19.csv* file in the folder containing the created directory.

The command line to launch the shallowHRD is :

```
Rscript /path/to/script/shallowHRD.R SAMPLE_NAME /path/to/folder/contaning/created/directory
```

## Outputs

![alt text](https://github.com/aeeckhou/shallowHRD/blob/master/summary_plot_example)

A : Genomic profile with LSTs in green (the entire processed segmentation is represented in red if there are no LST) <br/>
B : Density used to fix the difference between a copy level <br/>
C : Graphe representing the value of each final segment - if the segmentation is good, different should appear clearly <br/>
D : Table recapitulating different data and the final diagnostic for the HR status

## Nota Bene

1. The overall pipeline works also on BAM files with a coverage much more important.
2. The shallowHRD R script can be adapted to over kind of segmentation outputs ([QDNAseq](https://github.com/ccagc/QDNAseq) for instance)
3. The ploidy provided by ControlFREEC is just an indication and should only be taken as such

## Contact

alexandre.eeckhoutte@curie.fr <br/>
tatiana.popova@curie.fr <br/>
marc-henri.stern@curie.fr


## Publication

Submitted to BMC Bioinformatics.

"shallowHRD : a cheap and fast automatized method to evaluate Homologous Recombination Deficiency using shallow Whole Genome Sequencing" A. Eeckhoutte et al., 2019
