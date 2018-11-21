# Phoenix

This pipeline consists of subtemplates (jst files found throughout this 
directory) that will only be used if specific data types are present in the 
project. E.g. If RNA FASTQs are present, the RNA module will be used, but 
otherwise it will not be included. "Features" are enabled/disabled purely
by the content of the project. This pattern is discussed in more detail in the
[implementation](implementation) section.

# Features

## DNA Alignment

The DNA Alignment module looks for fastq data in the read_groups section of
the config. Read groups will be aligned, post-processed, and qc'ed.

## Somatic Variant Calling

This module will run several somatic variant callers on tumor/normal data pairs.

## Germline Variant Calling

Generates germline variant call files (VCF) with several callers. Additionally,
this will create a gVCF for each sample that can be used to jointly call large
cohorts.


# Notes

- Somatic variant calling Mutect must know the SM of the tumor/normal. If we do data merging
  based on sampleName (or sampleMergeKey), and sampleName is _different_ than rgsm. What will
  happen?

- For single samples, Broad using haplotypecaller -ERC GVCF followed by ValidateGVCF?

- duplicate variant calling steps these for exome and genome?

- GATK tools have bug which trims everything from contig name following a final `:`. Is this still
an issue? Are we excluding data by using their intervals for paralellizing tasks?

