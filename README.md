# Phoenix

## Install Guide
First, jetstream must be installed, available [here](https://github.com/tgen/jetstream/tree/master). There is an install guide for jetstream available at the link provided, but the gist of the guide is to install using pip3 using a command similar to:  
```bash
pip3 install --upgrade --user git+https://github.com/tgen/jetstream.git@master
```  
Then the recommended install method for installing the phoenix pipeline is:
```
cd ~
mkdir jetstream_pipelines
cd jetstream_pipelines
git clone https://github.com/tgen/phoenix
```
//TODO - still working on install guide  
From here and what I can tell, we should be able to run the pipeline within TGen, outsiders would almost definitely need to get every dependency we use and essentially replicate our structure to the T. I'm still working on ensuring that the pipeline is runnable, and from there I will likely be able to handle configurations to allow the pipeline to work anywhere, or at least give a guide to the bare minimum framework that is needed.

## Features

_Click to show details_

<details>
  <summary><b>DNA Alignment</b></summary>

  The DNA Alignment module takes fastqs to aligned BAMs, then runs BAM qc steps.
  `dnaAlignmentStyle` can be used to select from different alignment strategies.

</details>

<details>
  <summary><b>Somatic Variant Calling</b></summary>
  
  This module will run several somatic variant callers on tumor/normal 
  data pairs:

  - Strelka2 (Somatic mode)
  - Mutect2
  - Lancet

</details>

<details>
  <summary><b>Germline Variant Calling</b></summary>
  
  Generates germline variant call files (VCF) with several callers. Additionally,
  this will create a gVCF for each sample that can be used to jointly call large
  cohorts.

  - GATK HaplotypeCaller (gVCF mode)
  - Freebayes
  - Strelka2 (Germline mode)
  - Octopus (Individual)

</details>


## Config

### Data file attributes

  - *glType* [genome|genomephased|exome|rna|singlecellrna|singlecelldna|matepair|chip]
    
    Used for determining if the sample is DNA/RNA/etc. and adding the corresponding
    tasks to the final workflow. Each sample discovered will take this attribute from
    the first file encountered for that sample in the config file.

  - *glPrep* [genome|capture|rna|singlecellrna|singlecellenrichment|singlecellcdna|singlecelltargetamp|matepair|chip]

    Not used yet.

  - rnaStrandDirection [NotApplicable|]
    rnaStrandType [NotApplicable|]
    readOrientation [Inward|]
    assayCode
    fastqCode [R1|R2]
    fastqPath
    fileType
    numberOfReads
    
    rgcn
    rgid
    rgks
    rglb
    rgpl
    rgpm
    rgpu
    rgsm

    sampleMergeKey
    sampleName
    subGroup

# TGen Naming Convention
- The naming format used in the LIMS trys to ensures all specimens have unique file names based on specific metadata fields

STUDY_PATIENT_VISIT_SOURCE_FRACTION_SubgroupIncrement_ASSAY_LIBRARY<br/>
----------------- PatientID<br/>
---------------------- VisitID<br/>
--------------------------------- SpecimenID<br/>
----------------------------------------------------------------- SampleID<br/>
----------------------------------------------------------------- RG.SM (VCF file genotype column header)<br/>
------------------------------------------------------------------------- sampleMergeKey (BAM filename)<br/>

# sampleName vs RG.SM vs sample_mergeKey

sampleName is used for naming files and merging together multiple input read files (fastqs) into 
a single bam for the sample. RG SM is very similar but intentionally left separate to support 
several "samples" from a single "specimen". 

#### RG.SM

`STUDY_1234_1_PB_Whole_C1`
- This is the value we expect to be in the header of any VCF genotype column related to sample in question
- Allows rapid comparison between assays for the same patient

#### sampleMergeKey

`STUDY_1234_1_PB_Whole_C1_KHWGS`
- This key is used to merge data from multiple sequencing lanes or flowcells for data from the same specimen tested with the same assay
- This is the expected BAM filename

#### sampleName

`STUDY_1234_1_PB_Whole_C1_KHWGS_K12345`
- This is the expected base FASTQ filename


