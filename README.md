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
We're getting close to being able to easily run the pipeline now, and from this point you might be able to hack your way to make everything run. But it is recommended that you use similar settings to the ones detailed here in order to get the best support possible.  
By running the following command you should be able to see the settings that jetstream is currently using:
```
jetstream settings -v
```
The -v enables a verbose view. The important settings we need to change are slurm and pipelines home, by default they should look similar to this (please note that a bulk amount of settings were omitted from this block, only relevant settings are shown):
```
backend: local
backends:
  local:
    (): jetstream.backends.local.LocalBackend
    blocking_io_penalty: 30
    cpus: null
  slurm:
    (): jetstream.backends.slurm.SlurmBackend
    job_monitor_max_fails: 5
    sacct_fields:
    - JobID
    - JobName
    - Partition
    - Account
    - AllocCPUS
    - State
    - ExitCode
    - Elapsed
    - MaxRSS
    sacct_frequency: 10
.
.
.
pipelines:
  home: null
projects:
  history_filename: '%Y-%m-%dT%H-%M-%SZ.yaml'
  history_max_tries: 10
  id_format: p{id}
  lock_timeout: 60
runner:
  autosave_max: 60
  autosave_min: 5
  id_format: js{id}
  max_concurrency: null
  throttle: 0.1
tasks:
  summary_fields:
  - name
  - state.status
```
We need to change the backend to be slurm for running at TGen and we also need to change the home location of our pipelines to the parent directory of the phoenix pipeline that we downloaded earlier. To do this, we simply need to edit the config.yaml file for jetstream or create the config file if it does not exist already. The location for this file is, by default, located in the .config/jetstream directory of our home directory. The following commands will allow you to find and edit/create this file:
```
cd ~
ls -a
cd .config/jetstream/
vim config.yaml
```
cd ~ <- changes our directory to home 
ls -a <- shows all files present in our directory, including our .config folder hopefully
cd .config/jetstream/ <- change to location of jetstream config.yaml
vim config.yaml <- this will open prompt for editing the config file  
The config.yaml needs to be modified to look similar to:
```
# Jetstream Common User Settings
backend: slurm
pipelines:
  home: /home/tgenjetstream/jetstream_pipelines/:/home/tgenjetstream/jetstream_centro/
```
The home location will differ by install and should be the location of our downloaded pipelines. If you aren't sure where they are downloaded, use
```
cd ~
find . -name pipeline.yaml
```
This may return more than one result.

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


