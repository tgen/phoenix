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
jetstream settings -c -b "slurm" -P "/home/USERNAME/jetstream_pipelines/"
```  
Note that USERNAME is replaced by your username, e.g. "/home/jsmith/jetstream_pipelines/" for a user name John Smith within TGen. Also, the home location will differ by install and should be the location of our downloaded pipelines. If you aren't sure where they are downloaded, use
```
cd ~
find . -name pipeline.yaml
```
This may return more than one result. The one we are looking for should look similar to "./jetstream_pipelines/phoenix/pipeline.yaml". The ./ means that we have jetstream_pipelines in our current directory. Which means that the true path to our pipelines is /home/USER/jetstream_pipeline/. Note that USER is your username, within TGen this is generally your first initial and then last name. Like the John Smith example above.  

Congratulations! That's it! We now have jetstream and the phoenix pipeline installed. In the next section we will discuss how we can actually run the pipeline from the command line.  

## Running from command line
Now that we have phoenix installed and ready to run, we need to setup a project to run. To do so, we need to create a config file for our project. The general format is as follows:
```
{
  "project": "",
  "study": "",
  "email": "",
  "hpcAccount": "",
  "isilonPath": "",
  "pipeline": "phoenix",
  "dataFiles": []
}
```

Here is a larger example with actual data for running the phoenix pipeline on a NA12878 project:
<details><summary>NA12878 Example</summary>
  <p>
    
  ##### Some of this data has been modified to hide the identity of the original submitter(s)
  ```bash
  {
    "cram": true,
    "dataFiles": [
        {
            "assayCode": "TPFWG",
            "dnaRnaMergeKey": "GIAB_NA12878_1_CL_Whole",
            "fastqCode": "R1",
            "fastqPath": "/home/tgenref/homo_sapiens/control_files/giab/fastq/NA12878_140407_D00360_0016_ASUPERFQS01/Project_GIAB_NA12878_1_TPFWG/Sample_GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS01/GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS01_NoIndex_L001_R1_001.fastq.gz",
            "fileType": "fastq",
            "fraction": "Whole",
            "glPrep": "Genome",
            "glType": "Genome",
            "index1Length": 6,
            "index2Length": 0,
            "limsLibraryRecordId": 64391,
            "numberOfReads": 228228468,
            "read1Length": 148,
            "read2Length": 148,
            "readOrientation": "Inward",
            "rgcn": "TGen",
            "rgid": "SUPERFQS01_1_K18088",
            "rgks": "ATCACG",
            "rglb": "K18088",
            "rgpl": "ILLUMINA",
            "rgpm": "HiSeq2500",
            "rgpu": "SUPERFQS01_1",
            "rgsm": "GIAB_NA12878_1_CL_Whole_C1",
            "rnaStrandDirection": "NotApplicable",
            "rnaStrandType": "NotApplicable",
            "sampleMergeKey": "GIAB_NA12878_1_CL_Whole_C1_TPFWG",
            "sampleName": "GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088",
            "subGroup": "Constitutional",
            "umiInLine": "false",
            "umiLength": 0,
            "umiRead": false
        },
        {
            "assayCode": "TPFWG",
            "dnaRnaMergeKey": "GIAB_NA12878_1_CL_Whole",
            "fastqCode": "R2",
            "fastqPath": "/home/tgenref/homo_sapiens/control_files/giab/fastq/NA12878_140407_D00360_0016_ASUPERFQS01/Project_GIAB_NA12878_1_TPFWG/Sample_GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS01/GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS01_NoIndex_L001_R2_001.fastq.gz",
            "fileType": "fastq",
            "fraction": "Whole",
            "glPrep": "Genome",
            "glType": "Genome",
            "index1Length": 6,
            "index2Length": 0,
            "limsLibraryRecordId": 64391,
            "numberOfReads": 228228468,
            "read1Length": 148,
            "read2Length": 148,
            "readOrientation": "Inward",
            "rgcn": "TGen",
            "rgid": "SUPERFQS01_1_K18088",
            "rgks": "ATCACG",
            "rglb": "K18088",
            "rgpl": "ILLUMINA",
            "rgpm": "HiSeq2500",
            "rgpu": "SUPERFQS01_1",
            "rgsm": "GIAB_NA12878_1_CL_Whole_C1",
            "rnaStrandDirection": "NotApplicable",
            "rnaStrandType": "NotApplicable",
            "sampleMergeKey": "GIAB_NA12878_1_CL_Whole_C1_TPFWG",
            "sampleName": "GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088",
            "subGroup": "Constitutional",
            "umiInLine": "false",
            "umiLength": 0,
            "umiRead": false
        },
        {
            "assayCode": "TPFWG",
            "dnaRnaMergeKey": "GIAB_NA12878_1_CL_Whole",
            "fastqCode": "R1",
            "fastqPath": "/home/tgenref/homo_sapiens/control_files/giab/fastq/NA12878_140407_D00360_0016_ASUPERFQS01/Project_GIAB_NA12878_1_TPFWG/Sample_GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS01/GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS01_NoIndex_L002_R1_001.fastq.gz",
            "fileType": "fastq",
            "fraction": "Whole",
            "glPrep": "Genome",
            "glType": "Genome",
            "index1Length": 6,
            "index2Length": 0,
            "limsLibraryRecordId": 64391,
            "numberOfReads": 232329584,
            "read1Length": 148,
            "read2Length": 148,
            "readOrientation": "Inward",
            "rgcn": "TGen",
            "rgid": "SUPERFQS01_2_K18088",
            "rgks": "ATCACG",
            "rglb": "K18088",
            "rgpl": "ILLUMINA",
            "rgpm": "HiSeq2500",
            "rgpu": "SUPERFQS01_2",
            "rgsm": "GIAB_NA12878_1_CL_Whole_C1",
            "rnaStrandDirection": "NotApplicable",
            "rnaStrandType": "NotApplicable",
            "sampleMergeKey": "GIAB_NA12878_1_CL_Whole_C1_TPFWG",
            "sampleName": "GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088",
            "subGroup": "Constitutional",
            "umiInLine": "false",
            "umiLength": 0,
            "umiRead": false
        },
        {
            "assayCode": "TPFWG",
            "dnaRnaMergeKey": "GIAB_NA12878_1_CL_Whole",
            "fastqCode": "R2",
            "fastqPath": "/home/tgenref/homo_sapiens/control_files/giab/fastq/NA12878_140407_D00360_0016_ASUPERFQS01/Project_GIAB_NA12878_1_TPFWG/Sample_GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS01/GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS01_NoIndex_L002_R2_001.fastq.gz",
            "fileType": "fastq",
            "fraction": "Whole",
            "glPrep": "Genome",
            "glType": "Genome",
            "index1Length": 6,
            "index2Length": 0,
            "limsLibraryRecordId": 64391,
            "numberOfReads": 232329584,
            "read1Length": 148,
            "read2Length": 148,
            "readOrientation": "Inward",
            "rgcn": "TGen",
            "rgid": "SUPERFQS01_2_K18088",
            "rgks": "ATCACG",
            "rglb": "K18088",
            "rgpl": "ILLUMINA",
            "rgpm": "HiSeq2500",
            "rgpu": "SUPERFQS01_2",
            "rgsm": "GIAB_NA12878_1_CL_Whole_C1",
            "rnaStrandDirection": "NotApplicable",
            "rnaStrandType": "NotApplicable",
            "sampleMergeKey": "GIAB_NA12878_1_CL_Whole_C1_TPFWG",
            "sampleName": "GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088",
            "subGroup": "Constitutional",
            "umiInLine": "false",
            "umiLength": 0,
            "umiRead": false
        },
        {
            "assayCode": "TPFWG",
            "dnaRnaMergeKey": "GIAB_NA12878_1_CL_Whole",
            "fastqCode": "R1",
            "fastqPath": "/home/tgenref/homo_sapiens/control_files/giab/fastq/NA12878_140407_D00360_0017_BSUPERFQS02/Project_GIAB_NA12878_1_TPFWG/Sample_GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS02/GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS02_NoIndex_L001_R1_001.fastq.gz",
            "fileType": "fastq",
            "fraction": "Whole",
            "glPrep": "Genome",
            "glType": "Genome",
            "index1Length": 6,
            "index2Length": 0,
            "limsLibraryRecordId": 64391,
            "numberOfReads": 231481708,
            "read1Length": 148,
            "read2Length": 148,
            "readOrientation": "Inward",
            "rgcn": "TGen",
            "rgid": "SUPERFQS02_1_K18088",
            "rgks": "ATCACG",
            "rglb": "K18088",
            "rgpl": "ILLUMINA",
            "rgpm": "HiSeq2500",
            "rgpu": "SUPERFQS02_1",
            "rgsm": "GIAB_NA12878_1_CL_Whole_C1",
            "rnaStrandDirection": "NotApplicable",
            "rnaStrandType": "NotApplicable",
            "sampleMergeKey": "GIAB_NA12878_1_CL_Whole_C1_TPFWG",
            "sampleName": "GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088",
            "subGroup": "Constitutional",
            "umiInLine": "false",
            "umiLength": 0,
            "umiRead": false
        },
        {
            "assayCode": "TPFWG",
            "dnaRnaMergeKey": "GIAB_NA12878_1_CL_Whole",
            "fastqCode": "R2",
            "fastqPath": "/home/tgenref/homo_sapiens/control_files/giab/fastq/NA12878_140407_D00360_0017_BSUPERFQS02/Project_GIAB_NA12878_1_TPFWG/Sample_GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS02/GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS02_NoIndex_L001_R2_001.fastq.gz",
            "fileType": "fastq",
            "fraction": "Whole",
            "glPrep": "Genome",
            "glType": "Genome",
            "index1Length": 6,
            "index2Length": 0,
            "limsLibraryRecordId": 64391,
            "numberOfReads": 231481708,
            "read1Length": 148,
            "read2Length": 148,
            "readOrientation": "Inward",
            "rgcn": "TGen",
            "rgid": "SUPERFQS02_1_K18088",
            "rgks": "ATCACG",
            "rglb": "K18088",
            "rgpl": "ILLUMINA",
            "rgpm": "HiSeq2500",
            "rgpu": "SUPERFQS02_1",
            "rgsm": "GIAB_NA12878_1_CL_Whole_C1",
            "rnaStrandDirection": "NotApplicable",
            "rnaStrandType": "NotApplicable",
            "sampleMergeKey": "GIAB_NA12878_1_CL_Whole_C1_TPFWG",
            "sampleName": "GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088",
            "subGroup": "Constitutional",
            "umiInLine": "false",
            "umiLength": 0,
            "umiRead": false
        },
        {
            "assayCode": "TPFWG",
            "dnaRnaMergeKey": "GIAB_NA12878_1_CL_Whole",
            "fastqCode": "R1",
            "fastqPath": "/home/tgenref/homo_sapiens/control_files/giab/fastq/NA12878_140407_D00360_0017_BSUPERFQS02/Project_GIAB_NA12878_1_TPFWG/Sample_GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS02/GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS02_NoIndex_L002_R1_001.fastq.gz",
            "fileType": "fastq",
            "fraction": "Whole",
            "glPrep": "Genome",
            "glType": "Genome",
            "index1Length": 6,
            "index2Length": 0,
            "limsLibraryRecordId": 64391,
            "numberOfReads": 228988006,
            "read1Length": 148,
            "read2Length": 148,
            "readOrientation": "Inward",
            "rgcn": "TGen",
            "rgid": "SUPERFQS02_2_K18088",
            "rgks": "ATCACG",
            "rglb": "K18088",
            "rgpl": "ILLUMINA",
            "rgpm": "HiSeq2500",
            "rgpu": "SUPERFQS02_2",
            "rgsm": "GIAB_NA12878_1_CL_Whole_C1",
            "rnaStrandDirection": "NotApplicable",
            "rnaStrandType": "NotApplicable",
            "sampleMergeKey": "GIAB_NA12878_1_CL_Whole_C1_TPFWG",
            "sampleName": "GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088",
            "subGroup": "Constitutional",
            "umiInLine": "false",
            "umiLength": 0,
            "umiRead": false
        },
        {
            "assayCode": "TPFWG",
            "dnaRnaMergeKey": "GIAB_NA12878_1_CL_Whole",
            "fastqCode": "R2",
            "fastqPath": "/home/tgenref/homo_sapiens/control_files/giab/fastq/NA12878_140407_D00360_0017_BSUPERFQS02/Project_GIAB_NA12878_1_TPFWG/Sample_GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS02/GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088_SUPERFQS02_NoIndex_L002_R2_001.fastq.gz",
            "fileType": "fastq",
            "fraction": "Whole",
            "glPrep": "Genome",
            "glType": "Genome",
            "index1Length": 6,
            "index2Length": 0,
            "limsLibraryRecordId": 64391,
            "numberOfReads": 228988006,
            "read1Length": 148,
            "read2Length": 148,
            "readOrientation": "Inward",
            "rgcn": "TGen",
            "rgid": "SUPERFQS02_2_K18088",
            "rgks": "ATCACG",
            "rglb": "K18088",
            "rgpl": "ILLUMINA",
            "rgpm": "HiSeq2500",
            "rgpu": "SUPERFQS02_2",
            "rgsm": "GIAB_NA12878_1_CL_Whole_C1",
            "rnaStrandDirection": "NotApplicable",
            "rnaStrandType": "NotApplicable",
            "sampleMergeKey": "GIAB_NA12878_1_CL_Whole_C1_TPFWG",
            "sampleName": "GIAB_NA12878_1_CL_Whole_C1_TPFWG_K18088",
            "subGroup": "Constitutional",
            "umiInLine": "false",
            "umiLength": 0,
            "umiRead": false
        }
    ],
    "dnaAlignmentStyle": "tgen",
    "email": "example@tgen.org",
    "ethnicity": "Caucasian",
    "familyCode": "",
    "holdConfig": false,
    "hpcAccount": "tgen-#####",
    "isilonPath": "/example/giab/",
    "matchedNormal": true,
    "matchedNormalToUse": "",
    "maternalID": "",
    "patCode": "NA12878",
    "paternalID": "",
    "pipeline": "phoenix",
    "project": "GIAB_NA12878",
    "sex": "Female",
    "study": "GIAB",
    "submitter": "euser",
    "submitterEmail": "examplet@tgen.org",
    "tasks": {
        "Exome_base_recalibration_gatk": false,
        "Exome_deepvariant_deepvariant": false,
        "Exome_freebayes_freebayes": false,
        "Exome_genotype_hc_gvcf_gatk": false,
        "Exome_hc_gvcf_gatk": true,
        "Exome_hmmcopy_make_wig_HMM_Copy_Utils": true,
        "Exome_ichor_cna_ichor": true,
        "Exome_mark_duplicates_gatk": false,
        "Exome_mark_duplicates_samblaster": false,
        "Exome_mark_duplicates_samtools": true,
        "Exome_mutect2_gatk": true,
        "Exome_octopus_octopus": false,
        "Exome_strelka2_strelka2": true,
        "Genome_base_recalibration_gatk": false,
        "Genome_deepvariant_deepvariant": false,
        "Genome_freebayes_freebayes": false,
        "Genome_genotype_hc_gvcf_gatk": false,
        "Genome_hc_gvcf_gatk": true,
        "Genome_hmmcopy_make_wig_HMM_Copy_Utils": true,
        "Genome_ichor_cna_ichor": true,
        "Genome_mark_duplicates_gatk": false,
        "Genome_mark_duplicates_samblaster": false,
        "Genome_mark_duplicates_samtools": true,
        "Genome_mutect2_gatk": true,
        "Genome_octopus_octopus": false,
        "Genome_strelka2_strelka2": true
    },
    "varDB": false
}
```
  </p>
</details>


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


