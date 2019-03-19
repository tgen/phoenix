# Phoenix

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

  - rnaStrandDirection 
    rnaStrandType
    readOrientation

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
