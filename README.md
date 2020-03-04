# Phoenix Workflow - A JetStream Workflow

This workflow supports the analysis of human sequencing samples against the GRCh38 reference genome using the ensembl version 98 gene models

## Supported Analysis

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
    - We do NOT run it on genomes 
  - VarDictJava
  - Octopus (Only calls on primary contigs 1-22, X, Y)

</details>

<details>
  <summary><b>Constitutional Variant Calling</b></summary>
  
  Generates constitutional variant call files (VCF) with several callers. Additionally,
  this will create a gVCF for each sample that can be used to jointly call large
  cohorts.

  - Deepvariant (Only calls on primary contigs 1-22, X, Y)
  - GATK HaplotypeCaller (gVCF mode)
    - single sample calling using genotypeGVCFs and CNNscroreVariant
  - Freebayes
  - Strelka2 (Germline mode)
  - Octopus (Individual) (Only calls on primary contigs 1-22, X, Y)

</details>

<details>
  <summary><b>RNA Alignment</b></summary>

  - STAR

</details>

<details>
  <summary><b>Gene Expression Estimates</b></summary>

  - Salmon
  - Star-count - Recommended count method
  - HtSeq
  - FeatureCounts

</details>

<details>
  <summary><b>Fusion Transcript Detection</b></summary>

  - Salmon
  - Star-count - Recommended count method
  - HtSeq
  - FeatureCounts

</details>

<details>
  <summary><b>Single Cell Sequencing</b></summary>
  
  We support multiple things
    
  <details>
    <summary><b>scRNAseq</b></summary>
    XXX
  </details>
  
  <details>
    <summary><b>scDNAseq</b></summary>
    YYY
  </details>

</details>

### Required Software
_Click to show details_

<details>
<summary><b>Public Tools Used by the Workflow</b></summary>  

All tools build with public easybuild configuration files, available [here](https://github.com/easybuilders/easybuild).<br/>
*Last Updated March 4th, 2020*  

| Tool | Version Implemented | Current Version | Dependency and Notes | EasyBuild |
| :---: | :---: | :---: | :--- | :---: |
| [bcftools](https://github.com/samtools/bcftools/releases) | 1.10.1 | **1.10.2** | | Yes |
| [bedtools](https://github.com/arq5x/bedtools2/releases) | 2.29.0 | **2.29.2** | delly-filter, addmatchRNA, vardict, vcfmerger2 | Yes |
| [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) | 2.3.5.1 | **2.4.1** | star-fusion | Yes |
| [bwa](https://github.com/lh3/bwa/releases) | 0.7.17 | 0.7.17 | | Yes |
| [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) | 3.1.0 | 3.1.0 | | Yes |
| [deepvariant](https://github.com/google/deepvariant/releases) | 0.9.0 | 0.9.0 | singularity container | Yes |
| [delly](https://github.com/dellytools/delly/releases) | 0.7.6 | **0.8.2** | staying with 0.7.6 for compatibility reasons | Yes |
| [freebayes](https://github.com/ekg/freebayes/releases) | 1.3.1 | **1.3.2** | 1.3.2 ensures python3 compatibility | Yes |
| [gatk](https://github.com/broadinstitute/gatk//releases) | 4.1.4.0 | **4.1.5.0** |  | Yes |
| [gmap-gsnp](http://research-pub.gene.com/gmap/) | 2019-09-12 | 2019-09-12 | star-fusion | Yes |
| [gridss](https://github.com/PapenfussLab/gridss/releases) | 2.6.3 | **2.8.0** |  | |
| [hmmcopyutils](https://github.com/shahcompbio/hmmcopy_utils) | 1.0 | 1.0 | no official release | Yes |
| [htseq](https://github.com/simon-anders/htseq/releases) | 0.11.1 | **0.11.3** | | |
| [htslib](https://github.com/samtools/htslib/releases) | 1.10.1 | **1.10.2** | star-fusion(bgzip) | Yes |
| [ichor](https://github.com/broadinstitute/ichorCNA/releases) | 0.2.0 | 0.2.0 | package in R/3.6.1-phoenix module | Yes? |
| [jellyfish](https://github.com/gmarcais/Jellyfish/releases) | 2.3.0 | 2.3.0 | star-fusion | Yes |
| [lancet](https://github.com/nygenome/lancet/releases) | 1.1.0 | 1.1.0 | | Yes |
| [manta](https://github.com/Illumina/manta/releases) | 1.6.0 | 1.6.0 | | Yes |
| [multiQC](https://github.com/ewels/MultiQC/releases) | 1.7 | **1.8** | python3 pip | Yes |
| [octopus](https://github.com/luntergroup/octopus/releases) | 0.6.3-beta | 0.6.3-beta | | Yes |
| [pairoscope](https://github.com/genome/pairoscope/releases) | 0.4.2 | 0.4.2 | | Yes |
| [perl](https://github.com/Perl/perl5/releases) | 5.28.1 | **5.30.2-RC1** | star-fusion | Yes |
| [phaser](https://github.com/secastel/phaser/tree/master/phaser) | 1.1.1 | 1.1.1 | vcfmerger2 | Yes |
| [python2](https://www.python.org/downloads/) | 2.7.15 | **2.7.17** | | Yes |
| [python3](https://www.python.org/downloads/) | 3.7.2 | **3.8.2** | star-fusion, vcfmerger2 | Yes |
| [R](https://www.r-project.org/) | 3.6.1 | **3.6.3** | gatk cnv, varDict, vcfmerger2 | Yes |
| [sambamba](https://github.com/biod/sambamba/releases) | 0.7.0 | **0.7.1** | | |
| [samblaster](https://github.com/GregoryFaust/samblaster/releases) | 0.1.24 | 0.1.24 | | |
| [salmon](https://github.com/COMBINE-lab/salmon/releases) | 0.14.1 | **1.1.0** | self, star-fusion | Yes |
| [samtools](https://github.com/samtools/samtools/releases) | 1.10 | 1.10 | markdup needs unreleased patch | Yes |
| [singularity](https://github.com/sylabs/singularity/releases) | 3.5.2 | **3.5.3** | deepvariant | Yes, release-candidate |
| [snpEff](https://sourceforge.net/projects/snpeff/files/) | 4.3t | 4.3t | no updates since Nov 2017| Yes |
| [star](https://github.com/alexdobin/STAR/releases) | 2.7.3a | 2.7.3a | self, star-fusion | Yes |
| [star-fusion](https://github.com/STAR-Fusion/STAR-Fusion/releases) | 1.8.1 | 1.8.1 | | Yes |
| [strelka](https://github.com/Illumina/strelka/releases) | 2.9.10 | 2.9.10 | | Yes |
| [subread](https://sourceforge.net/projects/subread/) | 2.0.0 | 2.0.0 | part of subread package | Yes |
| [trinityrnaseq](https://github.com/trinityrnaseq/trinityrnaseq/releases) | 2.8.6 | **2.9.1** | star-fusion | Yes |
| [vardictJava](https://github.com/AstraZeneca-NGS/VarDictJava/releases) | 1.7.0 | 1.7.0 | | Yes |
| [vcfmerger2](https://github.com/tgen/vcfMerger2/releases) | 0.8.7 | 0.8.7 | | Yes |
| [vep](https://github.com/Ensembl/ensembl-vep/releases) | 98.3 | **99.2** | | Yes |
| [verifybamid2](https://github.com/Griffan/VerifyBamID/releases) | 1.0.6 | 1.0.6 | | Yes |
| [vt](https://github.com/atks/vt/releases) | 0_57721 | 0_57721 | | Yes |

</details>

<details>
<summary><b>Required PERL Modules</b></summary>

|	Module	|	Version	|
|	:---:	|	:---:	|
|	AutoLoader	|	5.74	|
|	Bio::EnsEMBL::XS	|	2.3.2	|
|	Carp::Assert	|	0.21	|
|	common::sense	|	3.74	|
|	constant	|	1.33	|
|	Cwd	|	3.75	|
|	Data::Dumper	|	2.173	|
|	Data::Dumper	|	2.173	|
|	DB_File	|	1.852	|
|	Encode	|	3.01	|
|	ExtUtils::CBuilder	|	0.280231	|
|	ExtUtils::CppGuess	|	0.19	|
|	ExtUtils::MakeMaker	|	7.36	|
|	ExtUtils::ParseXS	|	3.35	|
|	IO::File	|	1.39	|
|	IO::Select	|	1.39	|
|	IO::Socket	|	1.39	|
|	IPC::Cmd	|	1.02	|
|	JSON::XS	|	4.02	|
|	List::Util	|	1.5	|
|	Locale::Maketext::Simple	|	0.21	|
|	MIME::Base64	|	3.15	|
|	Module::CoreList	|	5.2019062	|
|	Module::Load	|	0.34	|
|	Module::Load::Conditional	|	0.68	|
|	Module::Metadata	|	1.000036	|
|	Net::Domain	|	3.11	|
|	Params::Check	|	0.38	|
|	parent	|	0.237	|
|	Perl::OSType	|	1.01	|
|	PerlIO::gzip	|	0.2	|
|	Pod::Escapes	|	1.07	|
|	Pod::Man	|	4.12	|
|	Pod::Simple	|	3.39	|
|	Scalar::Util	|	1.5	|
|	Set::IntervalTree	|	0.12	|
|	Socket	|	2.029	|
|	Storable	|	3.15	|
|	Test	|	1.26	|
|	Test::More	|	1.302164	|
|	Text::ParseWords	|	3.3	|
|	Text::Wrap	|	2013.0523	|
|	Time::Local	|	1.28	|
|	Types::Serialiser	|	1	|
|	URI::Escape	|	1.76	|
|	XSLoader	|	0.24	|

</details>

<details>
<summary><b>Required PYTHON2 Modules</b></summary>

|	Module	|	Version	|
|	:---:	|	:---:	|
|	alabaster	|	0.7.12	|
|	asn1crypto	|	0.24.0	|
|	atomicwrites	|	1.3.0	|
|	attrs	|	19.1.0	|
|	Babel	|	2.6.0	|
|	bcrypt	|	3.1.6	|
|	bitstring	|	3.1.5	|
|	blist	|	1.3.6	|
|	certifi	|	2019.3.9	|
|	cffi	|	1.12.2	|
|	chardet	|	3.0.4	|
|	Click	|	7	|
|	cryptography	|	2.6.1	|
|	Cython	|	0.29.6	|
|	deap	|	1.2.2	|
|	decorator	|	4.3.2	|
|	docopt	|	0.6.2	|
|	docutils	|	0.14	|
|	ecdsa	|	0.13	|
|	enum34	|	1.1.6	|
|	funcsigs	|	1.0.2	|
|	future	|	0.17.1	|
|	idna	|	2.8	|
|	imagesize	|	1.1.0	|
|	ipaddress	|	1.0.22	|
|	Jinja2	|	2.1	|
|	joblib	|	0.13.2	|
|	liac-arff	|	2.4.0	|
|	MarkupSafe	|	1.1.1	|
|	mock	|	2.0.0	|
|	more-itertools	|	5.0.0	|
|	netaddr	|	0.7.19	|
|	netifaces	|	0.10.9	|
|	nose	|	1.3.7	|
|	packaging	|	19	|
|	paramiko	|	2.4.2	|
|	pathlib2	|	2.3.3	|
|	paycheck	|	1.0.2	|
|	pbr	|	5.1.3	|
|	pip	|	19.0.3	|
|	pluggy	|	0.9.0	|
|	psutil	|	5.6.1	|
|	py	|	1.8.0	|
|	py_expression_eval	|	0.3.6	|
|	pyasn1	|	0.4.5	|
|	pycparser	|	2.19	|
|	pycrypto	|	2.6.1	|
|	Pygments	|	2.3.1	|
|	PyNaCl	|	1.3.0	|
|	pyparsing	|	2.3.1	|
|	pytest	|	4.3.1	|
|	python-dateutil	|	2.8.0	|
|	pytz	|	2018.9	|
|	requests	|	2.21.0	|
|	scandir	|	1.10.0	|
|	setuptools	|	40.8.0	|
|	setuptools_scm	|	3.2.0	|
|	singledispatch	|	3.4.0.3	|
|	six	|	1.12.0	|
|	snowballstemmer	|	1.2.1	|
|	Sphinx	|	1.8.5	|
|	sphinxcontrib-websupport	|	1.1.0	|
|	tabulate	|	0.8.3	|
|	typing	|	3.6.6	|
|	ujson	|	1.35	|
|	urllib3	|	1.24.1	|
|	virtualenv	|	16.4.3	|
|	wheel	|	0.33.1	|
|	xlrd	|	1.2.0	|

</details>

<details>
<summary><b>Required PYTHON3 Modules</b></summary>

* cyvcf2
* getopt
* igv_reports
* multiprocessing
* multiqc
* PIL
* pybedtools
* pysam
* requests
* sys (argv, exit)

</details>

<details>
<summary><b>Required R Libraries</b></summary>

* ichorCNA2.0
* tidyverse (ggplot)

</details>

## Install Guide
Please see the [wiki](https://github.com/tgen/phoenix/wiki) for a detailed install guide

## Running from command line
In order to run from the command line, we need to create a config file for our project. The general format is as follows:
```json
{
    "project": "",
    "study": "",
    "email": "",
    "hpcAccount": "",
    "isilonPath": "",
    "pipeline": "phoenix",
    "dataFiles": [],
    "dnaAlignmentStyle": "tgen",
    "email": "somebody@tgen.org",
    "isilonPath": "",
    "pipeline": "phoenix@version",
    "project": "Project_Name",
    "submitter": "somebody",
    "tasks": {},
}
```

Here is a larger example with actual data for running the phoenix pipeline on a NA12878 project:
<details><summary>NA12878 Example</summary>
  <p>
    
  ##### Some of this data has been modified to hide the identity of the original submitter(s)
  ```json
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
    "submitter": "user",
    "submitterEmail": "examplet@tgen.org",
    "tasks": {
        "Exome_alignment_base_recalibration_gatk": false,
        "Exome_alignment_mark_duplicates_gatk": false,
        "Exome_alignment_mark_duplicates_samblaster": false,
        "Exome_alignment_mark_duplicates_samtools": true,
        "Exome_constitutional_annotate_vcfs_bcftools_clinvar_20190715": true,
        "Exome_constitutional_annotate_vcfs_bcftools_cosmic_coding_v89": true,
        "Exome_constitutional_annotate_vcfs_bcftools_cosmic_noncoding_v89": true,
        "Exome_constitutional_annotate_vcfs_bcftools_dbsnp_v152": true,
        "Exome_constitutional_annotate_vcfs_bcftools_gnomad_exome_v2_1_1_liftover": true,
        "Exome_constitutional_annotate_vcfs_bcftools_gnomad_genome_v3_0": true,
        "Exome_constitutional_annotate_vcfs_snpEff_ann": true,
        "Exome_constitutional_annotate_vcfs_vep": true,
        "Exome_constitutional_genotype_hc_gvcf_gatk_GenotypeGVCFs": false,
        "Exome_constitutional_snp_indel_caller_deepvariant": true,
        "Exome_constitutional_snp_indel_caller_freebayes": false,
        "Exome_constitutional_snp_indel_caller_gatk_HaplotypeCaller": true,
        "Exome_constitutional_snp_indel_caller_octopus": false,
        "Exome_constitutional_snp_indel_caller_strelka2": false,
        "Exome_constitutional_structural_caller_manta": true,
        "Exome_quality_control_constitutional_contamination_check_VerifyBamID": true,
        "Exome_quality_control_constitutional_sex_check_freebayes": true,
        "Exome_quality_control_genotype_concordance_snpSniffer": true,
        "Exome_quality_control_stats_gatk_CollectHsMetrics": true,
        "Exome_quality_control_stats_gatk_CollectMultipleMetrics": true,
        "Exome_quality_control_stats_gatk_ConvertSequencingArtifactToOxoG": true,
        "Exome_quality_control_stats_samtools_flagstat": true,
        "Exome_quality_control_stats_samtools_idxstats": true,
        "Exome_quality_control_stats_samtools_stats": true,
        "Exome_somatic_annotate_vcfs_bcftools_clinvar_20190715": true,
        "Exome_somatic_annotate_vcfs_bcftools_cosmic_coding_v89": true,
        "Exome_somatic_annotate_vcfs_bcftools_cosmic_noncoding_v89": true,
        "Exome_somatic_annotate_vcfs_bcftools_dbsnp_v152": true,
        "Exome_somatic_annotate_vcfs_bcftools_gnomad_exome_v2_1_1_liftover": true,
        "Exome_somatic_annotate_vcfs_bcftools_gnomad_genome_v3_0": true,
        "Exome_somatic_annotate_vcfs_snpEff_ann": true,
        "Exome_somatic_annotate_vcfs_vep": true,
        "Exome_somatic_cna_caller_gatk": true,
        "Exome_somatic_merge_vcfs_vcfMerger2": true,
        "Exome_somatic_snp_indel_caller_VarDict": true,
        "Exome_somatic_snp_indel_caller_gatk_mutect2": true,
        "Exome_somatic_snp_indel_caller_lancet": true,
        "Exome_somatic_snp_indel_caller_octopus": true,
        "Exome_somatic_snp_indel_caller_strelka2": true,
        "Exome_somatic_structural_caller_delly": false,
        "Exome_somatic_structural_caller_manta": true,
        "Exome_tumor_mm_igtx_pairoscope": true,
        "Genome_alignment_base_recalibration_gatk": false,
        "Genome_alignment_mark_duplicates_gatk": false,
        "Genome_alignment_mark_duplicates_samblaster": false,
        "Genome_alignment_mark_duplicates_samtools": true,
        "Genome_constitutional_annotate_vcfs_bcftools_clinvar_20190715": true,
        "Genome_constitutional_annotate_vcfs_bcftools_cosmic_coding_v89": true,
        "Genome_constitutional_annotate_vcfs_bcftools_cosmic_noncoding_v89": true,
        "Genome_constitutional_annotate_vcfs_bcftools_dbsnp_v152": true,
        "Genome_constitutional_annotate_vcfs_bcftools_gnomad_exome_v2_1_1_liftover": true,
        "Genome_constitutional_annotate_vcfs_bcftools_gnomad_genome_v3_0": true,
        "Genome_constitutional_annotate_vcfs_snpEff_ann": true,
        "Genome_constitutional_annotate_vcfs_vep": true,
        "Genome_constitutional_cna_caller_ichor": true,
        "Genome_constitutional_genotype_hc_gvcf_gatk_GenotypeGVCFs": false,
        "Genome_constitutional_snp_indel_caller_deepvariant": true,
        "Genome_constitutional_snp_indel_caller_freebayes": false,
        "Genome_constitutional_snp_indel_caller_gatk_HaplotypeCaller": true,
        "Genome_constitutional_snp_indel_caller_octopus": false,
        "Genome_constitutional_snp_indel_caller_strelka2": false,
        "Genome_constitutional_structural_caller_manta": true,
        "Genome_quality_control_constitutional_contamination_check_VerifyBamID": true,
        "Genome_quality_control_constitutional_sex_check_freebayes": true,
        "Genome_quality_control_genotype_concordance_snpSniffer": true,
        "Genome_quality_control_stats_gatk_CollectMultipleMetrics": true,
        "Genome_quality_control_stats_gatk_CollectRawWgsMetrics": true,
        "Genome_quality_control_stats_gatk_CollectWgsMetrics": true,
        "Genome_quality_control_stats_gatk_CollectWgsMetricsWithNonZeroCoverage": true,
        "Genome_quality_control_stats_gatk_ConvertSequencingArtifactToOxoG": true,
        "Genome_quality_control_stats_samtools_flagstat": true,
        "Genome_quality_control_stats_samtools_idxstats": true,
        "Genome_quality_control_stats_samtools_stats": true,
        "Genome_somatic_annotate_vcfs_bcftools_clinvar_20190715": true,
        "Genome_somatic_annotate_vcfs_bcftools_cosmic_coding_v89": true,
        "Genome_somatic_annotate_vcfs_bcftools_cosmic_noncoding_v89": true,
        "Genome_somatic_annotate_vcfs_bcftools_dbsnp_v152": true,
        "Genome_somatic_annotate_vcfs_bcftools_gnomad_exome_v2_1_1_liftover": true,
        "Genome_somatic_annotate_vcfs_bcftools_gnomad_genome_v3_0": true,
        "Genome_somatic_annotate_vcfs_snpEff_ann": true,
        "Genome_somatic_annotate_vcfs_vep": true,
        "Genome_somatic_cna_caller_gatk": true,
        "Genome_somatic_merge_vcfs_vcfMerger2": true,
        "Genome_somatic_snp_indel_caller_VarDict": true,
        "Genome_somatic_snp_indel_caller_gatk_mutect2": true,
        "Genome_somatic_snp_indel_caller_lancet": false,
        "Genome_somatic_snp_indel_caller_octopus": true,
        "Genome_somatic_snp_indel_caller_strelka2": true,
        "Genome_somatic_structural_caller_delly": true,
        "Genome_somatic_structural_caller_manta": true,
        "Genome_tumor_mm_igtx_pairoscope": true,
        "RNA_alignment_rna_alignment_STAR": true,
        "RNA_quality_control_genotype_concordance_snpSniffer": true,
        "RNA_quality_control_stats_gatk_CollectMultipleMetrics": true,
        "RNA_quality_control_stats_gatk_CollectRnaSeqMetrics": true,
        "RNA_quality_control_stats_gatk_ConvertSequencingArtifactToOxoG": true,
        "RNA_quality_control_stats_samtools_bedcov": true,
        "RNA_quality_control_stats_samtools_flagstat": true,
        "RNA_quality_control_stats_samtools_idxstats": true,
        "RNA_quality_control_stats_samtools_stats": true,
        "RNA_transcriptome_fusion_caller_Keep_STAR_Fusion_BAM": false,
        "RNA_transcriptome_fusion_caller_STAR_Fusion": true,
        "RNA_transcriptome_quantify_expression_HTseq": false,
        "RNA_transcriptome_quantify_expression_RSEM": false,
        "RNA_transcriptome_quantify_expression_salmon_bam": false,
        "RNA_transcriptome_quantify_expression_salmon_fastqs": true,
        "RNA_transcriptome_quantify_expression_subread_featureCounts": false,
        "SingleCellRNA_VDJ_Assembly_cellranger_vdj": true,
        "SingleCellRNA_transcriptome_quantify_expression_STARsolo": true,
        "SingleCellRNA_transcriptome_quantify_expression_cellranger_count": true
    },
    "varDB": false
}
```

Looking at the block of objects between the datafiles and the tasks one might notice some objects not mentioned in the minimal example provided in the general example. Some of these might be specific to the project and your environment. The common ones that we use in our primary use case are:
```json
    "dnaAlignmentStyle": "",
    "email": "",
    "ethnicity": "",
    "familyCode": "",
    "holdConfig": false,
    "hpcAccount": "",
    "isilonPath": "",
    "matchedNormal": true,
    "matchedNormalToUse": "",
    "maternalID": "",
    "patCode": "",
    "paternalID": "",
    "pipeline": "",
    "project": "",
    "sex": "",
    "study": "",
    "submissionSource": "",
    "submitter": "",
```
  </p>
</details>


Once we have a config file for the project we're ready to initialize and launch the project. We can initialize a project via
```
$ jetstream init -h
usage: jetstream init [-h] [-l] [-p PROJECT] [-f] [--project-id PROJECT_ID]
                      [-c TYPE:KEY VALUE] [-C PATH]
                      [path]

Create or reinitialize a project This command is used to create a new
Jetstream project directory. If no path is given, the current directory will
be initialized. If config data options are given (-c/--config/--config-file),
they will be added to the project config file.

positional arguments:
  path                  Path to a initialize a project

optional arguments:
  -h, --help            show this help message and exit
  -l , --logging        set the logging profile
  -p PROJECT, --project PROJECT
                        path to a Jetstream project directory
  -f, --force           Force overwrite of project.yaml
  --project-id PROJECT_ID
                        Force a project ID instead of using letting it be
                        generated automatically

template variables:
  These options are used to add data that is available for rendering
  templates. These arguments should follow the syntax "-c <[type:]key>
  <value>". They can be used multiple times.

  -c TYPE:KEY VALUE, --config TYPE:KEY VALUE
                        add a single template variable
  -C PATH, --config-file PATH
                        load template variables from a file
                        
$ jetstream init GIAB -C GIAB_NA12878_24582bb3f7.json
```

This creates a jetstream project with the title of GIAB. Now in order to run the Phoenix pipeline on this project, we need to use:
```
$ jetstream pipelines -h
usage: jetstream pipelines [-h] [-l] [-p PROJECT] [-o OUT] [-b] [-r]
                           [--backend {local,slurm}]
                           [--format {template,module,workflow}]
                           [--reset-method {retry,resume,reset}]
                           [--existing-workflow EXISTING_WORKFLOW]
                           [--template-dir [SEARCH_PATH]] [-c TYPE:KEY VALUE]
                           [-C PATH] [--pipelines-home PIPELINES_HOME] [-L]
                           [path]

Run a pipeline. Pipelines are Jetstream templates that have been documented
with version information and added to the jetstream pipelines directory. This
command allows pipelines to be referenced by name and automatically includes
the pipeline scripts and constants in the run. Run Jetstream from a template,
module, or workflow

positional arguments:
  path                  path to a template, module, or workflow file. (if
                        using "pipelines" command, the name of the pipeline)

optional arguments:
  -h, --help            show this help message and exit
  -l , --logging        set the logging profile
  -p PROJECT, --project PROJECT
                        path to a Jetstream project directory
  -o OUT, --out OUT     path to save the workflow progress (this will be set
                        automatically if working with a project) [None]
  -b, --build-only      just render the template, build the workflow, and stop
  -r, --render-only     just render the template and stop
  --backend {local,slurm}
                        runner backend name used for executing tasks [slurm]
  --format {template,module,workflow}
                        workflow format - if this is None, it will be inferred
                        from the extension of the path [None]
  --reset-method {retry,resume,reset}
                        controls which tasks are reset prior to starting the
                        run - "retry": pending and failed, "resume": pending,
                        or "reset": all [retry]
  --existing-workflow EXISTING_WORKFLOW
                        path to an existing workflow file that will be merged
                        into run (this will be set automatically if working
                        with a project)
  --template-dir [SEARCH_PATH]
                        directory to add to search path for loading templates,
                        this can be used multiple times

template variables:
  These options are used to add data that is available for rendering
  templates. These arguments should follow the syntax "-c <[type:]key>
  <value>". They can be used multiple times.

  -c TYPE:KEY VALUE, --config TYPE:KEY VALUE
                        add a single template variable
  -C PATH, --config-file PATH
                        load template variables from a file

pipeline options:
  --pipelines-home PIPELINES_HOME
                        override path to the pipelines home
  -L, --list            show a list of all the pipelines installed

$ jetstream pipelines phoenix -p GIAB
```  
Now we wait for the pipeline to finish!


___


## Required Configuration Variables

For each of our data files/fastqs we have some required data, many of which are self explained,
but we will explain the more unique variables. Here is an example:

```json
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
    }
```

## Data file attributes

There are restrictions on what some of these variables can be assigned to, these will be denoted in the [ ]'s. 
If the attribute isn't strictly required then it is not included in this list.

  - *assayCode* [*WGS|*WGL|]  
    Used for determining if the sample is DNA/RNA/etc. and adding the corresponding
    tasks to the final workflow. Each sample discovered will take this attribute from
    the first file encountered for that sample in the config file.

  - *dnaRnaMergeKey*  
    Used during STAR alignment of RNA.
    
  - *fastqCode* [R1|R2]  
    Assigns the read number of the fastq.
    
  - *fastqPath*   
    Assigns the path to the fastq.

  - *fileType*  
    Assigns the file type.

  - *glPrep* [genome|capture|rna|singlecellrna|singlecellenrichment|singlecellcdna|singlecelltargetamp|matepair|chip]  
    Used for determining the prep used to create the sample and then modifying how the
    pipeline runs depending on the prep. This is used for single cell as well as CHIP preps.

  - *glType* [genome|genomephased|exome|rna|singlecellrna|singlecelldna|matepair|chip]  
    Used for determining if the sample is DNA/RNA/etc. and adding the corresponding
    tasks to the final workflow. Each sample discovered will take this attribute from
    the first file encountered for that sample in the config file.

  - *limsLibraryRecordId*  
    Generated by our LIMS, this allows for the input of data back into the LIMS.

  - *numberOfReads*  
    Used for assigning the number of chunks during alignment.
    
  - *read1Length / read2Length*   
    Used to assign the STAR indexes.
    
  - *readOrientation* [inward|outward]  
    Used to set the strandedness of RNA. Used in conjunction with rnaStrandDirection and rnaStrandType.

  - *rg values*  
    These are standards set in the [SAM/BAM Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf):  
    rgcn - Name of sequencing center producing the read  
    rgid - Read group identifier.  
    rgks - The array of nucleotide bases that correspond to the key sequence of each read.  
    rglb - Library.  
    rgpl - Platform/technology used to produce the reads.  
    rgpm - Platform model. Free-form text providing further details of the platform/technology used.  
    rgpu - Platform unit (e.g., flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.  
    rgsm - Sample. Use pool name where a pool is being sequenced.  
    
  - *fraction*  
    Relevant to the TGen naming scheme. See TGen Naming Convention.

  - *rnaStrandDirection* [notapplicable|forward|reverse]  
    Used during STAR alignment of RNA.
    
  - *rnaStrandType* [unstranded|stranded]  
    Assigns the strandedness of RNA
    
  - *sampleMergeKey*   
    This is the expected BAM filename and is used to merge data from multiple sequencing 
    lanes or flowcells for data from the same specimen tested with the same assay

  - *sampleName*  
    This is the expected base FASTQ filename.

  - *subGroup*   
    Sets where the data file is for tumour or constitutional, changes the analysis of the data file as well as sets
    the distinction of files during somatic analysis.

## TGen Naming Convention
- The naming format used in the LIMS tries to ensure that all specimens have unique file names based on specific metadata fields

STUDY_PATIENT_VISIT_SOURCE_FRACTION_SubgroupIncrement_ASSAY_LIBRARY<br/>
-------- PatientID<br/>
---------------- VisitID<br/>
------------------------ SpecimenID<br/>
-------------------------------- SampleID<br/>
---------------------------------------- RG.SM (VCF file genotype column header)<br/>
------------------------------------------------ sampleMergeKey (BAM filename)<br/>
