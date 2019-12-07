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

All tools build with public easybuild configuration files - **REPO-X**.<br/>
Last Update November 14, 2019  


| Tool | Version Implemented | Current Version | Dependancy and Notes | EasyBuild |
| :---: | :---: | :---: | :--- | :---: |
| [bcftools](https://github.com/samtools/bcftools/releases) | 1.9 | 1.9 | | Yes |
| [bedtools](https://github.com/arq5x/bedtools2/releases) | 2.29.0 | 2.29.0 | delly-filter, addmatchRNA, vardict, vcfmerger2 | Yes |
| [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) | 2.3.5.1 | 2.3.5.1 | star-fusion | Yes |
| [bwa](https://github.com/lh3/bwa/releases) | 0.7.17 | 0.7.17 | | Yes |
| [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) | 3.1.0 | 3.1.0 | | Yes |
| [deepvariant](https://github.com/google/deepvariant/releases) | 0.8.0 | 0.8.0 | singularity container | Yes |
| [delly](https://github.com/dellytools/delly/releases) | 0.7.6 | **0.8.1** | staying with 0.7.6 for compatibility reasons | Yes |
| [freebayes](https://github.com/ekg/freebayes/releases) | 1.3.1 | 1.3.1 | update allows skipping of high coverage regions | Yes |
| [gatk](https://github.com/broadinstitute/gatk//releases) | 4.1.4 | 4.1.4 |  | Yes |
| [gmap-gsnp](http://research-pub.gene.com/gmap/) | 2019-09-12 | 2019-09-12 | star-fusion | Yes |
| [gridss](https://github.com/PapenfussLab/gridss/releases) | 2.6.3 | 2.7.2 |  | |
| [hmmcopyutils](https://github.com/shahcompbio/hmmcopy_utils) | 1.0 | 1.0 | no official release | Yes |
| [htseq](https://github.com/simon-anders/htseq/releases) | 0.11.1 | 0.11.1 | | |
| [htslib](https://github.com/samtools/htslib/releases) | 1.9 | 1.9 | star-fusion(bgzip) | Yes |
| [ichor](https://github.com/broadinstitute/ichorCNA/releases) | 0.2.0 | 0.2.0 | package in R/3.6.1-phoenix module | Yes? |
| [jellyfish](https://github.com/gmarcais/Jellyfish/releases) | 2.3.0 | 2.3.0 | star-fusion | Yes |
| [lancet](https://github.com/nygenome/lancet/releases) | 1.0.7 | 1.0.7 | | Yes |
| [manta](https://github.com/Illumina/manta/releases) | 1.6.0 | 1.6.0 | | Yes |
| [multiQC](https://github.com/Illumina/manta/releases) | 1.7 | 1.7 | python3 pip | Yes |
| [octopus](https://github.com/luntergroup/octopus/releases) | 0.6.3-beta | 0.6.3-beta | | Yes |
| [perl](https://github.com/Illumina/manta/releases) | 5.28.1 | 5.30.0 | star-fusion | Yes |
| [phaser](https://github.com/secastel/phaser/tree/master/phaser) | 1.1.1 | 1.1.1 | vcfmerger2 | Yes |
| [python2](https://www.python.org/downloads/) | 2.7.15 | 2.7.15 | | Yes |
| [python3](https://www.python.org/downloads/) | 3.7.2 | 3.8.0 | star-fusion, vcfmerger2 | Yes |
| [R](https://www.r-project.org/) | 3.6.1 | 3.6.1 | gatk cnv, varDict, vcfmerger2 | Yes |
| [sambamba](https://github.com/biod/sambamba/releases) | 0.7.0 | 0.7.0 | | |
| [samblaster](https://github.com/GregoryFaust/samblaster/releases) | 0.1.24 | 0.1.24 | | |
| [salmon](https://github.com/COMBINE-lab/salmon/releases) | 0.14.1 | 0.14.2 | self, star-fusion | Yes |
| [samtools](https://github.com/samtools/samtools/releases) | 1.9 & 1.9-168-gb1e2c78 | 1.9 | markdup needs unreleased patch | Yes |
| [singularity](https://github.com/COMBINE-lab/salmon/releases) | 3.5.0-rc2 | 3.4.2 | deepvariant | Yes, release-candidate |
| [snpEff](https://sourceforge.net/projects/snpeff/files/) | 4.3t | 4.3t | no updates since Nov 2017| Yes |
| [star](https://github.com/alexdobin/STAR/releases) | 2.7.3a | 2.7.3a | self, star-fusion | Yes |
| [star-fusion](https://github.com/STAR-Fusion/STAR-Fusion/releases) | 1.8.1 | 1.8.1 | | Yes |
| [strelka](https://github.com/Illumina/strelka/releases) | 2.9.10 | 2.9.10 | | Yes |
| [subread](https://sourceforge.net/projects/subread/) | 2.0.0 | 2.0.0 | part of subread package | Yes |
| [trinityrnaseq](https://github.com/trinityrnaseq/trinityrnaseq/releases) | 2.8.6 | 2.8.6 | star-fusion | Yes |
| [vardictJava](https://github.com/AstraZeneca-NGS/VarDictJava/releases) | 1.7.0 | 1.7.0 | | Yes |
| [vcfmerger2](https://github.com/tgen/vcfMerger2/releases) | 0.8.3 | 0.8.3 | | Yes |
| [vep](https://github.com/Ensembl/ensembl-vep/releases) | 98.2 | 98.2 | | Yes |
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

### Single Cell RNA Sequencing Analysis
Two options exist for alignment and gene expression count generation per cell.
* 10x Genomics Cell Ranger 3.1.1
* starSOLO 2.7.3a

**Cell Barcode Whitelist Notes:**
PATH = /packages/cellranger/3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes
Single Cell 5' Gene Expression, transcript on R1 only = 737K-august-2016.txt
Single Cell 5' Gene Expression, transcript on R2 only = 737K-august-2016.txt
Single Cell 5' Gene Expression, paired-end = 737K-august-2016.txt
Single Cell V(D)J, transcript on R2 only = 737K-august-2016.txt
Single Cell V(D)J = 737K-august-2016.txt
CHEMISTRY_SC3P_V1 = 737K-april-2014_rc.txt
CHEMISTRY_SC3P_V2 = 737K-august-2016.txt
CHEMISTRY_SC3P_V3 = 3M-february-2018.txt.gz

## JJK Notes

dataFiles : JSON/YAML grouping that includes all fastq metadata fields
  * In main.main one of the first steps it so define each fastq meta data block as "file"
    * so afterwards "file" references each fastq block, be it an R1 or R2 read
    * then a series of additional metadata fields are added
      * name - this is the same as sampleMergeKey unless not defined in which case it will be sampleName
      * basename - proper fastq name not full path (ie. file.R1.fastq)
      * gltype - lower(glType)
      * glprep - lower(glPrep)
      
Single cell workflows are passed dataFiles object (ie all fastq sets)
  * in the single_cell/main
    * we group the dataFiles by "name" and test if the first file object is single cell or not
      * this creates unneeded code to then separate the libraries from within the group_by names records
      * IMPROVED code would group by library, rglb, as this is the analyzed unit for single cell
      * then the list of files submitted would be the needed unit of library summarized data files

## Install Guide
First, jetstream must be installed, available [here](https://github.com/tgen/jetstream/tree/master). There is an install guide for jetstream available at the link provided, but the gist of the guide is to install using pip3 using a command similar to:  
```
$ pip3 install --upgrade --user git+https://github.com/tgen/jetstream.git@master
```  
Then the recommended install method for installing the phoenix-development pipeline is:
```
$ cd ~
$ mkdir jetstream_pipelines
$ cd jetstream_pipelines
$ git clone --single-branch --branch develop https://github.com/tgen/phoenix
$ mv phoenix phoenix-development
```  
We're getting close to being able to easily run the pipeline now, and from this point you might be able to hack your way to make everything run. But it is recommended that you use similar settings to the ones detailed here in order to get the best support possible.

By running the following command you should be able to see the settings that jetstream is currently using:
```
$ jetstream settings -v
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
$ jetstream settings -c -b "slurm" -P "/home/USERNAME/jetstream_pipelines/"
```  
Note that USERNAME is replaced by your username, e.g. "/home/jsmith/jetstream_pipelines/" for a user name John Smith within TGen. Also, the home location will differ by install and should be the location of our downloaded pipelines. If you aren't sure where they are downloaded, use
```
$ cd ~
$ find . -name pipeline.yaml
```
This may return more than one result. The one we are looking for should look similar to "./jetstream_pipelines/phoenix/pipeline.yaml". The ./ means that we have jetstream_pipelines in our current directory. Which means that the true path to our pipelines is /home/USER/jetstream_pipeline/. Note that USER is your username, within TGen this is generally your first initial and then last name. Like the John Smith example above.  

We're nearly done now. At the time of writing, jetstream/phoenix is not entirely environment agnostic. The phoenix pipeline currently looks for reference data within our /home/tgenref/ directory. If we have/want to use data not within /home/tgenref/ we simply need to modify the pipeline.yaml for phoenix. 
We can view the pipeline.yaml by changing directories to where we downloaded the phoenix pipeline:  
```
$ cd ./jetstream_pipelines/phoenix-development/
$ less pipeline.yaml
```  
The areas that we are interested in are:
```
__pipeline__:
  name: phoenix
  main: main.jst
  description: Human GRCh38 genomics suite
  version: development
constants:
  install_path:
    path_to_phoenix_repo: /home/tgenjetstream/jetstream_pipelines/phoenix-development
.
.
.
phoenix:
    species: Homo sapiens
    genome_build: grch38_hg38
    reference_fasta: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa
    reference_fai: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.fa.fai
    reference_dict: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy_alts_hla.dict
    dbsnp_v152: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/dbsnp/b152/dbSNP_b152_hg38tgen.bcf
    gnomad_exome_v2_1_1_liftover: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38_NoINFO.bcf
    gnomad_genome_v3_0: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r3.0/gnomad.genomes.r3.0.sites.pass.AnnotationReference.bcf
    gnomad_exome_v2_1_1_mutect_germlinereference: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutect.vcf.gz
    gnomad_genome_v3_0_mutect_germlinereference: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r3.0/gnomad.genomes.r3.0.sites.pass.AnnotationReference.vcf.gz
    gnomad_exome_v2_1_1_mutect_contamination: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r2.1.1/gnomad.exomes.r2.1.1.sites.liftover_grch38_ForMutectContamination.vcf.gz
    gnomad_genome_v3_0_mutect_contamination: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/gnomad/r3.0/gnomad.genomes.r3.0.sites.pass.ForMutectContamination.vcf.gz
    cosmic_coding_v90: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/cosmic/v90/CosmicCodingMuts_v90_hg38tgen.bcf
    cosmic_noncoding_v90: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/cosmic/v90/CosmicNonCodingMuts_v90_hg38tgen.bcf
    clinvar_20190715: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/clinvar/20190715/clinvar_20190715_hg38tgen.bcf
    bcftools_annotate_contig_update_ucsc2ensembl: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/tool_resources/bcftools/GRCh38_PrimaryContigs_UCSC_2_Ensembl_CrossMap.txt
    black_list: /home/tgenref/homo_sapiens/grch38_hg38/public_databases/encode/Blacklist-2.0/lists/hg38-blacklist.v2.bed.gz
    delly_annotation: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/delly/delly_anno_Homo_sapiens.GRCh38.98.ucsc.bed
    delly_exclusions: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/delly/hg38.excl
    delly_addRC_to_Delly_VCF_script: addRC_to_Delly_VCF.py
    delly_svtop_delly_sv_annotation_parellel_script: svtop.delly.sv_annotation.parallel.py
    plotCNVplus_Rscript: plotCNVplus_06b34ff.R
    stats2json: samStats2json.py
    stats2lims: uploadStats2Lims.py
    cellranger_reference: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/cellranger_3.1.0/GRCh38_hg38tgen.98
    cellranger_vdj_reference: /home/tgenref/homo_sapiens/grch38_hg38/tool_specific_resources/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0
    scrna_chemistry_options:
      X3SCR:
        chemistry_name: SC3Pv1
        umi_length: 10
        cell_barcode_whitelist_file: /packages/cellranger/3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-april-2014_rc.txt
      XCSCR:
        chemistry_name: SC3Pv2
        umi_length: 10
        cell_barcode_whitelist_file: /packages/cellranger/3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt
      X3SC3:
        chemistry_name: SC3Pv3
        umi_length: 12
        cell_barcode_whitelist_file: /packages/cellranger/3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
      X5SCR:
        chemistry_name: SC5P-R2
        umi_length: 10
        cell_barcode_whitelist_file: /packages/cellranger/3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt
      unknown:
        chemistry_name: auto
        umi_lenth: 10
        cell_barcode_whitelist_file: /packages/cellranger/3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt
    gatk_known_sites:
      - /home/tgenref/homo_sapiens/grch38_hg38/public_databases/broad_resource_bundle/Homo_sapiens_assembly38.dbsnp138.vcf
      - /home/tgenref/homo_sapiens/grch38_hg38/public_databases/broad_resource_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz
      - /home/tgenref/homo_sapiens/grch38_hg38/public_databases/broad_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    gatk_cnn_resources:
      - /home/tgenref/homo_sapiens/grch38_hg38/public_databases/broad_resource_bundle/hapmap_3.3.hg38.vcf.gz
      - /home/tgenref/homo_sapiens/grch38_hg38/public_databases/broad_resource_bundle/Homo_sapiens_assembly38.known_indels.vcf.gz
      - /home/tgenref/homo_sapiens/grch38_hg38/public_databases/broad_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    bwa_index: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/tool_resources/bwa_0.7.17/GRCh38tgen_decoy_alts_hla.fa
    gtf: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/Homo_sapiens.GRCh38.98.ucsc.gtf
    ref_flat: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/Homo_sapiens.GRCh38.98.ucsc.refFlat.txt
    ribo_locations: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/Homo_sapiens.GRCh38.98.ucsc.ribo.interval_list
    gatk_cnv_primary_contigs_female: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/tool_resources/gatk_cnv/Homo_sapiens.GRCh38.primary.contigs.female.interval_list
    gatk_cnv_primary_contigs_male: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/tool_resources/gatk_cnv/Homo_sapiens.GRCh38.primary.contigs.male.interval_list
    transcriptome_fasta: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/Homo_sapiens.GRCh38.98.ucsc.transcriptome.fasta
    salmon_index: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/salmon_0.14.1/salmon_quasi_75merPlus
    sex_check_targets: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/tool_resources/tgen_gender_check/chrx_common_dbSNPv152_snv_exons.bed
    sex_check_vcf: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/tool_resources/tgen_gender_check/chrx_common_dbSNPv152_snv_exons.vcf.gz
    lymphocyteReceptor_loci_bed: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/tool_resources/tgen_lymphocyteReceptor_counts/lymphocyteReceptor_loci.bed
    star_fasta: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/genome_reference/GRCh38tgen_decoy.fa
    star_indices:
      # Multiple STAR references defined here in order to accommodate 
      # RNA data with different read lengths.
      75bpReads: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/75bpReads
      82bpReads: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/82bpReads
      83bpReads: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/83bpReads
      86bpReads: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/86bpReads
      87bpReads: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/87bpReads
      89bpReads: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/89bpReads
      100bpReads: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/100bpReads
      109bpReads: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/109bpReads
      150bpReads: /home/tgenref/homo_sapiens/grch38_hg38/hg38tgen/gene_model/ensembl_v98/tool_resources/star_2.7.3a/150bpReads
    starfusion_index: /home/tgenref/homo_sapiens/grch38_hg38/tool_specific_resources/STAR-fusion/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir
.
.
.
```
In order to change the location that phoenix looks for reference data, one can either manually modify each individual line, or as long as we have not left the phoenix directory, we can use:
```
sed -i 's|/home/tgenref|/home/newLocation|g' pipeline.yaml
```
To change all /home/tgenref text to /home/newLocation where /home/newLocation is the location of where our new references are. We can also use sed to replace more of the paths to reference data if needed simply by replicating the pattern above.  

Congratulations! That's it! We now have jetstream and the phoenix pipeline installed. In the next section we will discuss how we can actually run the pipeline from the command line.  

## Running from command line
Now that we have phoenix installed and ready to run, we need to setup a project to run. To do so, we need to create a config file for our project. The general format is as follows:
```json
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


