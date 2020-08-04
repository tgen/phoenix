# Required Scripts
This directory contains scripts for at least one of the following reasons
 - From a private source with no intention of being made public 
 - Scripts where making a module does not make sense.

### Example 1
runIchorCNA
 - External R library distributed by CRAN but is executed from a provided RunIchorCNA.R file, that is not executable/detect within R session
### Example 2
plotCNVplus
 - This is a stand alone script with two R library dependencies (optparse, data.table). Making a module for this script that would put it in the users path would be wasteful of institutional resources.  
## Development Guidelines
 - Copy file into the external_scripts directory with the git commit ID embedded into the filename.
 - Clearly update the Documentation section of this README.md with the following information about the script.
    - How was the script acquired (git clone path, wget path)
    - What the Commit id is (git rev-parse --short HEAD)
    
## Documentation
### plotCNVplus
```
$ git clone git@github.com:tgen/plotCNVplus.git
$ cd plotCNVplus
$ git rev-parse --short HEAD
06b34ff
$ cp plotCNVplus.R /path/to/phoenix/required_scripts/plotCNVplus_06b34ff.R
```

## Script Source Locations
[addRC_to_Delly_VCF_f4d178e.py](https://github.com/tgen/jetstream_resources/commit/f4d178e2b8982ff49025d42cb7c18d7b12053f42)  
[annotSeg_8820499.pl](https://github.com/tgen/jetstream_resources/commit/8820499e113a387fee98044112951fa534ad6f8e)  
[manta_prepare_sv_vcf_ba6879d.py](https://github.com/tgen/jetstream_resources/commit/ba6879d05fad13b8e6c7aecd86498bc5c8e1b76e)  
[manta_sv_annotation_parallel_8820499.py](https://github.com/tgen/jetstream_resources/commit/8820499e113a387fee98044112951fa534ad6f8e)  
[mm_igtx_pairoscope_calling_b38_f3010c3.py](https://github.com/tgen/mm_IgTx_Calling/commit/f3010c358970f4c25cefddc824636f60a19842e1)  
[plotCNVplus_4d89cb4.R](https://github.com/tgen/plotCNVplus/commit/4d89cb4d8f35e48b916d660c82c52b8725ade16f)  
[runIchorCNA_47ce8db.R](https://github.com/broadinstitute/ichorCNA/commit/47ce8db4d81ada2d3ce09280661d1240f3dcd530#diff-79cb887cc56cef135b77c5b7a725975c)  
[samStats2json_deb00a0.py](https://github.com/tgen/samStats2json/commit/deb00a0199dca764da0542584616ad4333b9e8b5)  
[seg_extend_229b8c7.py](https://github.com/tgen/jetstream_resources/commit/229b8c7641dd505789664aab88c1662d1f97e429)  
[summarize_samstats_69627ba.R](https://github.com/tgen/plot_samstats/commit/69627baab42126a991097700c0ae0ef91f22b586)  
[svtop.delly.sv_annotation.parallel_8820499.py](https://github.com/tgen/jetstream_resources/commit/8820499e113a387fee98044112951fa534ad6f8e)  
[uploadStats2Lims_98319ea.py](https://github.com/tgen/uploadStats2Lims/commit/98319eaa56cad71c726df308e68262156c152caf)  