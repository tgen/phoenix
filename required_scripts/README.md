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
5f25b11
$ cp plotCNVplus.R /path/to/phoenix/required_scripts/plotCNVplus_5f25b11.R
```