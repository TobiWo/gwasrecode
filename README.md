# gwasrecode is a program for extracting and recoding specific markers from PLINK associated .tped-files

## Functionality and description
- extract specific markers from .tped-files based on:
  - a list of markers saved in a file (each SNP on a separate line)
  - a reference SNP and a definition of bases downstream of the reference SNP - all markers in this range will be extracted
- the .tped files can have common rs-identifier or Illumina exm-identifier
  - for the exome-data scenario, a Illumina auxillery file is necessary which contains the annotation from exm- to rs-identifiers
    - the auxillery file should has two columns: exm-identifier (column 1) and rs-identifer (column 2)
	- the user has to bring the file into the above-mentioned format (if it is not the existing format) 
  - the user has to supply common rs-identifiers irrespective whether it is exome- or common SNP-chip data
- OS independent
