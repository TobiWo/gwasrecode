# gwasrecode is a program for extracting and recoding specific markers from PLINK associated .tped-files

**It was designed for non programmers

**However, one need minimal command line knowledge to use to program properly (see README-FILE)

## Functionality and description
- extract specific markers from .tped-files based on:
  - a list of markers saved in a file (each SNP on a separate line)
  - a reference SNP and a definition of bases downstream of the reference SNP - all markers in this range will be extracted
- the .tped files can have common rs-identifier or Illumina exm-identifier (exome data)
  - for the exome-data scenario, a Illumina auxillery file is necessary which contains the annotation from exm- to rs-identifiers
    - the auxillery file should has two columns: exm-identifier (column 1) and rs-identifer (column 2)
	- the user has to bring the file into the above-mentioned format (if it is not the existing format) 
  - the user only has to supply common rs-identifiers, irrespective whether it is exome- or common SNP-chip data
- gwasrecode was tested under Linux and Windows (version 7 or higher)
- **A DETAILED DESCRIPTION OF THE USAGE, ALSO FOR THE GENERAL COMMAND LINE USE, ARE WRITTEN IN THE README-FILE

## Requirements
Python 2.7.x

## Installing
Just download the .zip-file and extract it in a specific folder

The ReadMe has its own container since there are two pictures which refer to the python installation

## Bugs
Because the program was needed it is not tested for all possible bugs or errors. 
If you find errors which make no sense or you have a problem in general, please contact me.

## Notes
In the near future I want to add extensive comments for the different functions and also to clean the code.

I have to admit that the program is not coded in best coding practice in version 0.34.