# gwasrecode is a program for extracting and recoding specific markers from PLINK associated .tped-files

**It was designed for non programmers**

**However, one need minimal command line knowledge to use to program properly (see README-FILE)**

## Functionality and description
- extract specific markers from .tped-files based on:
  - a list of markers saved in a file (each SNP on a separate line)
  - a reference SNP and a definition of bases downstream of the reference SNP - all markers in this range will be extracted
- the .tped files can have common rs-identifier or Illumina exm-identifier (exome data)
  - for the exome-data scenario, a Illumina auxillery file is necessary which contains the annotation from exm- to rs-identifiers
    - [Download page for Illumina support files e.g. auxillery files for specific arrays](http://support.illumina.com/downloads.html)
    - the auxillery file should has two columns: exm-identifier (column 1) and rs-identifer (column 2)
	- the user has to bring the file into the above-mentioned format (if it is not the existing format) 
  - the user only has to supply common rs-identifiers, irrespective whether it is exome- or common SNP-chip data

  - gwasrecode was tested under Linux and Windows (version 7 or higher)

- **A DETAILED DESCRIPTION OF THE USAGE, ALSO FOR THE GENERAL COMMAND LINE USE, IS WRITTEN IN THE README-FILE**

## Requirements
Python 2.7.x

## Installing
Just download the .zip-file and extract it in a specific folder.

The ReadMe has its own container since there are two pictures which refer to the python installation.

## Bugs
Because the program was needed it is not tested for all possible bugs or errors. 
If you find errors which make no sense or you have a problem in general, please contact me.

## Notes
In the near future I want to add extensive comments for the different functions and also to clean the code.

I have to admit that the program is not coded in best coding practice in version 0.34.

## Example usage (this is also written in the ReadMe)

###### Example 1 (based on Windows cmd)
You have exome data and want to extract all SNPs 40000 base pairs down stream of a reference SNP.

Beside -t, -o and -b  you have to supply the options -a -s -d -c


- Assume the "test.tped" and "test.tfam" file are stored in the folder "D:\Test\GWAS"
- Assume the auxillery file "data_aux.txt" is stored in the folder "D:\Test\Exomeinfo"
- Assume the reference snp rs1234 on chromosome 2
- Assume we want to extract all SNPs 40000 base pairs down stream of the reference snp
- Assume we only want the recoded models for all extracted SNPs and not the single output file per SNP

**The command would be:**
```
python gwasrecode.py -t d:\Test\GWAS\test.tped -o outtest.txt -b no -a d:\Test\Exomeinfo\data_aux.txt -s rs1234 -d 40000 -c 2
```


###### Example 2 (based on Windows cmd)
You have common SNP-chip data and want to extract a list of SNPs.

Beside -t, -o and -b you only have to supply the option -l.


- Assume the "test.tped" and "test.tfam" file are stored in the folder "D:\Test\GWAS"
- Assume the SNP-list-file "snps.txt", stored in the folder "D:\Test\GWAS\List" (this file contains rs-numbers you want to extract, every identifier on a new line)
- Assume we want both, the recoded models for all extracted SNPs and the single output file per SNP

**The command would be:**
```
python gwasrecode.py -t d:\Test\GWAS\test.tped -o outtest2.txt -b yes -l d:\Test\GWAS\List\snps.txt
```