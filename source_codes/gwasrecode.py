from optparse import OptionParser
from Extract_SNPS_Tped_viaDistance_v035 import ExtractSNP
import allel_recoding_v035 as ac
import os
import platform
import sys

################################################################################################################
# OPTPARSER Arguments

usage = """\nProgram to extract and recode SNPs from tped-file!
The program can handle four scenarios!\nHereinafter the scenarios and necessary options are listed (Note: The options -t, -b and -o are always mandatory):

Exome-SNPs via Distance: -a -s -d -c
Exome-SNPs via list of SNPs: -a -l
Common SNPs via Distance: -s -d -c
Common SNPs via list of SNPs: -l

NOTE: THE OPTIONS LISTED ABOVE ARE MANDATORY FOR THE RESPECTIVE SCENARIOS.
      IF YOU DO NOT SUPPLY THEM TO THE PROGRAM IT WILL CRASH!!!"""

parser = OptionParser(usage=usage)

parser.add_option("-t", "--tped", dest="file", action="store", type="string", help="Complete path to tped file (Note: Corresponding tfam-file should be in the same directory)")
parser.add_option("-o", "--out", action="store", dest="output", type="string", help="Name of output file which will contain the raw_data of the desired SNPS!")
parser.add_option("-b", "--single", action="store", dest="single_file", type="string", help="Type 'yes' or 'no' to define, whether information on single SNPs should be created")
parser.add_option("-a", "--aux", action="store", dest="auxillery", type="string", help="Optional: Path to illumina auxillery file")
parser.add_option("-s", "--snp", action="store", dest="snpname", type="string", help="Optional: rs-number of reference SNP")
parser.add_option("-d", "--dist", action="store", dest="distance", type="int", help="Optional: Distance in base pairs which is used to extract all SNPs downstream of the reference SNP")
parser.add_option("-c", "--chr", action="store", dest="chromosome", type="int", help="Optional: Chromosome")
parser.add_option("-l", "--list", action="store", dest="snplist", type="string", help="Optional: Path to file which contains SNPs (every SNP on a new line)")

(options, args) = parser.parse_args()

################################################################################################################


ExtractObject = ExtractSNP(out=options.output, tppath=options.file, aux=options.auxillery, snpls=options.snplist, refsnp=options.snpname,
                           chr=options.chromosome, dist=options.distance)

if options.snpname is not None and options.snplist is None:
    if options.auxillery is not None:
        ExtractObject.exome_refsnp_extraction()
    elif options.auxillery is None:
        ExtractObject.commonsnp_refsnp_extraction()
elif options.snpname is None and options.snplist is not None:
    if options.auxillery is not None:
        ExtractObject.exome_snplist_extraction()
    elif options.auxillery is None:
        ExtractObject.commonsnp_snplist_extraction()

rec_out = options.file

if platform.system() == "Linux" or platform.system() == "Darwin":
    rec_index = rec_out.rindex("/")
elif platform.system() == "Windows":
    rec_index = rec_out.rindex("\\")
else:
    rec_index = None
    sys.exit("Operating system %s not supported" % platform.system())

rec_out = rec_out[0:rec_index]
rec_out = os.path.join(rec_out, options.output)
final_list = ac.AllelRecode(tppath=rec_out, tppath_orig=options.file)

if options.single_file.lower() == "yes":
    final_list.recoder()
    final_list.one_model_recode()
elif options.single_file.lower() == "no":
    final_list.one_model_recode()
elif options.single_file is None:
    sys.exit("-b flag is mandatory - Please enter 'yes' or 'no'")
