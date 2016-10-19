import os
import platform
import sys
import re


class ExtractSNP(object):

    def __init__(self, out, tppath, chr=None, refsnp=None, dist=None, aux=None, snpls=None, ):

        self.path_tped = tppath
        self.chrom = chr
        self.snpname = refsnp
        self.distance = dist
        self.auxillery = aux
        self.snplist = snpls
        self.out = out

        self.path_out = tppath

        if platform.system() == "Linux" or platform.system() == "Darwin":
            path_index = self.path_out.rindex("/")
        elif platform.system() == "Windows":
            path_index = self.path_out.rindex("\\")
        else:
            path_index = None
            sys.exit("Operating system %s not supported" % platform.system())

        self.path_out = os.path.join(self.path_out[0:(path_index + 1)], self.out)

    def out_exist(self):

        self.path_out = os.path.abspath(self.path_out)

        if os.path.exists(self.path_out):
            print "The file '%s' already exists! If you not continue, a new file with suffix 'new' will be created!" % self.out
            file_check = raw_input("Do you want to continue? If yes, the file will be overwritten (type: yes or no): ")
            file_check = file_check.upper()
            if file_check == "YES":
                outfile = open(self.path_out, "w+")
            elif file_check == "NO":
                print "File was created with suffix 'new'!"
                self.path_out = os.path.abspath(self.path_out + "_new.txt")
                outfile = open(self.path_out, "w+")
        else:
            outfile = open(self.path_out, "w+")

        return outfile

    def tped_handler(self):
        self.path_tped = os.path.abspath(self.path_tped)
        tpedfile = open(self.path_tped)
        return tpedfile

    def aux_handler(self):
        self.auxillery = os.path.abspath(self.auxillery)
        exome = open(self.auxillery)
        return exome

    def get_refsnp(self):

        ref_snp = None
        for line in self.aux_handler():
            line = line.rstrip()
            if self.snpname in line and len(re.search("rs\S+", line).group(0)) == len(self.snpname):
                fline = line
                first_space = fline.find("\t")
                fline = fline[0:first_space]
                ref_snp = fline

        if ref_snp is None:
            sys.exit(self.snpname + " not covered on chip!")

        return ref_snp

    def exome_refsnp_extraction(self):

        ref_snp = self.get_refsnp()
        tpedfile = self.tped_handler()
        outfile = self.out_exist()

        for line in tpedfile:
            line = line.rstrip()
            if self.chrom == int(line[0:2]):
                if ref_snp in line and len(re.search("exm\S+", line).group(0)) == len(ref_snp):
                    pos = re.search("\s[0-9]{2,}\s", line)
                    start_pos = pos.group(0)
                    start_pos = int(start_pos)
                    end_pos = start_pos + self.distance

        if "end_pos" not in locals():
            sys.exit("Couldn't find snp %s in tped" % self.snpname)

        tpedfile.seek(0)

        collect_exm = list()

        exm_list = list()
        for line in tpedfile:
            line = line.rstrip()
            pline = line[0:2]
            pline = int(pline)
            if pline == self.chrom:
                p = re.search("\s[0-9]{2,}\s", line)
                pos = p.group(0)
                pos = int(pos)
                if pos >= start_pos and pos <= end_pos:
                    collect_exm.append(line)
                    snp_m = re.search("exm\S+", line)
                    if snp_m != None:
                        exm_list.append(snp_m.group(0))

        # first one should extract the rs-numbers
        # loop over the elements of the exm-number list and then loop over the auxillery file to find the exm-number
        # Note: Again I use the length of the exm-number to get a perfect match in combination with the 'in' statement
        # Note2: Some exm-numbers have two or more rs-numbers --> only the first rs-number will be extracted
        # IMPORTANT#####IMPORTANT## If a rs-number does not make sense you have to check this!!!! ######IMPORTANT######IMPORTANT#############
        rs_collect = list()
        exome = self.aux_handler()
        for s in exm_list:
            for line in exome:
                line = line.rstrip()
                if s in line and len(re.search("exm\S+", line).group(0)) == len(s):
                    rs_num = re.search("rs[0-9]+", line).group(0)
                    rs_collect.append(rs_num)
            exome.seek(0)  # because the loop over the auxillery file should be repeated for every exm-number...
            # ...I need to jump back to the start of the file after every complete auxillery loop

        # create a new empty list for the final lines (with replaced rs-numbers)
        # and do the replacing
        # the exvar-variable is only a iteration variable to get access to the exm-number-list and the rs-number-list and use these for search and replace
        collect_final = list()
        exvar = 0
        for item in collect_exm:
            item = item.replace(exm_list[exvar], rs_collect[exvar])
            exvar += 1
            collect_final.append(item)

        for item in collect_final:
            item.rstrip()
            outfile.write(item + "\n")

        outfile.close()
        tpedfile.close()
        exome.close()

    def commonsnp_refsnp_extraction(self):

        ref_snp = self.snpname
        tpedfile = self.tped_handler()
        outfile = self.out_exist()

        for line in tpedfile:
            line = line.rstrip()
            if ref_snp in line and len(re.search("rs\S+", line).group(0)) == len(ref_snp):
                start_pos = re.search("\s[0-9]{2,}\s", line).group(0)
                start_pos = int(start_pos)
                end_pos = start_pos + self.distance
        tpedfile.seek(0)

        # main extraction of the desired snps
        rs_list = list()
        for line in tpedfile:
            line = line.rstrip()
            pline = line[0:2]
            pline = int(pline)
            if pline == self.chrom:
                pos = re.search("\s[0-9]{2,}\s", line).group(0)
                pos = int(pos)
                if pos >= start_pos and pos <= end_pos:
                    rs_list.append(line)

        for item in rs_list:
            item.rstrip()
            outfile.write(item + "\n")

        outfile.close()
        tpedfile.close()

    def exome_snplist_extraction(self):

        snplist = open(os.path.abspath(self.snplist), "r")
        exome = self.aux_handler()
        tpedfile = self.tped_handler()
        outfile = self.out_exist()

        exm_save = list()
        rs_save = list()
        for snp in snplist:
            snp = snp.rstrip()
            for line in exome:
                line = line.rstrip()
                if snp in line and len(re.search("rs\S+", line).group(0)) == len(snp):
                    exm = line[0:line.find("\t")]
                    exm_save.append(exm)
                    rs = line[(line.find("\t") + 1):]
                    rs_save.append(rs)
                    break
            exome.seek(0)

        tped_save = list()
        for snp in exm_save:
            for line in tpedfile:
                line = line.rstrip()
                if snp in line and len(re.search("exm\S+", line).group(0)) == len(snp):
                    tped_save.append(line)
                    break
            tpedfile.seek(0)

        collect_final = list()
        exvar = -1
        while len(collect_final) < len(exm_save) and exvar < len(exm_save):
            exvar += 1
            for item in tped_save:
                save_exm = re.search("exm\S+", item).group(0)
                if save_exm == exm_save[exvar]:
                    item = item.replace(exm_save[exvar], rs_save[exvar])
                    exvar += 1
                    collect_final.append(item)

        for item in collect_final:
            item.rstrip()
            outfile.write(item + "\n")

        outfile.close()
        tpedfile.close()
        exome.close()
        snplist.close()

    def commonsnp_snplist_extraction(self):

        tpedfile = self.tped_handler()
        outfile = self.out_exist()
        snplist = open(os.path.abspath(self.snplist), "r")

        collect_final = list()
        for snp in snplist:
            snp = snp.rstrip()
            for line in tpedfile:
                line = line.rstrip()
                if snp in line and len(re.search("rs\S+", line).group(0)) == len(snp):
                    collect_final.append(line)
                    break
            tpedfile.seek(0)

        for item in collect_final:
            item.rstrip()
            outfile.write(item + "\n")

        outfile.close()
        tpedfile.close()
        snplist.close()

# test = ExtractSNP(out = "test.txt", tppath="/media/tobias/HD-PNTU3/GWAS-Daten/Exome_Chips/Sorbs/test2.tped", chr=10, refsnp="rs241", dist=1000,
#                   aux="/media/tobias/HD-PNTU3/GWAS-Daten/Exome_Chips/Exome_Chip_details/HumanExome-12-v1-1-B-auxilliary-file.txt")
#
# test.exome_refsnp_extraction()