import os
import platform
import re


class AllelRecode(object):

    def __init__(self, tppath, tppath_orig):
        self.tfam_path = tppath_orig
        self.tfam_path = self.tfam_path.replace(".tped", ".tfam")
        self.tped_path = tppath

    def rec_output(self):

        if platform.system() == "Linux":
            path_index = self.tped_path.rindex("/")
        elif platform.system() == "Windows":
            path_index = self.tped_path.rindex("\\")
        working_path = os.path.abspath(self.tped_path[0:(path_index + 1)])
        if not os.path.exists(os.path.join(working_path, "recoded")):
            working_path = os.path.join(working_path, "recoded")
            os.makedirs(working_path)
            return working_path
        else:
            working_path = os.path.join(working_path, "recoded")
            return working_path

    def allel_counter(self, x):
        """

            :rtype: counts allels and returns a list of two tuples
                    tuples consist of allel letter and the number present in the provided population
        """

        allel_list = list()

        allel_a = x.count("A")
        allel_c = x.count("C")
        allel_t = x.count("T")
        allel_g = x.count("G")
        allel_zero = x.count("0")

        if allel_a != 0:
            allel_list.append(("A", allel_a))
        if allel_c != 0:
            allel_list.append(("C", allel_c))
        if allel_t != 0:
            allel_list.append(("T", allel_t))
        if allel_g != 0:
            allel_list.append(("G", allel_g))
        if allel_zero != 0:
            allel_list.append(("0", allel_zero))

        return allel_list

    def tfam_ls(self):
        """

            :rtype: just read in a tfam-file to extract the ids and store it in a list for later use
        """
        tfam = open(self.tfam_path, "r")

        tfam_list = list()

        for line in tfam:
            line = line.rstrip()
            space = line.find(" ")
            ids = line[0:space]
            tfam_list.append(ids)

        return tfam_list

    def recoder(self):

        tped = open(self.tped_path, "r")

        for line in tped:

            id_counter = 0
            recode_list = ["id genotype recode add dom rez"]

            line = line.rstrip()
            if line == '':
                continue
            rs = re.search("rs[0-9]+", line).group(0)
            re_c = re.compile("[0-9]+ [A-Z]")
            pos = re_c.search(line).end()
            line2 = line[(pos-1):]
            line2 = line2.replace(" ", "")

            allels = self.allel_counter(line2)

            if len(allels) < 2:
                continue
            else:
                if allels[0][1] > allels[1][1]:
                    ma_allel = allels[0][0]
                    mi_allel = allels[1][0]
                    ma_allel_num = float(allels[0][1])
                    mi_allel_num = float(allels[1][1])
                else:
                    ma_allel = allels[1][0]
                    mi_allel = allels[0][0]
                    ma_allel_num = float(allels[1][1])
                    mi_allel_num = float(allels[0][1])

            recode_list.insert(0, "Cohort.based.Frequencies::::MajorAllel:%s_MinorAllel:%s" % (round(((ma_allel_num)/(ma_allel_num+mi_allel_num)), 4),
                                                                                                   round(((mi_allel_num)/(ma_allel_num+mi_allel_num)), 4)))
            recode_list.insert(0, "MajorAllel:%s____MinorAllel:%s" % (ma_allel, mi_allel))

            out = os.path.join(self.rec_output(), rs+"_recode"+".txt")
            out = open(out, "w")

            tfam_list = self.tfam_ls()

            for i in range(0, len(line2), 2):
                allels = line2[i:(i+2)]
                if allels[0] == ma_allel and allels[1] == ma_allel:
                    allels = tfam_list[id_counter] + " " + allels + " 11 0 1 1"
                    id_counter += 1
                elif allels[0] == mi_allel and allels[1] == mi_allel:
                    allels = tfam_list[id_counter] + " " + allels + " 22 2 2 2"
                    id_counter += 1
                elif allels[0] == str(0) and allels[1] == str(0):
                    allels = tfam_list[id_counter] + " NA NA NA NA NA"
                    id_counter += 1
                else:
                    allels = tfam_list[id_counter] + " " + allels + " 12 1 1 2"
                    id_counter += 1

                recode_list.append(allels)

            for item in recode_list:
                item.rstrip()
                out.write(item+"\n")

            out.close()

        return None

    def one_model_recode(self):

        tped = open(self.tped_path, "r")

        out_snps = ["SNP Homozygous_allele"]
        out2 = os.path.join(self.rec_output(), "homozygous_snps" + ".txt")
        out2 = open(out2, "w")

        for m in ["add", "dom", "rez"]:
            model = m
            recode_list = list()

            id_list = self.tfam_ls()
            id_list.insert(0, "ID")

            recode_list.append(id_list)

            for line in tped:

                line = line.rstrip()
                if line == '':
                    continue
                rs = re.search("rs[0-9]+", line).group(0)
                re_c = re.compile("[0-9]+ [A-Z]")
                pos = re_c.search(line).end()
                line2 = line[(pos - 1):]
                line2 = line2.replace(" ", "")

                allels = self.allel_counter(line2)

                if len(allels) < 2:
                    if model == "add":
                        out_snps_line = rs + " " + allels[0][0]
                        out_snps.append(out_snps_line)
                    continue
                else:
                    if allels[0][1] > allels[1][1]:
                        ma_allel = allels[0][0]
                        mi_allel = allels[1][0]
                    else:
                        ma_allel = allels[1][0]
                        mi_allel = allels[0][0]

                out = os.path.join(self.rec_output(), model + "_all_snps.txt")
                out = open(out, "w")

                collection = list()

                for i in range(0, len(line2), 2):
                    allels = line2[i:(i + 2)]
                    if allels[0] == ma_allel and allels[1] == ma_allel:
                        allels = "0 1 1"
                    elif allels[0] == mi_allel and allels[1] == mi_allel:
                        allels = "2 2 2"
                    elif allels[0] == str(0) and allels[1] == str(0):
                        allels = "NANANA"
                    else:
                        allels = "1 1 2"

                    if model == "add":
                        allels = allels[0]
                    elif model == "dom":
                        allels = allels[2]
                    elif model == "rez":
                        allels = allels[4]

                    collection.append(allels)

                collection.insert(0, rs)
                recode_list.append(collection)

                recode_list_t = [list(x) for x in zip(*recode_list)]

            if len(recode_list) == 1:
                raise IndexError("Only homozygous SNPs in data-set")

            for item in recode_list_t:
                item = " ".join(item)
                out.write(item + "\n")

            out.close()

            if model == "add":
                for item in out_snps:
                    item.rstrip()
                    out2.write(item + "\n")

                out2.close()

            tped.seek(0)

        return None


#liste = AllelRecode(tppath="/media/tobias/Daten/Daten/Projekte/GWAS_with_R/GWAS_Daten_Sorbs_raw/Daten_MPI/Entpackt/tansposed/extracted.txt", tppath_orig="/media/tobias/Daten/Daten/Projekte/GWAS_with_R/GWAS_Daten_Sorbs_raw/Daten_MPI/Entpackt/tansposed/test.tped")
#
#liste.one_model_recode()
