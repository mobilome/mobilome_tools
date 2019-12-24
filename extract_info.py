#coding=utf-8
import  sys
input_file_name = sys.argv[1]
output_file_name = sys.argv[2]
# input_file_name = "3-SINERTIP-merged-all.bed"
# output_file_name = "3-SINERTIP-merged-all.bed.out"

REF_column = 3 #
Breed_column = 9 #

flag = "-"

BE = flag
BM = flag
BS = flag
CB = flag
EG = flag
GT = flag
HS = flag
JH = flag
LD = flag
LW = flag
MS = flag
PT = flag
RC = flag
TB = flag
WZ = flag


with open(input_file_name,"r") as f:
    input_file = f.readlines()
dict = {}
for i in input_file:
    REF_Number = i.split("\t")[REF_column]
    Breed_Number = i.split("\t")[Breed_column]
    if REF_Number in dict:
        dict[REF_Number] = dict[REF_Number] + "\t" + Breed_Number
    else:
        dict[REF_Number] = Breed_Number



with open(output_file_name,"w") as f:
    f.write("REF" + "\t" + "BE"+"\t"+ "BM"+"\t"+"BS"+"\t"+"CB"+"\t"+"EG"+"\t"+"GT"+"\t"+"HS"+"\t"+"JH"+"\t"+"LD"+"\t"+"LW"+"\t"+"MS"+"\t"+"PT"+"\t"+"RC"+"\t"+"TB"+"\t"+"WZ"+"\n")
    for k, v in dict.items():
        for i in v.split():
            # s = i.split("-")[0] + '= "+"'
            s = i.split("-")[0]
            # exec(s)
            globals()[s] = "+"


        f.write(k+"\t"+ BE+"\t"+BM+"\t"+BS+"\t"+CB+"\t"+EG+"\t"+GT+"\t"+HS+"\t"+JH+"\t"+LD+"\t"+LW+"\t"+MS+"\t"+PT+"\t"+RC+"\t"+TB+"\t"+WZ+"\n")
        BE = flag
        BM = flag
        BS = flag
        CB = flag
        EG = flag
        GT = flag
        HS = flag
        JH = flag
        LD = flag
        LW = flag
        MS = flag
        PT = flag
        RC = flag
        TB = flag
        WZ = flag

