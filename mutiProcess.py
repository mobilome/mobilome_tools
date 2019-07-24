import os
from multiprocessing import Pool

def Compare(input):
    with open(input,"r") as f1:
        f1 = f1.readlines()
    with open("/seagate/Python/D.bed","r") as f2:
        f2 = f2.readlines()
    test_outfile = open(input+".out","w")
    for i in f1:
        Name1 = i.split("\t")[0]
        Start1 = i.split("\t")[1]
        End1 = i.split("\t")[1]
        for j in f2:
            if Name1 == j.split("\t")[0]:
                if abs(int(Start1)-int(j.split("\t")[1]))<100 or abs(int(End1)-int(j.split("\t")[2]))<100 :
                    i = i.strip('\n')
                    i = i.strip('\r')
                    test_outfile.writelines(i+"\t"+j)
    test_outfile.close()
pathss=[]
for root, dirs, files in os.walk("/seagate/Python/split"):
    path = [os.path.join(root, name) for name in files]
        #print(path)
    pathss.extend(path)



if __name__ ==  '__main__':
    p = Pool()
    p.map(Compare,pathss)
    p.close()
    p.join()
