# coding=utf-8

import  sys
input_file = sys.argv[1]
output_file = sys.argv[2]


# 将parent的值最为key，每个RNA创建一个字典的键值对，正负链信息用0和1表示，0代表负链，1代表正链
# 字典结构，key为每个RNA的编号，值为列表，列表包含两个列表，第一个列表包含以下信息，正负链值为0和1，染色体和注释信息
# 第二个列表为exon的起止位置

all_rna = {}
with open(input_file,"r") as f:
    for line in f:
        Chr = line.split()[0]
        parent_rna = line.split(";")[1].split("=")[1]
        start_pos = int(line.split()[1])
        end_pos = int(line.split()[2])
        strand = line.split()[5]
        if strand == "+":
            strand = 1
        else:
            strand = 0
        id = line.split(";",1)[0].split()[6].split("=")[1]
        annotation = line.split(";",1)[1]
        # record = [start_pos,end_pos]
        if parent_rna in all_rna:
            # all_rna[parent_rna] = [all_rna[parent_rna][0],all_rna[parent_rna][1].append([start_pos,end_pos])]
            all_rna[parent_rna][1].append(start_pos)
            all_rna[parent_rna][1].append(end_pos)
        else:
            all_rna[parent_rna] = [[strand,Chr,annotation],[start_pos,end_pos]]
            # all_rna[parent_rna] = [strand,start_pos,end_pos]
# print(all_rna)

# 判断字典里面每一个rna记录值的正负链
with open(output_file,"w") as f:
    for value in all_rna.values():
        Chr = value[0][1]
        annotation = value[0][2]
        sort_list = value[1]
        sort_list.sort()
        #判断有几个exon
        exon_num = len(sort_list) / 2
        if value[0][0]:
            strand = "+"
            #此时序列是正链
            # print(sort_list)
            pos_num = -1
            for i in range(int(exon_num)-1):
                intron_num = i + 1
                pos_num += 2
                intron_start_pos = sort_list[pos_num] + 1
                intron_end_pos = sort_list[pos_num + 1] - 1
                Line = "%s %s %s %s %s %s" % (Chr, intron_start_pos, intron_end_pos, strand, intron_num, annotation)
                f.write(Line)
                # print(Line)
        else:
                #此时是负链
            strand = "-"
            # print(sort_list)
            pos_num = -1
            for i in range(int(exon_num) - 1):
                intron_num = i + 1
                pos_num -= 2
                intron_start_pos = sort_list[pos_num] + 1
                intron_end_pos = sort_list[pos_num + 1] - 1
                Line = "%s %s %s %s %s %s" % (Chr, intron_start_pos, intron_end_pos, strand, intron_num, annotation)
                f.write(Line)
