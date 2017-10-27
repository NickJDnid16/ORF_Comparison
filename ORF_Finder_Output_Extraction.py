import collections
import fileinput
import linecache
import sys
import numpy
from Bio import SeqIO


TransDecoderInput = open('./longest_orfs.cds', mode='rb')
TransDecoderBest = open('./TransDecoder_Best_Candidates.txt', mode='rb')
#MetaGeneMarkOutput = open('./Sample_Data.txt.lst', mode='rb')
TransDecoderOutput = open('./TransDecoder_ORFs.csv', mode='wb')
MetaGeneMarkOutput = open('./MetaGeneMark_ORFs.csv', mode='wb')
Alignments = open('./Sample_Data.txt', mode='rb')

TransDecoderBestOutput = open('./TransDecoder_Best_ORFs.txt', mode='wb')


TransDecoderORFs = collections.OrderedDict()
MetaGeneMarkORFs = collections.OrderedDict()

TransDecoderLengths = []


def N50_(numlist):
    """
    Abstract: Returns the N50 value of the passed list of numbers.
    Usage:    N50(numlist)
    Based on the Broad Institute definition:
    https://www.broad.harvard.edu/crd/wiki/index.php/N50
    """
    numlist.sort()
    newlist = []
    for x in numlist :
        newlist += [x]*x
    # take the mean of the two middle elements if there are an even number
    # of elements.  otherwise, take the middle element
    if len(newlist) % 2 == 0:
        medianpos = len(newlist)/2
        return float(newlist[medianpos] + newlist[medianpos-1]) /2
    else:
        medianpos = len(newlist)/2
    return newlist[medianpos]


for line in TransDecoderInput:
    if ">Gene" in line :
        Len = line.split(' ')[2]
        Len = Len.split(':')[1]
        TransDecoderLengths.append(int(Len))

        temp = line.split('_')[:2]
        ContigName = '_'.join(temp)

        Cords = line.split(':')[9].split('(')[:1]
        Cords = str(Cords).strip('[]')

        CDS = TransDecoderInput.next().strip('\n')





        Data = Cords + ',' + Len + ',' + CDS
        Median =  numpy.median(TransDecoderLengths)
        Ninput = sorted(TransDecoderLengths)
        N50 = N50_(Ninput)
        TransDecoderORFs.update({ContigName:Data})
        TransDecoderOutput.write(ContigName+','+Data+','+str(Median)+','+str(N50)+'\n')
#
#
# # list = []
# # list.append(91)
# # list.append(77)
# # list.append(70)
# # list.append(69)
# # N55 = sorted(list)[len(list)/2]
# N50 = sorted(TransDecoderLengths)[len(TransDecoderLengths)//2]
# print N50
# sortedd = sorted(TransDecoderLengths)
#
#
# with open('./Out.txt','w') as lengths:
#     for i in TransDecoderLengths:
#         lengths.write(str(i)+"\n")
#
#
# print N50_(sortedd)
#
# print numpy.median(TransDecoderLengths)
# #
# #
# for line in fileinput.input('./Sample_Data.txt.lst'):
#     if ">gene_" in line:
#         Len = line.split('|')[2].split('_')[0]
#         temp = int(Len)
#         Len = temp/3
#         temp = line.split('|')[5].split('\t')[1]
#         ContigName = temp.split('_',4)[:4]
#         ContigName = '_'.join(ContigName)
#         L = line.split('|')[4]#.split(':')[:2]
#         Strand = line.split('|')[3]
#         R = line.split('|')[5].split('\t')[0]
#         Cords = L + '-' + R + '(' + Strand + ')'
#         CDS = ''
#         l_no = fileinput.lineno()+1
#
#         while linecache.getline('./Sample_Data.txt.lst',l_no) != '\n':
#             CDS += linecache.getline('./Sample_Data.txt.lst',l_no).strip('\n')
#             l_no += 1
#         Data = Cords + ',' + str(Len) + ',' + CDS
#         MetaGeneMarkORFs.update({ContigName: Data})
#         MetaGeneMarkOutput.write(ContigName+','+Data+'\n')
# #




# for line in TransDecoderBest:
#     if "gene" in line:
#         Contig = line.split('_')[:2]
#         Contig = '_'.join(Contig)
#         Contig +="_"
#         strand = line.split('\t')[:6]
#         start = line.split('\t')[3]
#         end = line.split('\t')[4]
#
#         Seq = ""
#
#         for line in fileinput.input('./Sample_Data.txt'):
#             if Contig in line:
#                 l_no = fileinput.lineno()
#                 l_no += 1
#                 while ">NODE" not in linecache.getline('./Sample_Data.txt', l_no):
#                     Seq += linecache.getline('./Sample_Data.txt',l_no).strip('\n')
#                     l_no +=1
#                     if ">NODE" in linecache.getline('./Sample_Data.txt', l_no):
#                         break
#
#                 ORF = Seq[int(start):int(end)]
#                 TransDecoderBestOutput.write(Contig+","+Seq+"\n")
#                 fileinput.close()
#                 break
#
#
#









# # import re
# import HTSeq
#
# #
# # # #
# RNA_Alignments = open('./RNA_Alignments.txt', mode='rb')
# DNA_Alignments = open('./Contigs_1_GFF.gff', mode='rb')
#
#
#
#
# output = open('./Unaligned_Contigs_RNA_Aligned.txt', mode='wb')
#
# Interesting_Contigs = {}
#
#
#
#
# RNA_Nodes = {}
# DNA_Nodes = {}
#
#
# for line in RNA_Alignments:
#     if "NODE" in line:
#
#         line = line.split("\t")
#         Alignments = line[2]
#
#         if int(Alignments) != 0:
#             Node = line[0]
#             Node = Node.split("_", 2)
#             Node = Node[1]
#             RNA_Nodes.update({str(Node):0})
#
#
#
# for line in DNA_Alignments:
#     if "NODE" in line:
#
#         line = line.split("\t")
#         Node = line[0]
#         Node = Node.split("_", 2)  # Strip on nth occurance
#         Node = Node[1]
#         DNA_Nodes.update({str(Node):0})
#
#
# for node in RNA_Nodes:
#     if node not in DNA_Nodes:
#         output.write(str(node)+"\n")
#         Interesting_Contigs.update({str(node): 0})
#
#
#
# output.close()
#
# RNA_Aligned = open('./Unaligned_Contigs_RNA_Aligned.txt', mode='rb')
#
# output = open('./Interesting_Contigs.txt', mode='wb')
#
#
# RNA_Alignments = []
#
#
# for node in RNA_Aligned:
#     node = node.replace("\n","")
#     RNA_Alignments.append(node)
#
#
# for read in HTSeq.FastaReader("Contigs_1.fa"):
#
#     name = str(read.name)
#     name = name.split("_", 2)  # Strip on nth occurance
#     name = name[1]
#
#
#
#     if name in RNA_Alignments:
#             output.write(str(read.name)+"\n"+str(read)+"\n")
#
#             continue
#
#
#
#
#
#
#
#
#
