import numpy as np
import sys
import collections
from DNA_Reverse_Compliment import revCompIterative


def orfComparison(Genes,ORFs,Start_Codons,Genome):
    ###################################################





    ORFsFiltered = collections.OrderedDict()


    ####################################################
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)
    Under_Predicted = []
    inFrame = 0
    GCovered = []
    Start_Precision = []
    Stop_Precision = []
    Perfect_Matches = 0
    Perfect_Starts = 0
    Perfect_Stops = 0
    Lengths = []
    ###################################################

    Switch = 0

    #####
    FP_ORFs = open('./FP/FP_ORFs_Myco.txt', mode='wb')

    #####
    Gene_A = np.zeros((Genome_Size), dtype=np.int)
    SingleORF = np.zeros((Genome_Size), dtype=np.int)





    for positions in Genes.itervalues():
        G_Start = int(positions.split(',')[0])
        G_Stop = int(positions.split(',')[1])
        Gene_A[G_Start:G_Stop] = [1] * (G_Stop - G_Start)  # Changing all between the two positions to 1's

    for O_po, O_St in ORFs.items():
        A = np.zeros((Genome_Size), dtype=np.int)
        O_Start = int(O_po.split(',')[0])
        O_Stop = int(O_po.split(',')[1])
        length = O_Stop - O_Start

        A[O_Start:O_Stop] = [1] * (O_Stop - O_Start)

        Gene_AND = A & Gene_A

        tempy = np.count_nonzero(Gene_AND)

        if np.count_nonzero(Gene_AND) <= (length*0.5):
            SC_S = Start_Codons.get(O_po)
            if "+" in SC_S:
                FP_ORFs.write(">"+str(O_Start)+'+'+str(O_Stop)+'\n')
                FP_ORFs.write(Genome[O_Start:O_Stop+3]+'\n')

            elif "-" in SC_S:
                FP_ORFs.write(">" + str(O_Start) + '-' + str(O_Stop) + '\n')
                FP_ORFs.write(Genome_rev[O_Start:O_Stop+3]+'\n')




    sys.exit()
#
#
#     for positions in Genes.itervalues():
#         G_Start = int(positions.split(',')[0])
#         G_Stop = int(positions.split(',')[1])
#
#
#
#
#
#
#
#
#
#         Found = False
#
#         O_po, O_St = min(ORFs.items(), key=lambda (_, v): abs(v - G_Stop))  # Gets closest
#
#         O_Start = int(O_po.split(',')[0])
#         O_Stop = int(O_po.split(',')[1])
#         SingleGene = np.zeros((Genome_Size), dtype=np.int)
#         SingleORF = np.zeros((Genome_Size), dtype=np.int)
#
#         gg = []
#         oo = []
#
#         SingleGene[G_Start:G_Stop] = [1] * (G_Stop - G_Start)  # Changing all between the two positions to 1's
#         gg.append(abs(G_Stop - G_Start))
#
#
#         SingleORF[O_Start:O_Stop] = [1] * (O_Stop - O_Start)  # Changing all between the two positions to 1's
#         oo.append(abs(O_Stop - O_Start))
#         AND = SingleGene & SingleORF
#         AND_Count_Stop = np.count_nonzero(AND)
#
#         O_po, O_St = min(ORFs.items(), key=lambda (_, v): abs(v - G_Start))  # Gets closest
#         O_Start = int(O_po.split(',')[0])
#         O_Stop = int(O_po.split(',')[1])
#
#
#         SingleGene = np.zeros((Genome_Size), dtype=np.int)
#         SingleORF = np.zeros((Genome_Size), dtype=np.int)
#
#         gg = []
#         oo = []
#
#         SingleGene[G_Start:G_Stop] = [1] * (G_Stop - G_Start)  # Changing all between the two positions to 1's
#         gg.append(abs(G_Stop - G_Start))
#
#         SingleORF[O_Start:O_Stop] = [1] * (O_Stop - O_Start)  # Changing all between the two positions to 1's
#         oo.append(abs(O_Stop - O_Start))
#         AND = SingleGene & SingleORF
#         AND_Count_Start = np.count_nonzero(AND)
#
#         if AND_Count_Stop >= AND_Count_Start:
#             O_po, O_St = min(ORFs.items(), key=lambda (_, v): abs(v - G_Stop))  # Gets closest
#             O_Start = int(O_po.split(',')[0])
#             O_Stop = int(O_po.split(',')[1])
#         else:
#             Switch +=1
#             print "Switch"
#             print Switch
#
#
#
#
#
#
#
#
#
#
#
#
# #############################
#         if G_Start == O_Start: ##Perfect Starts
#             Perfect_Starts += 1
#
#         if G_Stop == O_Stop: ##Perfect Stops
#             Perfect_Stops += 1
#
#         if G_Start == O_Start and G_Stop == O_Stop: ##Perfect Match
#             ORFsFiltered.update({O_Start: O_Stop})
#             GCovered.append(G_Stop)
#
#
#
#
#             Perfect_Matches += 1
#             Found = True
#
#         elif O_Start <= G_Start and O_Stop > G_Start:
#             ORFsFiltered.update({O_Start: O_Stop})
#             Found = True
#             GCovered.append(G_Stop)
#
#
#
#         elif O_Start >= G_Start and O_Start < G_Stop:
#             ORFsFiltered.update({O_Start: O_Stop})
#             Found = True
#             GCovered.append(G_Stop)
#
#
#         elif O_Start < G_Start and O_Stop > G_Stop:
#             ORFsFiltered.update({O_Start: O_Stop})
#             Found = True
#             GCovered.append(G_Stop)
#
#
#
#         if Found == True:
#             Start_Precision.append(O_Start - G_Start)
#             Stop_Precision.append((O_Stop - G_Stop))
#             l = (int(O_Stop) - int(O_Start))
#             if l < 0:
#                 print "sssss"
#             Length = int(O_Stop) - int(O_Start)
#             if Length < 1:
#                 print "DERP"
#             Lengths.append(Length)
#             frame_tmp = O_Stop - G_Stop
#             if frame_tmp % 3 == 0 or frame_tmp == 0:
#                 inFrame += 1
#
#
#
#
#
#         ##################################
#     for O_po, O_St in ORFs.items():
#         O_Start = int(O_po.split(',')[0])
#         O_Stop = int(O_po.split(',')[1])
#         if O_Start not in ORFsFiltered and ORFsFiltered[O_Start] != O_Stop:
#             FP_ORFs.write(Genome[O_Start:O_Stop])
#             FP_ORFs.write('\n')
#             FP_ORFs_RC.write(Genome_rev[O_Start:O_Stop])
#             FP_ORFs_RC.write('\n')
#
#
#     ######################################
#
#     m = min(Lengths)
#     for positions in Genes.itervalues():
#         G_Start = int(positions.split(',')[0])
#         G_Stop = int(positions.split(',')[1])  # If There is an un-covered Gene
#         if G_Stop not in GCovered and G_Stop not in Under_Predicted:
#             Under_Predicted.append(G_Stop)
#
#     Over_Predicted = len(ORFs) - len(GCovered)
#
#     print len(GCovered)
#     print len(Under_Predicted)
#     print Over_Predicted
#
#     ######################################################
#
#     ORF_Num = len(ORFsFiltered)
#
#     Start_Codon = []
#
#     #####################################################
#
#     i = 0
#     tempy = []
#     for po, direction in Start_Codons.iteritems():
#
#         scStart = int(po.split(',')[0])
#         scStop = int(po.split(',')[1])
#         f_start, f_stop = min(ORFsFiltered.items(), key=lambda (_, v): abs(v - scStop)) #Gets closest
#
#
#
#         if str(scStop) == str(f_stop):# and str(scStart) == str(f_start):
#             i = i+1
#             if f_stop in tempy:
#                 print "derp"
#             tempy.append(f_stop)
#             if '-' in direction:
#                 Start_Codon.append(Genome_rev[scStart:scStart + 3])
#                 if "TAG"not in Codon or "TGA" not in Codon or "TAA" not in Codon:
#                     Codon = (Genome[scStart:scStart + 3])
#                     if ('TAG' not in Codon) and ('TGA' not in Codon) and ('TAA' not in Codon):
#                         print Codon
#                 #print Genome_rev[scStart:scStart+ 3]
#
#             elif '+' in direction:
#                 Start_Codon.append(Genome[scStart:scStart + 3])
#                 Codon =  (Genome[scStart:scStart+3])
#                 if ('TAG'not in Codon) and ('TGA' not in Codon) and ('TAA' not in Codon):
#                     print Codon
#                 #sys.exit()
#                 #print Genome[scStart-1:scStart-1 + 3]
#
#
#
#     print"This should match filtered"
#     print i
#     print"filtered"
#     print len(ORFsFiltered)
#
#
#
#     ########################################################
#
#     Median_Start_Precision = np.median(Start_Precision)
#     Median_Stop_Precision = np.median(Stop_Precision)
#     Quarterly_Start_Precison = np.percentile(Start_Precision,25)
#     Quarterly_Stop_Precison = np.percentile(Stop_Precision,25)
#
#     ######################################################
#     Min_Length = 9999999999999
#     Max_Length = 0
#     Average_Length = sum(Lengths)/len(Lengths)
#     for length in Lengths:
#         if length < Min_Length:
#             Min_Length = length
#         if length > Max_Length:
#             Max_Length = length
#     print Max_Length
#     print "Max"
#     #############################################
#     TAG = 0
#     TAA = 0
#     TGA = 0
#
#     Other = 0
#     print ORF_Num
#
#     for codon in Start_Codon:
#         if codon == 'TAG':
#             TAG += 1
#         elif codon == 'TAA':
#             TAA += 1
#         elif codon == 'TGA':
#             TGA += 1
#         else:
#             Other += 1
#             print codon
#
#
#
#     TAG_Percentage = float((TAG) * float(100) / float(ORF_Num))
#     TAA_Percentage = float((TAA) * float(100) / float(ORF_Num))
#     TGA_Percentage = float((TGA) * float(100) / float(ORF_Num))
#
#     Other_Percentage = float(Other) * float(100) / float(ORF_Num)
#
#
#
#
#     GeneArray = np.zeros((Genome_Size), dtype=np.int)
#     ORFArray = np.zeros((Genome_Size), dtype=np.int)
#
#
#
#     gg = []
#     oo = []
#     for positions in Genes.itervalues():
#         G_Start = int(positions.split(',')[0])
#         G_Stop = int(positions.split(',')[1])
#         GeneArray[G_Start:G_Stop] = [1] * (G_Stop - G_Start) #Changing all between the two positions to 1's
#         gg.append(abs(G_Stop-G_Start))
#
#     for O_Start, O_Stop in ORFsFiltered.iteritems():
#         ORFArray[O_Start:O_Stop] = [1] * (O_Stop - O_Start)  # Changing all between the two positions to 1's
#         oo.append(abs(O_Stop-O_Start))
#
#     print np.median(gg)
#     print np.median(oo)
#
#
#
#     print np.count_nonzero(GeneArray)
#     print np.count_nonzero(ORFArray)
#
#     print "Numbers"
#     AND = GeneArray & ORFArray
#     AND_Count = np.count_nonzero(AND)
#     print AND_Count
#
#     NOT_O = np.logical_not(ORFArray) + [0 for i in xrange(len(ORFArray))]
#     NOT_O_AND_Gene = NOT_O & GeneArray
#     NOT_O_AND_Gene_Count = np.count_nonzero(NOT_O_AND_Gene)
#     print NOT_O_AND_Gene_Count
#
#     NOT_G = np.logical_not(GeneArray) + [0 for i in xrange(len(GeneArray))]
#     NOT_G_AND_ORF = NOT_G & ORFArray
#     NOT_G_AND_ORF_Count = np.count_nonzero(NOT_G_AND_ORF)
#     print NOT_G_AND_ORF_Count
#
#     NOT_NOT = NOT_G & NOT_O
#     NOT_NOT_Count = np.count_nonzero(NOT_NOT)
#     print NOT_NOT_Count
#
#     print Genome_Size
#
#     print "Stats"
#     NT_TP = (float(AND_Count) * 100 / float(Genome_Size))
#     NT_FN = (float(NOT_O_AND_Gene_Count) * 100 / float(Genome_Size))
#     NT_FP = (float(NOT_G_AND_ORF_Count) * 100 / float(Genome_Size))
#     NT_TN = (float(NOT_NOT_Count) * 100 / float(Genome_Size))
#
#     print NT_TP
#     print NT_FN
#     print NT_FP
#     print NT_TN
#
#     print "Precision"
#     Precision = NT_TP/(NT_TP+NT_FP)
#     Recall = NT_TP/(NT_TP+NT_FN)
#     print Precision
#     print Recall
#
#
#
#
#     ############################################
#     #Equations for calculation of one-to-one statis-tics. TP stands for true positive, FN stands for false negative, FP stands for false positive.
#
#     TP = len(GCovered)
#     FN = len(Under_Predicted)
#     FP = Over_Predicted
#     TPFN = float(TP)+float(FN)
#     TPFP = float(TP)+float(FP)
#     Sensitivity = float(TP/TPFN)
#     Specificity = float(TP/TPFP)
#     print Sensitivity
#     print Specificity
#     #############################################
#     Frame_Percentage = inFrame*100/ORF_Num
#     #############################################
#     Output = [len(ORFs),ORF_Num,len(GCovered),Average_Length, Min_Length, Max_Length, Median_Start_Precision, Quarterly_Start_Precison, Median_Stop_Precision
#               , Quarterly_Stop_Precison, Perfect_Matches, Frame_Percentage
#               ,"Stops",format(TAG_Percentage, '.2f'), format(TAA_Percentage, '.2f'), format(TGA_Percentage, '.2f'),format(Other_Percentage, '.2f'),
#               len(Under_Predicted), Over_Predicted, format(Sensitivity, '.2f'), format(Specificity, '.2f'),"Stats",format(AND_Count, '.2f'),format(NT_TP, '.2f'),format(NOT_O_AND_Gene_Count, '.2f'),format(NT_FN, '.2f'),format(NOT_G_AND_ORF_Count, '.2f'),format(NT_FP, '.2f'),format(NOT_NOT_Count, '.2f'),format(NT_TN, '.2f'),format(Precision, '.2f'),format(Recall, '.2f')]
#     return Output,Start_Precision,Stop_Precision
#
#     ##################################################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
