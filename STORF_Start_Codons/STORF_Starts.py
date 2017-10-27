import re

FP_ORFs = open('../STORF/FP/FP_ORFs_Myco.txt', mode='rb')
FP_ORFs_RC = open('../STORF/FP/FP_ORFs_RC_Myco.txt', mode='rb')
FP_With_Starts = open('./FP_With_Starts_Myco.txt', mode='wb')


A = "ATG"
G = "GTG"
T = "TTG"
AT = "ATT"
C = "CTG"



Count = 0


for line in FP_ORFs:
    line = line.split(',')
    line.remove('\n')
    if int(line[2]) - int(line[1]) >100:

        Codons = re.findall('...', line[0])
        Starts = []
        #for Codon in Codons:

        if A in Codons:
            Starts.append(Codons.index(A))
        if G in Codons:
            Starts.append(Codons.index(G))
        if T in Codons:
            Starts.append(Codons.index(T))
        if AT in Codons:
            Starts.append(Codons.index(AT))
        if C in Codons:
            Starts.append(Codons.index(C))
        if len(Starts) > 0:
            FP_With_Starts.write('>'+'Myco:'+str(Count)+'\n')
            FP_With_Starts.write(str(', '.join(line)) + ',' + str(min(Starts)) + '\n')
            Count +=1

#################################################################
Count = 0
for line in FP_ORFs_RC:
    line = line.split(',')
    line.remove('\n')
    if int(line[2]) - int(line[1]) >100:

        Codons = re.findall('...', line[0])
        Starts = []
        #for Codon in Codons:

        if A in Codons:
            Starts.append(Codons.index(A))
        if G in Codons:
            Starts.append(Codons.index(G))
        if T in Codons:
            Starts.append(Codons.index(T))
        if AT in Codons:
            Starts.append(Codons.index(AT))
        if C in Codons:
            Starts.append(Codons.index(C))

        if len(Starts) > 0:
            FP_With_Starts.write('>'+'Myco_Rev:'+str(Count)+'\n')
            FP_With_Starts.write(str(', '.join(line)) + ',' + str(min(Starts)) + '\n')
            Count +=1
