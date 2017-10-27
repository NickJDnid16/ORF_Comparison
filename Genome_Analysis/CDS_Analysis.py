Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.gff', mode='rb')
Lengths= []
Longest = ""
LongNum = 0
for line in Input:
    if "Chromosome	ena	gene" in line:
        Length = int(line.split()[4]) - int(line.split()[3])
        Lengths.append(Length)
        if Length > LongNum:
            LongNum = Length
            Longest = line
print sum(Lengths) / len(Lengths)
print min(Lengths)
print max(Lengths)

print LongNum
print Longest


GInput = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.fa', mode='rb')


Genome = []
gene = ''
for line in GInput:
    line = line.replace("\n","")
    if ">" not in line:
        gene += str(line)
    else:
        Genome.append(gene)
        gene = ''
ATG = 0
GTG = 0
TTG = 0
ATT = 0
CTG = 0
Other = 0


for gene in Genome:
    codon = gene[:3]
    if len(gene) != 0:
        if codon == 'ATG':
            ATG += 1
        elif codon == 'GTG':
            GTG += 1
        elif codon == 'TTG':
            TTG += 1
        elif codon == 'ATT':
            ATT += 1
        elif codon == 'CTG':
            CTG += 1
        else:
            Other += 1
            print gene



ATG_Percentage = float((ATG) * float(100) / float(3737))
GTG_Percentage = float((GTG) * float(100) / float(3737))
TTG_Percentage = float((TTG) * float(100) / float(3737))
ATT_Percentage = float((ATT) * float(100) / float(3737))
CTG_Percentage = float((CTG) * float(100) / float(3737))
Other_Percentage = float(Other) * float(100) / float(3737)
print format(ATG_Percentage, '.2f')
print format(GTG_Percentage, '.2f')
print format(TTG_Percentage, '.2f')
print format(ATT_Percentage, '.2f')
print format(CTG_Percentage, '.2f')
print format(Other_Percentage, '.2f')


TAG = 0
TAA = 0
TGA = 0
count = 0
print "STOPS"
for gene in Genome:
    codon = gene[-3:]
    if len(gene) != 0:
        count +=1
        if codon == 'TAG':
            TAG += 1
        elif codon == 'TAA':
            TAA += 1
        elif codon == 'TGA':
            TGA += 1
        else:
            print codon
print "TAG = "+str(TAG)
print "TAA = "+str(TAA)
print "TGA = "+str(TGA)
print "NumGenes = "+str(count)






















