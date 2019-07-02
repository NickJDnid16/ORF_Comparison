import collections
import click
import csv
from DNA_Reverse_Compliment import revCompIterative


import numpy as np
from OrfFinderComparison import orfComparison


@click.command()
@click.option('--tool', default='GFF', help='Which tool/format would you like to analyse?')
@click.option('--tool_parameters', default='NULL', help='Parameters of tool if needed, e.g. prediction model')
@click.option('--input_to_analyse', help='Location of input file')
@click.option('--input_fasta', help='Location of fasta')
@click.option('--input_gff', help='Location of gff')

def comparator(tool,tool_parameters,input_to_analyse,input_fasta,input_gff):


    Genome = ""
    with open(input_fasta, 'r') as genome:
        for line in genome:
            line = line.replace("\n","")
            if ">" not in line:
                Genome += str(line)


##############################################
    Genes = collections.OrderedDict()

    count = 0
    with open(input_gff,'r') as genome_gff:
        for line in genome_gff:
            if "Chromosome	ena	CDS" in line:
                count += 1
                Start = int(line.split('\t')[3])
                Stop = int(line.split('\t')[4])
                Gene = str(Start)+','+str(Stop)
                Genes.update({count:Gene})
    print(len(Genes))
    #
# ##################################

    from .tool import tool




# Start_Codons = collections.OrderedDict()
# Stop_Codons = collections.OrderedDict()
# Prev_Stop = 0
# OL = False
# Genome_Size = len(Genome)
# Genome_rev = revCompIterative(Genome)
#
# for line in AugustusInput:
#     if "Chromosome	AUGUSTUS	CDS" in line :
#         Start = int(line.split('\t')[3])
#         Stop = int(line.split('\t')[4])
#         if Start < Prev_Stop and OL == False:
#             print ("Overlapping")
#             OL = True
#         Prev_Stop = Stop
#         reverse = '-'
#         forward = '+'
#         if '-' in line.split('\t')[6]: #Reverse Compliment starts and stops to confirm to our deffinition
#             #Switched to match Sense Strand
#             r_start = Genome_Size - Stop
#             r_stop = Genome_Size - Start
#             po = str(r_start)+','+str(Stop)
#             Start_Codons.update({po:reverse})
#             po = str(Start) + ',' + str(r_stop)
#             Stop_Codons.update({po: reverse})
#         elif '+' in line.split('\t')[6]:
#             po = str(Start)+','+str(Stop)
#
#             Start_Codons.update({po:forward})
#             Stop_Codons.update({po: forward})
#
#         AugustusORFs.update({Start:Stop})
#
#
#
#
#
# Output,Start_Precision,Stop_Precision = orfComparison(Genes,AugustusORFs,Start_Codons,Stop_Codons,Genome)
#
#
# # with open('./Augustus_Metrics.csv', mode='w') as outfile:
# #     outfile_writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
# #     outfile_writer.writerow
#
# AugustusOutput.write(str(Output)+'\n')
# AugustusOutput.write(str(Start_Precision)+'\n')
# AugustusOutput.write(str(Stop_Precision)+'\n')
#





















if __name__ == "__main__":
    comparator()








