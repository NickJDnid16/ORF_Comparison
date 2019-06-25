import sys
import pdb

def revComp(seq):
    seq = seq.upper()  # Makes seq uppercase
    seq = seq[::-1]  # Reverses seq
    seq = seq.replace('A', 't')  # Replace ACGT with lowercase complement
    seq = seq.replace('C', 'g')
    seq = seq.replace('G', 'c')
    seq = seq.replace('T', 'a')
    seq = seq.upper()  # Make seq uppercase again

    isitempty = seq
    isitempty = isitempty.replace('A', "")
    isitempty = isitempty.replace('C', "")
    isitempty = isitempty.replace('G', "")
    isitempty = isitempty.replace('T', "")
    isitempty = isitempty.replace('N', "")
    if isitempty != "":
        print ("Careful, improper characters!")

    return seq


#####################################
#  Iterative method
#####################################

def revCompIterative(watson):
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    watson = watson.upper()
    watsonrev = watson[::-1]

    crick = ""

    print (len(watson))
    print (len(watsonrev))
    for nt in watsonrev:
        crick += complements[nt]
    print (len(crick))
    return crick
