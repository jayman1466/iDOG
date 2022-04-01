# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 16:51:46 2015

@author: Jaymin
"""
import random
import csv
import ast
import itertools
from RNAfold_wrapper import RNAfold

def orf_recode(sequence,species='escherichia_coli',restriction_sites=[]):
    sequence = sequence.upper()
    species = species
    
    #pull dictionary for reverse translation with weights adjusted by codon usage in highly expressed genes of species (removed TTG, and GTG and TTA codons)
    #open file containing codon table and rRNA_sequence for species
    with open("codon_tables/{}.csv".format(species)) as csvfile:
        reader = csv.reader(csvfile)
        
        #output as a dictionary linking to codons and weights for each amino acid. Don't read the last line, which is the rRNA_sequence
        codon_table = {row[0] : [ast.literal_eval(row[1]) , ast.literal_eval(row[2])] for row in itertools.islice(reader,20)}

    #split into N and C terminal sequences
    sequence_n = sequence[0:12]
    sequence_c = sequence[12:]
    
    while True:
        output=''
        for c in sequence_n:
            output+=random.choices(codon_table[c][0],codon_table[c][1])[0]                            
        #ensure that the first 30bp structure is >-3.5 kcal/mol
        this_folding_energy = RNAfold(output[:30]).folding[1]
        if this_folding_energy > -3.5:
            break
    for c in sequence_c:
        output+=random.choices(codon_table[c][0],codon_table[c][1])[0]
    output+='TAA'
    return output