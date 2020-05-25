# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 16:51:46 2015

@author: Jaymin
"""
import random
import csv
import ast
import itertools

def orf_recode(sequence,species):
    sequence = sequence.upper()
    species = species
    
    #indicate which percent of the N-terminal 12 codons should be pulled from the low_struct dictionary. 2=50% 3=66% 4=75%
    low_struct_frac = 3
    
    #pull dictionary for reverse translation with weights adjusted by codon usage in highly expressed genes of species (removed TTG, and GTG and TTA codons)
    #open file containing codon table and rRNA_sequence for species
    with open("codon_tables/{}.csv".format(species)) as csvfile:
        reader = csv.reader(csvfile)
        
        #output as a dictionary linking to codons and weights for each amino acid. Don't read the last line, which is the rRNA_sequence
        codon_table = {row[0] : [ast.literal_eval(row[1]) , ast.literal_eval(row[2])] for row in itertools.islice(reader,20)}
    
    low_struct = {"A":"GCA","C":"TGT","D":"GAT","E":"GAA","F":"TTT","G":"GGA","H":"CAC","I":"ATA","K":"AAA","L":"CTT","M":"ATG","N":"AAT","P":"CCC","Q":"CAA","R":"AGG","S":"TCA","T":"ACA","V":"GTA","Y":"TAT","W":"TGG","*":"TAA"};  
    
    output=''

    i=1
    for c in sequence:
        if i<13:
            
            k=random.randint(1,low_struct_frac)
            if k==1:
                output+=random.choices(codon_table[c][0],codon_table[c][1])[0]            
            else:
                output+=low_struct[c]        
        
        else:
            output+=random.choices(codon_table[c][0],codon_table[c][1])[0]
        i=i+1
    output+='TAA'
    return output
