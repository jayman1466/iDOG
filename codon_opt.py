# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 14:31:39 2017

@author: jaymin patel
"""

import csv
import itertools
from dG_tot_calc import int_RBS_remove
from orf_recode import orf_recode
from create_RBS import create_RBS
from transterm_wrapper import transterm
from detailed_output import GC_plot, structure_plot

def codon_opt(orf_aa, target_tir, int_RBS = True, terminators = True, species = "escherichia_coli", save_url = "", detailed_output = True):
    #where to save any files created:
    save_url = save_url

    #orf sequence. Input is String
    orf_aa = orf_aa

    #target RBS strength. Input is Integer
    target_tir = target_tir

    #screen and remove internal RBSs? Input is Boolean
    int_RBS = int_RBS

    #screen and remove internal transcriptional terminators? Input is Boolean
    terminators = terminators

    #species optimizing for. Input is string
    species = species

    #open file containing codon table and rRNA_sequence for species. Open just the last line, which is the rRNA_sequence
    with open("codon_tables/{}.csv".format(species)) as csvfile:
        reader = csv.reader(csvfile)
        species_table = {row[0]:row[1] for row in itertools.islice(reader,20,21)}

    #pull rRNA_sequence
    rRNA_sequence = species_table['rRNA_sequence']

    #cycle through this pipeline until no terminators are found, if that option is set to 1
    while True:

        #recode sequence
        orf_nu = orf_recode(orf_aa,species)

        #convert to RNA
        orf_nu = orf_nu.replace("T","U")

        #create RBS
        (RBS_seq,RBS_tir) = create_RBS(orf = orf_nu[0:35], target_tir = target_tir, rRNA_sequence = rRNA_sequence, detailed_output = detailed_output, save_url = save_url)

        #remove internal RBSs
        if int_RBS == True:
            int_output = int_RBS_remove(RBS_seq,orf_nu,rRNA_sequence)
            orf_nu = int_output[1]
            RBS_seq = int_output[0]

        #check for internal terminators. Break while loop if clear of terminators
        if terminators == True:
            terminator_out = transterm(orf_nu)
            if terminator_out == "clear of terminators":
                break
        else:
            break

    #convert to DNA
    orf_nu = orf_nu.replace("U","T")
    RBS_seq = RBS_seq.replace("U","T")

    #Call the detailed output function
    if detailed_output == True:
        GC_plot(RBS_seq, orf_nu, 30, 15, save_url)
        structure_plot(RBS_seq, orf_nu, 40, 15, save_url)

    #return output
    return RBS_seq,orf_nu,RBS_tir
