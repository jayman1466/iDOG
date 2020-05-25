# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:54:30 2017

@author: root
"""

from dG_tot_calc import dG_tot_calc
import random
from operator import itemgetter
from detailed_output import folding_after_rRNA_binding

def create_RBS(orf, target_tir, rRNA_sequence = "ACCUCCUUA", detailed_output = True, save_url=""):

    max_attempts = 200


    #convert to RNA sequence
    orf = orf.upper()
    orf = orf.replace("T","U")

    #RBS is defined by N17,6x A/U,AGGAG, N4, AAA
    #build library of RBSs up to max_attempts. If any are within 1000 of target_tir, stick with that one

    RBSs = []

    i = 0
    while i < max_attempts:

        #Make the RBS sequence
        RBS = random.choice(["A","U","C"])

        j=1
        while j < 17:

            if RBS[j-1] == "U":
                #avoid making internal RBSs
                RBS += random.choice(["A","U","C"])
            else:
                RBS += random.choice(["A","U","G","C"])

            j = j+1

        while j < 23:
            RBS += random.choice(["A","U"])

            j = j+1

        RBS += "AGGAG"
        j = j+5

        while j<32:

            if RBS[j-1] == "U":
                #avoid making internal RBSs
                RBS += random.choice(["A","U","C"])
            else:
                RBS += random.choice(["A","U","G","C"])

            j = j+1

        RBS += "AAA"

        RBS_calc = dG_tot_calc(RBS,orf,rRNA_sequence,0)
        RBS_tir = RBS_calc[8]
        difference = abs(RBS_tir-target_tir)

        RBSs += [(difference,RBS,RBS_tir)]


        i = i+1

        #see if you reached within 10% of your target TIR
        if difference <= 0.1*target_tir:
            break
    #once you have all you RBS candidates, rank them by difference
    RBSs_sorted = sorted(RBSs, key = itemgetter(0))

    RBS_seq = RBSs_sorted[0][1]
    RBS_seq.replace("U","T")
    RBS_tir = RBSs_sorted[0][2]

    #Create the RNA folding plots if detailed_output_yes = True
    if detailed_output == True:
        folding_after_rRNA_binding(RBS_seq, orf, rRNA_sequence, save_url)

    #return sequence and TIR of RBS
    return RBS_seq,RBS_tir
