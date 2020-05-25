#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 12:47:20 2020

@author: jayman1466
"""

from RNAfold_wrapper import RNAfold, RNAduplex, RNAcofold
from dG_tot_calc import dG_tot_calc

import os
from os import path
import subprocess

import matplotlib.pyplot as plt
from PIL import Image
from glob import glob

#this script expects the RBS to be 35bp for now.
def folding_after_rRNA_binding(RBS_seq, orf_nu, rRNA_sequence = "ACCUCCUUA", save_url = ""):
    orf_seq = orf_nu[0:3].upper() + orf_nu[3:35].lower() #we only want 35bp of the orf and only have the start codon uppercase
    rRNA_sequence = rRNA_sequence.lower()
    #truncate RBS to 35bp
    if len(RBS_seq)>34:
        RBS_seq = RBS_seq[-35:]
    RBS_seq = RBS_seq.lower()
    save_url = save_url

    #use dG_tot_calc to determine how the rRNA binds to the mRNA
    dG = dG_tot_calc(RBS_seq,orf_seq,rRNA_sequence,0)

    #pull out the bases on the mRNA upstream of the rRNA binding site and fold. The standby site is 4bp immediately upsteam of the rRNA binding site, so we won't fold that
    five_prime_base_number = dG[7]
    five_prime_RNA = RBS_seq[0:five_prime_base_number - 4]

    #pull out the bases of the rRNA binding site
    rRNA_length = len(rRNA_sequence)
    rRNA_binding_RNA = RBS_seq[five_prime_base_number:five_prime_base_number + rRNA_length]

    #pull out the bases downstream of the ribosome footprint (15bp) on the orf
    downstream_RNA = orf_seq[15:]

    #fold each of the individual components
    five_prime_RNA_folding = RNAfold(five_prime_RNA).folding[0]
    downstream_RNA_folding = RNAfold(downstream_RNA).folding[0]

    RNAduplex_input = '{}\n{}'.format(rRNA_binding_RNA,rRNA_sequence)
    RNAduplex_output = RNAduplex(RNAduplex_input)[0].split()
    RNAduplex_output_folding = RNAduplex_output[0].split('&')
    RNAduplex_output_mRNAcoords = RNAduplex_output[1].split(',')
    RNAduplex_output_rRNAcoords = RNAduplex_output[3].split(',')

    #fill in the unbound dots in the RNAduplex output
    mRNA_5_append = ''
    mRNA_3_append = ''
    rRNA_5_append = ''
    rRNA_3_append = ''
    if int(RNAduplex_output_mRNAcoords[0]) > 1:
        mRNA_5_append = '.' * (int(RNAduplex_output_mRNAcoords[0]) - 1)
    if int(RNAduplex_output_mRNAcoords[1]) < rRNA_length:
        mRNA_3_append = '.' * (rRNA_length - int(RNAduplex_output_mRNAcoords[1]))
    if int(RNAduplex_output_rRNAcoords[0]) > 1:
        rRNA_5_append = '.' * (int(RNAduplex_output_rRNAcoords[0]) - 1)
    if int(RNAduplex_output_rRNAcoords[1]) < rRNA_length:
        rRNA_3_append = '.' * (rRNA_length - int(RNAduplex_output_rRNAcoords[1]))

    #get the mRNA componenet of the duplex output.
    rRNA_binding_mRNA_folding = mRNA_5_append + RNAduplex_output_folding[0] + mRNA_3_append
    rRNA_folding = rRNA_5_append + RNAduplex_output_folding[1] + rRNA_3_append

    ribosome_footprint = "x" * 15
    standby_site = "x" *4

    #get the spacing between the rRNA binding site and the start codon
    spacing_length = len(RBS_seq) - len(five_prime_RNA) - len(standby_site) - rRNA_length
    spacing = '.' * spacing_length

    #compile the final folded mRNA for input into RNAcofold
    mRNA_folding = five_prime_RNA_folding + standby_site + rRNA_binding_mRNA_folding + spacing + ribosome_footprint + downstream_RNA_folding
    RNAcofold_input_folding = mRNA_folding + '&' + rRNA_folding
    #apply the constraint rules for RNAcofold
    RNAcofold_input_folding = RNAcofold_input_folding.replace('(','<')
    RNAcofold_input_folding = RNAcofold_input_folding.replace(')','>')

    #compile the input for RNAcofold
    RNAcofold_input_filename_post = ">post_binding"
    RNAcofold_input_sequence = RBS_seq + orf_seq + '&' + rRNA_sequence

    RNAcofold_input_post = '{}\n{}\n{}'.format(RNAcofold_input_filename_post, RNAcofold_input_sequence, RNAcofold_input_folding)

    #also compile the input for pre-binding folding
    pre_input_filename = ">pre_binding"
    pre_input_sequence = RBS_seq + orf_seq

    RNAcofold_input_pre = '{}\n{}'.format(pre_input_filename, pre_input_sequence)

    #change the working directory based on save_url
    current_working_dir = os.getcwd()
    #create the directory if it doesn't specified and doesn't exist
    directory_name = save_url
    if directory_name is not "":
        if not path.exists(directory_name):
            os.makedirs(directory_name, exist_ok=True) #remember to add mode later
        os.chdir(directory_name)

    #do the final foldings
    RNAcofold(RNAcofold_input_pre)
    RNAcofold(RNAcofold_input_post)

    #run relplot.pl to colorize the pre plot by accessibility
    replot_input = ["perl", current_working_dir + "/relplot.pl", "-a", "pre_binding_ss.ps", "pre_binding_dp.ps"]

    with open('pre.ps', "w") as outfile:
        subprocess.run(replot_input, stdout=outfile)

    #remove extra files
    os.remove("post_binding_dp.ps")
    os.remove("pre_binding_ss.ps")
    os.remove("pre_binding_dp.ps")

    #for some reason svg fails so have to use pil to convert

    psfiles = glob('*.ps')
    for u in psfiles:
        out = u.replace('ps','png')
        img=Image.open(u)
        img.save(out, dpi=(600,600))

    os.remove("post_binding_ss.ps")
    os.remove("pre.ps")

    #change the directory back
    os.chdir(current_working_dir)


#Sliding window of GC content
def GC_plot(RBS_seq, orf_nu, sliding_window = 30, step_size = 15, save_url = ""):
    this_RBS = RBS_seq
    this_orf = orf_nu
    this_mRNA = this_RBS + this_orf
    sliding_window = sliding_window #30 appears to be a decent value
    step_size = step_size #15 appears to be a good number
    save_url = save_url

    #trim the ends to avoid problems with the sliding window
    margin = sliding_window//2 + 1

    if margin < len(this_RBS):
        x = margin
    else:
        x=len(RBS_seq)


    GC_content = []
    nucleotide_position = []


    #run a sliding window of GC contect, Starting from the start codon. Go in steps of 10
    while x < len(this_mRNA) - margin:
        #get the subsequence
        this_query_string = this_mRNA[x - (sliding_window//2) : x+ (sliding_window//2)]

        #GC content of subsequence
        this_GC_count = this_query_string.count('G') + this_query_string.count('C') + this_query_string.count('g') + this_query_string.count('c')
        this_GC_content = this_GC_count/sliding_window

        #update the arrays
        GC_content.append(this_GC_content)
        nucleotide_position.append(x-len(this_RBS)+1)

        #Step forward
        x = x+step_size

    #plot it
    fig, ax = plt.subplots()

    ax.plot(nucleotide_position, GC_content, linewidth=0)
    ax.set(ylim=(0,1), xlim=(1-len(RBS_seq),len(this_orf)))
    ax.spines["top"].set_visible(False)
    #ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    #ax.spines["left"].set_visible(False)
    ax.tick_params(axis='y', which='major', pad=2, labelsize=12)
    ax.tick_params(axis='x', which='major', pad=2, labelsize=12)
    ax.set_ylabel('GC Content', fontsize=16, labelpad=6)
    ax.set_xlabel('Nucleotide Position', fontsize=16, labelpad=6)
    plt.setp(ax.spines.values(), linewidth=2)
    plt.tight_layout()
    plt.fill_between(nucleotide_position, GC_content, color=(0.38, 0.65, 0.87))
    ax.axvline(x=0, color="black", linestyle=':')

    #create the directory if specified and doesn't exist
    directory_name = save_url
    if directory_name is not "":
        if not path.exists(directory_name):
            os.makedirs(directory_name, exist_ok=True) #remember to add mode later
    plt.savefig("{}GC_content".format(directory_name))


#Sliding window of RNA structure
def structure_plot(RBS_seq, orf_nu, sliding_window = 40, step_size = 15, save_url = ""):
    this_RBS = RBS_seq
    this_orf = orf_nu
    this_mRNA = this_RBS + this_orf
    sliding_window = sliding_window #40 appears to be a decent value
    step_size = step_size #15 appears to be a good number
    save_url = save_url

    #trim the ends to avoid problems with the sliding window
    margin = sliding_window//2 + 1

    if margin < len(this_RBS):
        x = margin
    else:
        x=len(RBS_seq)

    folding_energies = []
    nucleotide_position = []


    #run a sliding window of GC contect, Starting from the start codon. Go in steps of 10
    while x < len(this_mRNA) - margin:
        #get the subsequence
        this_query_string = this_mRNA[x - (sliding_window//2) : x+ (sliding_window//2)]

        #GC content of subsequence
        this_folding_energy = RNAfold(this_query_string).folding[1]

        #update the arrays
        folding_energies.append(this_folding_energy)
        nucleotide_position.append(x-len(this_RBS)+1)

        #Step forward
        x = x+step_size

    #plot it
    fig, ax = plt.subplots()

    ax.plot(nucleotide_position, folding_energies, linewidth=0)
    ax.set(ylim=(0,-20), xlim=(1-len(RBS_seq),len(this_orf)))
    ax.spines["top"].set_visible(False)
    #ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    #ax.spines["left"].set_visible(False)
    ax.tick_params(axis='y', which='major', pad=2, labelsize=12)
    ax.tick_params(axis='x', which='major', pad=2, labelsize=12)
    ax.set_ylabel('Folding Energy', fontsize=16, labelpad=6)
    ax.set_xlabel('Nucleotide Position', fontsize=16, labelpad=6)
    plt.setp(ax.spines.values(), linewidth=2)
    plt.tight_layout()
    plt.fill_between(nucleotide_position, folding_energies, color=(0.38, 0.65, 0.87))
    ax.axvline(x=0, color="black", linestyle=':')


    #create the directory if it doesn't exist
    directory_name = save_url
    if directory_name is not "":
        if not path.exists(directory_name):
            os.makedirs(directory_name, exist_ok=True) #remember to add mode later
    plt.savefig("{}folding_energy".format(directory_name))
