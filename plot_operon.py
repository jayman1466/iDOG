# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 14:31:39 2017

@author: jaymin patel
"""
import dnaplotlib as dpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
import subprocess
from os import path

#plot generic layout of operon
def plot_operon_structure(design,fig_width,directory_name=""):
    # Redender the DNA
    fig = plt.figure(figsize=(fig_width/85,0.5))
    gs = gridspec.GridSpec(1, 1)
    ax_dna = plt.subplot(gs[0])

    dr = dpl.DNARenderer()
    part_renderers = dr.SBOL_part_renderers()
    start, end = dr.renderDNA(ax_dna, design, part_renderers)

    # Set bounds and display options for the DNA axis
    dna_len = end-start

    ax_dna.set_xlim([start-20, end+20])
    ax_dna.set_ylim([-15,20])
    ax_dna.plot([start-20,end+20], [0,0], color=(0,0,0), linewidth=1.0, zorder=1)
    ax_dna.axis('off')
    plt.subplots_adjust(hspace=.08, left=.01, right=.99, top=0.99, bottom=0.02)

    #create the directory if it doesn't exist
    if not path.exists(directory_name):
        os.makedirs(directory_name, exist_ok=True) #remember to add mode later
    plt.savefig("{}sbol_plot".format(directory_name), dpi = 300)
    plt.close('all')






#plot NuPoP nucleosome occupancy
def plot_nupop(design_nupop,full_sequence,directory_name=""):

    #wrap promoter with 420bp buffer sequence on each side
    pre = "aaaagactctaacaaaatagcaaatttcgtcaaaaatgctaagaaataggttattactgagtagtatttatttaagtattgtttgtgcacttgcctgcaagccttttgaaaagcaagcataaaagatctaaacataaaatctgtaaaataacaagatgtaaagataatgctaaatcatttggctttttgattgattgtacaggaaaatatacatcgcagggggttgacttttaccatttcaccgcaatggaatcaaacttgttgaagagaatgttcacaggcgcatacgctacaatgacccgattcttgctagccgaattccagtcaggctgctagcaccagagctacgtgaccgcaggactagctccagctgagcgacacgcgaacaggtgtagagtcagtcgtgctgcaaggtcgcac"
    post = "TTAAAGAGGAGAAAGGTACCATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCACTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATG"
    buffered_sequence = pre + full_sequence + post

    #run NuPoP on R
    NuPoP_out = subprocess.check_output(["Rscript", "NuPoP_Rscript.R", buffered_sequence])
    nuc_occu = NuPoP_out.decode("utf-8").split()

    nuc_occu = list(map(float, nuc_occu))

    operon_length = len(full_sequence)

    # Create the figure and all axes to draw to (width of figure should be size of operon divided by 200)
    fig = plt.figure(figsize=(operon_length/400,1.5))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.5,0.3])
    ax_dna = plt.subplot(gs[0])

    # Render the DNA
    dr = dpl.DNARenderer(scale=1, linewidth=0.9)
    start, end = dr.renderDNA(ax_dna, design_nupop, dr.trace_part_renderers())
    # Set bounds and display options for the DNA axis
    ax_dna.set_xlim([0,operon_length])
    ax_dna.set_ylim([-8,8])
    ax_dna.plot([0,operon_length], [0,0], color=(0,0,0), linewidth=1.0, zorder=1)
    ax_dna.axis('off')

    #plot the NuPoP output
    nuc_coords = range(-420,operon_length+420)
    ax_nupop = plt.subplot(gs[1])
    ax_nupop.plot(nuc_coords,nuc_occu, linewidth=1)
    ax_nupop.set_xlim([0,operon_length])
    ax_nupop.set_ylim([0,1.1])
    ax_nupop.set_xlabel('Nucleotide Position', fontsize=8, labelpad=6)

    ax_nupop.spines["top"].set_visible(False)
    #ax.spines["bottom"].set_visible(False)
    ax_nupop.spines["right"].set_visible(False)
    #ax.spines["left"].set_visible(False)

    ax_nupop.tick_params(axis='y', which='major', pad=2, labelsize=8)
    ax_nupop.tick_params(axis='x', which='major', pad=2, labelsize=8)

    plt.setp(ax_nupop.spines.values(), linewidth=1)
    plt.fill_between(nuc_coords, nuc_occu, color=(0.38, 0.65, 0.87))
    plt.subplots_adjust(hspace=.08, left=.03, right=.99, top=0.99, bottom=0.30)

    #create the directory if it doesn't exist
    if not path.exists(directory_name):
        os.makedirs(directory_name, exist_ok=True) #remember to add mode later
    plt.savefig("{}nupop_plot".format(directory_name), dpi = 300)
