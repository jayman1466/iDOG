# iDOG

inter-Domain Operon Generator

## Getting Started

This growing suite of tools enables biologists to generate genetic constructs capable of expressing in highly diverse microbial chassis, as outlined in [Jaymin's Paper]. A web interface for these tools is available at [Website]. This github repository provides local command line versions of the RBS generator and Codon Optimizer.  

### Prerequisites

These scripts are designed to work on UNIX systems, though they have only been tested in Linux.
Requirements:

* python (tested on version 3.7)
* python packages:
  * biopython (tested on version 1.7.6)
  * pillow (tested on version 6.2.0)
* Vienna RNA Package (tested on version 2.4.14) **Individual programs need to be in path**
* A compiled TranstermHP is already provided. This should work for x86 Linux builds. But for arm64 processors, you'll have to recompile from [source](http://transterm.ccb.jhu.edu/).

### Create a Ribosome Binding Site

Within create_RBS.py, the function create_RBS will design a custom broad host-range RBS

The basic command is as follows:
```
create_RBS(orf = "", target_tir = "")
```
where:
* orf is the Nucleotide Sequence of the Open Reading Frame, which must begin with NTG and must be >34bp long
* target_tir is the target Translation Initiation Rate as an integer from 1-999,999. 50,000 is a starting point for overexpression.

optional parameters:
* rRNA_sequence: nucleotide sequence of the presumed 3' end of 16S rRNA. Default: "ACCUCCUUA"
* detailed_output: = True/False. Whether you would like RNA folding plots of the RBS to be outputted. Default: True
* save_url: Directory location where you'd like the detailed output to be stored. Default: ""

### Codon Optimize a Gene

Within codon_opt.py, the function codon_opt with recode your gene for broad host-range expression as well as create a companion RBS for the gene

The basic command is as follows:
```
codon_opt(orf_aa = "", target_tir = "")
```
where:
* orf_aa is the Amino Acid Sequence of the Open Reading Frame, which must begin with M and must be >13AA lowest_energy
* target_tir is the target Translation Initiation Rate as an integer from 1-999,999. 50,000 is a starting point for overexpression.

optional parameters:
* int_RBS: = True/False. Remove internal RBSs? Default: True
* terminators: = True/False. Remove internal rho-independent transcriptional terminators? Default: True
* species: Base species to use for codon distribution. Default: "escherichia_coli". See section below to add additional organisms 
* save_url: Directory location where you'd like the detailed output to be stored. Default ""
* detailed_output: = True/False. Whether you would like RNA folding plots of the RBS, and %GC + RNA Structure plots of the CDS to be outputted. Default: True)

## Author
Jaymin Patel (jayman1466@gmail.com)
