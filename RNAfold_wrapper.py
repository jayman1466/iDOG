# -*- coding: utf-8 -*-
#!/usr/bin/env python3
'''A dumb wrapper for RNAfold from the Vienna RNA Package for Python 3.
Uses subprocess.Popen to open a stdin/stdout/stderr stream to RNAfold, with
the --noPS option, allowing inline use of RNAfold on RNA strings with output
wrapped in a simple object to separate sequence, energy and structure.
'''

'''The RNAfold adnd RNASubopt wrappers are entirely the work of github user cathalgarvey (RNAfoldWrap.py), for whom I am very greatful. Jaymin only contributed the wrapper for RNAduplex and Kinfold in this document''' 

import collections, math
import subprocess
import random

# Warn if Vienna not found.
from shutil import which
if not which("RNAfold"):
    from sys import stderr
    print("WARNING: Could not locate the RNAfold binary from the ViennaRNA package.", file=stderr)
del(which)

# Used for testing:
randseq = lambda n: ''.join([random.choice("ACGU") for x in range(0,n)])
RNAStructure = collections.namedtuple("RNAStructure",["structure","energy"])

class RNAFoldError(BaseException):
    'Used to wrap and raise messages from stderr when calling RNAfold.'
    pass

class RNAFoldOutput:
    'Wraps the two-line output from RNAfold and extracts sequence, structure and energy.'
    def __init__(self, rnafold_output):
        output_lines = rnafold_output.strip().splitlines()
        self.sequence = output_lines[0]
        structure = output_lines[1].split(None,1)[0].strip()
        energy = float(output_lines[1].rsplit("(",1)[1].strip("()").strip())
        self.folding = RNAStructure(structure, energy)

#This version is just for the dot-bracket notation output and no image is produced
def RNAfold(sequence, *args):
    # Note that RNAfold auto-converts "T" to "U" so this is unnecessary in Python
    # This behaviour can be overridden if calculation of DNA is desired.
    rnaf = subprocess.Popen(["RNAfold","--noPS"] + list(args),
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            # Universal Newlines effectively allows string IO.
                            universal_newlines=True)
    foldout, folderr = rnaf.communicate(sequence)
    if folderr:
        raise RNAFoldErr(folderr)
    output = RNAFoldOutput(foldout)
    #print("Debug: Energy was", output.folding.energy,"- RNAfold output was:\n",foldout) # Debug
    return output
    
class RNASuboptError(BaseException):
    'Used to wrap and raise messages from stderr when calling RNAsubopt.'
    pass

class RNASuboptOutput:
    'Wraps the muti-line output from RNAsubopt and extracts sequence, structures and energies.'
    def __init__(self, rnafold_output):
        output_lines = rnafold_output.strip().splitlines()
        # Top line of RNAsubopt has sequence and two numbers; split, take first.
        self.sequence = output_lines.pop(0).strip().split()[0]
        self.foldings = []
        for structure in output_lines:
            # Create a namedtuple of structure/energy by splitting line.
            structure, energy = structure.strip().split()
            self.foldings.append([structure, energy])

def RNAsubopt(sequence, *args):
    rnaf = subprocess.Popen(["RNAsubopt"] + list(args),
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            # Universal Newlines effectively allows string IO.
                            universal_newlines=True)        
    foldout, folderr = rnaf.communicate(sequence)
    if folderr:
        raise RNASuboptError(folderr)
    return RNASuboptOutput(foldout)
    
def RNAduplex(sequence, *args):
    duplex_out = subprocess.Popen(["RNAduplex"] + list(args),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                # Universal Newlines effectively allows string IO.
                                universal_newlines=True
                                )
    foldout, folderr = duplex_out.communicate(sequence)
    return foldout.splitlines()
  
  
def Kinfold(sequence, args):
    times_median = float(args[1])*math.log(2)
    Kinfold_out = subprocess.Popen(["Kinfold","--cut", args[0], "--time", args[1], "--num", args[2]],    
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                # Universal Newlines effectively allows string IO.
                                universal_newlines=True
                                )
    Kinout, Kinerr = Kinfold_out.communicate(sequence)
    trajectories = Kinout.splitlines()
    times = 0
    median_reached = 0
    for trajectory in trajectories:
        trajectory = trajectory.split()
        this_time = trajectory[2]
        #create running sum
        times = times + float(this_time)
        #to avoid excessive folding calculations, assume that distribution of foldings will assume an exponential distribution with a median of mean*ln(2). Then if at least half of the foldings meet this threshold, assume that mean would have been met.
        if float(this_time) >= times_median:
            median_reached = median_reached + 1        
    time_ave = times/float(args[2])
    
    #if the median was reached more than half of the time, then output the max as average folding time
    if median_reached >= float(args[2])/2:
        time_ave = float(args[1])
    return time_ave

def parse_dot(dot_seq):
    #parse each hairpin using dot output as input
    parsed = []
    current = 1
    backwaslast = 0
    
    for (position,term) in enumerate(dot_seq):        
        if term is ".":
            parsed.append(0)
        
        if term is "(":
            parsed.append("xx")
            if backwaslast == 1:
                current = current + 1
                backwaslast = 0
            
        if term is ")":
            parsed.append(current)
            #search parsed for forward bracket locations. Output most 3' value
            parsed_thusfar=[i for i,x in enumerate(parsed) if x=="xx"]
            last_for=max(parsed_thusfar)
            #convert identified forward bracket location to current hairpin identifier
            parsed[last_for]=current
            
            #keep track of when we switch into another hairpin
            backwaslast = 1
    return parsed

def RNAplot(sequence, *args):
    duplex_out = subprocess.Popen(["RNAplot", "--filename-full"] + list(args),
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                # Universal Newlines effectively allows string IO.
                                universal_newlines=True
                                )
    foldout, folderr = duplex_out.communicate(sequence)
    
def RNAcofold(sequence, *args):
    RNA_cofoldout = subprocess.Popen(["RNAcofold", "-C", '-p', "--filename-full"] + list(args),
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            # Universal Newlines effectively allows string IO.
                            universal_newlines=True
                            )
    foldout, folderr = RNA_cofoldout.communicate(sequence)