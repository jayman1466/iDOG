# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 17:29:22 2016

@author: Jaymin Patel
"""

from RNAfold_wrapper import RNAfold, RNAduplex, Kinfold, parse_dot
#from rRNA_bindingsite import rRNA_bindingsite

from Bio.Seq import Seq
import math
import random
import re


def dG_tot_calc(pre_sequence,orf_sequence,rRNA_sequence,kinetic):
    #parameters that are set in this script: limit of rRNA binding is kept to 3-12 residues upstream of start codon. MFE +/- 1.5 is queried for folding. first 15 nucleotides of orf are unfolded by the ribosome. Standby site is 4bases upstream of SD. For kinetics, 6 bases of RBS must be in a hairpin and dG_complex must be < 0. A warning is triggered if hairpin folding time in Kinfold AU is >200    
    
    
    #Find energy of initial mRNA
    full_mRNA = pre_sequence+orf_sequence
    initial_mRNA = RNAfold(full_mRNA).folding
    dG_initial = initial_mRNA[1]
    
    #Bind mRNA to rRNA, assuming mRNA is completely unfolded and output all possible structures within 1.5kcal/mol within the mfe. In doing so, only bind the 3' end of the mRNA with a length = rRNA_sequence+12, since spacing can't be greater than 12 
    
    #Truncate mRNA presequence if needed to limit the range of binding to 3 to 13 bases upstream of start codon      
    pre_sequence_truncated_size=13+len(rRNA_sequence)
    if pre_sequence_truncated_size < len(pre_sequence):
        pre_sequence_truncated=pre_sequence[-pre_sequence_truncated_size:-2]
        truncated_amount=len(pre_sequence)-pre_sequence_truncated_size
    else:
        pre_sequence_truncated=pre_sequence[0:-2]
        truncated_amount=0
    
    #run RNAduplex to bind rRNA to mRNA     
    folding_input='{}\n{}'.format(pre_sequence_truncated,rRNA_sequence)   
    structures=RNAduplex(folding_input," -e 1.5")
    
    
    lowest_energy=[100,100,100,100,100,"",100,100,0,[0,100,0,0]]
    
    for structure in structures:
        
        structure = structure.split()
        
        #correction for weird bug. Don't know why it works, but it does
        if structure[4] == "(":
            structure[4] = structure [5]
        #From each outputted structure, determine the rRNA binding site
        
        #find the region of the rRNA that bound and the bases overhanging at  both ends
        mRNA_bind_index = structure[1].split(',')
        rRNA_bind_index = structure[3].split(',')        

        rRNA_length = len(rRNA_sequence)

        dot_struct = structure[0].split('&')
        dot_struct_mRNA = dot_struct[0]
        dot_struct_rRNA = dot_struct[1]        

        #determine how many bases overhang from each side of the bound rRNA
        if dot_struct_rRNA[0] == '.':
            rRNA_5overhang = int(rRNA_bind_index[0])
        else:                
            rRNA_5overhang = int(rRNA_bind_index[0]) - 1        

        if dot_struct_rRNA[-1] == '.':
            rRNA_3overhang = rRNA_length - int(rRNA_bind_index[1]) + 1
        else:                
            rRNA_3overhang = rRNA_length - int(rRNA_bind_index[1])

        #determine 5' coordinate of rRNA binding site on the mRNA             
        if dot_struct_mRNA[0] == '.':
            mRNA_5coord = int(mRNA_bind_index[0]) + truncated_amount - rRNA_3overhang 
        else:                
            mRNA_5coord = int(mRNA_bind_index[0]) + truncated_amount - rRNA_3overhang - 1        
     
        #determine the spacing 
        if dot_struct_mRNA[-1] == '.':
            s_dist = pre_sequence_truncated_size - int(mRNA_bind_index[1]) - rRNA_5overhang + 1
        else:
            s_dist = pre_sequence_truncated_size - int(mRNA_bind_index[1]) - rRNA_5overhang
        
        #make sure the rRNA binding site is reasonable (spacing 3-13 bp from start codon)    
        if 3<=s_dist<=13:
            #Fold rRNA bound to mRNA as two different bindings
            #print(rRNA_bound)
        
            #fold the first half of the mRNA: everying thing before the standby site in front of rRNA binding site
            if mRNA_5coord > 4:    
                pt1=float(RNAfold(pre_sequence[0:mRNA_5coord-4]).folding[1])       
            else:
                pt1=0
                
            #the binding energy between SD and rRNA
            pt2=float(structure[4].strip('()'))

            #calculate the folding energy of the 15th-35th nucleotide of the orf
            if len(orf_sequence) > 16:
                pt3 = float(RNAfold(orf_sequence[16:]).folding[1])
            else:
                pt3 = 0

            #calculate the delta Gs
            s_penalty={2:9.63 , 3:1.53 , 4:0.01 , 5:0.00 , 6:0.29 , 7:0.67 , 8:1.15 , 9:1.73 , 10:2.40 , 11:3.17 , 12:4.03 , 13:4.99 , 14:6.05 , 15:7.20}
            start_binding={"AUG":-1.194 , "GUG":-0.0748 , "UUG":-0.0435, "CUG":-0.03406, "ATG":-1.194 , "GTG":-0.0748 , "TTG":-0.0435, "CTG":-0.03406}
            
            dG_mRNA_rRNA = pt1+pt2+pt3
            dG_spacing = s_penalty[s_dist]
            dG_start = start_binding[orf_sequence[0:3]]
            dG_final = dG_mRNA_rRNA + dG_spacing + dG_start
            dG_complex = pt2 + dG_spacing + dG_start 
            
            hairpin_trapped = 0
            slow_folding = 0
            folding_time = 0              
            #check if this structure results in a lower dG_final than the previous structure            
            if dG_final < lowest_energy[4]:
              
                if kinetic == 1:

                    #Check if dG_complex is strong enough for this to be an issue
                    if dG_complex < -3.0:
                        #if testing for kinetic issues, check if the RBS is situated within a hairpin
                        initial_folding = initial_mRNA[0]
                        #relevant portion is from the 5' end of the SD to the end of the start codon
                        initial_folding_relevant = initial_folding[mRNA_5coord:len(pre_sequence)+3]
                        #count the number of residues in a hairpin
                        hairpin_residues = initial_folding_relevant.count('(') + initial_folding_relevant.count(')')
                        #set threshold for number of residues that have to be sequestered in a hairpin
                        if hairpin_residues >= 6:
                            hairpin_trapped = 1                       

                            #parse the dot output of the mRNA
                            parsed = parse_dot(initial_folding)
                            
                            #determine which hairpins are in the relevant region. If there are multiple, choose the one with the most residues in the relevant region
                            relevant_parsed = [parsed[index] for index in range(mRNA_5coord,len(pre_sequence)+3)]
                            
                            relevant_hairpin = 1
                            residues_in_hairpin = 0                            
                            for hairpin in set(relevant_parsed):
                                if hairpin > 0:
                                    thiscount = relevant_parsed.count(hairpin)
                                    if thiscount >= residuemRNA_bind_indexs_in_hairpin:
                                        residues_in_hairpin = thiscount
                                        relevant_hairpin = hairpin
                            
                            #isolate the full sequence of the relevant hairpin 
                            indexes = [i for i,x in enumerate(parsed) if x==relevant_hairpin]
                            
                            first_index = min(indexes)
                            last_index = max(indexes)
                            
                            relevant_hairpin_sequence = full_mRNA[first_index:last_index]
                            
                            #get the MFE of the relevant hairpin
                            MFE = str(RNAfold(relevant_hairpin_sequence).folding[1])
    
                            #max time
                            time_max="250"
                            
                            #number of trajectories
                            tra_num="50"
                            
                            args = [MFE,time_max,tra_num]
                            
                            #use Kinfold to determine folding time                        
                            folding_time = Kinfold(relevant_hairpin_sequence,args)
                            
                            #check if this folding time is too long
                            if folding_time >= 250:
                                slow_folding = 1
                            
                #calculate translation initiation rate
                dG_tot=dG_final-dG_initial
                TIR=2500*math.exp(-0.45*dG_tot)
                lowest_energy=[dG_tot,dG_mRNA_rRNA,dG_spacing,dG_start,dG_final,dG_initial,s_dist,mRNA_5coord,TIR, [hairpin_trapped,dG_complex,folding_time,slow_folding]]
                #[dG_tot,dG_mRNA_rRNA,dG_spacing,dG_start,dG_final,dG_initial,spacing,length of mRNA 5' of rRNA binding site,Translation initiation rate,dG of initiation complex, assuming unfolded mRNA]
    
    return lowest_energy


    
def relevant_mRNA(mRNA_this,start_this): #outputs +/- 35 bp around the start codon
    if start_this>=35:
        pre_this=mRNA_this[start_this-35:start_this]
    else:
        pre_this=mRNA_this[0:start_this]
    if len(mRNA_this)-start_this>38:
        orf_this=mRNA_this[start_this:start_this+35]
    else:
        orf_this=mRNA_this[start_this:]
    return pre_this, orf_this





def int_RBS_remove(leader,orf,rRNA_sequence):
    #sequence parameters
    leader=leader
    orf=orf 
    ignore_first=0 #ignore start codons within this number of the first start codon (does not apply to in-frame GTG, CTG, and TTG)
    int_RBS_thresh=50 #TIR threshold for internal RBSs
    codon_change_dist=[-15]*1+[-12]*2+[-9]*3+[-6]*3+[-3]*3+[0]*1  #distribution of which codon will be recoded to remove internal RBS
    rRNA_sequence=rRNA_sequence #anti shine delgarno sequence
    attempts_max=20 #maximum number of attempts made at each start codon
    first_codon_thresh=0.8 #maximum allowable decrease in the first start codon's strength
    
    #dictionary for reverse translation - removed ttg, ctg, and gtg from list to prevent putting in excessive new start codons. and removed TTA for streptomyces
    A_heg=['GCA']+['GCC']+['GCG']+['GCU']
    C_heg=['UGC']+['UGU']
    D_heg=['GAC']+['GAU']
    E_heg=['GAA']+['GAG']
    F_heg=['UUC']+['UUU']
    G_heg=['GGA']+['GGC']+['GGG']+['GGU']
    H_heg=['CAC']+['CAU']
    I_heg=['AUA']+['AUC']+['AUU']
    K_heg=['AAA']+['AAG']
    L_heg=['CUA']+['CUC']+['CUU']
    M_heg=['AUG']
    N_heg=['AAC']+['AAU']
    P_heg=['CCA']+['CCC']+['CCG']+['CCU']
    Q_heg=['CAA']+['CAG']
    R_heg=['AGA']+['AGG']+['CGA']+['CGC']+['CGG']+['CGU']
    S_heg=['AGC']+['AGU']+['UCA']+['UCC']+['UCG']+['UCU']
    T_heg=['ACA']+['ACC']+['ACG']+['ACU']
    V_heg=['GUA']+['GUC']+['GUU']
    W_heg=['UGG']
    Y_heg=['UAC']+['UAU']    
    
    first_start_index=len(leader) #identify first start codon
    
    mRNA_sequence=(leader+orf).upper().replace("T","U")
        
    #determine initial dG_tot of first start codon 
    first_start_rel_mRNA=relevant_mRNA(mRNA_sequence,first_start_index)
    first_start_TIR=dG_tot_calc(first_start_rel_mRNA[0],first_start_rel_mRNA[1],rRNA_sequence,0)[8]
    
    #first, get rid of in frame GTG, TTG codons 
    start_indexes=[m.start()-1 for m in re.finditer('(?=UG)', mRNA_sequence)]
        
    for start_index in start_indexes:
        if start_index>first_start_index and (start_index-first_start_index)%3==0:
            if mRNA_sequence[start_index:start_index+3]=="GUG":
                mRNA_sequence=mRNA_sequence[0:start_index]+random.choice(V_heg)+mRNA_sequence[start_index+3:]
            
            '''
            if mRNA_sequence[start_index:start_index+3]=="CUG":
                mRNA_sequence=mRNA_sequence[0:start_index]+random.choice(L_heg)+mRNA_sequence[start_index+3:]
            '''
            
            if mRNA_sequence[start_index:start_index+3]=="UUG":
                mRNA_sequence=mRNA_sequence[0:start_index]+random.choice(L_heg)+mRNA_sequence[start_index+3:]
                
    #then, get rid of all other internal start codons
                
    counter=[0]*(len(mRNA_sequence)-1) #start a counter for number of attempts at each start_index
    last_codon_opt=0
    while True:
        i=0  
        #Find all start codons    
        start_indexes=[m.start()-1 for m in re.finditer('(?=UG)', mRNA_sequence)]
            
        #iterate through start codons to find internal start codons 
        for start_index in start_indexes:
            
            #check if start codon is in the targetted range and that it hasn't reached its max attempts and isn't more than 60bp upstream of the last fixed RBS
            if start_index>first_start_index+ignore_first and counter[start_index]<=attempts_max and start_index>last_codon_opt-60:
        
                rel_mRNA=relevant_mRNA(mRNA_sequence,start_index)
                this_dG_tot=dG_tot_calc(rel_mRNA[0],rel_mRNA[1],rRNA_sequence,0)
                 
                if this_dG_tot[8] > int_RBS_thresh:
                    #identify this as the last optimized codon                
                    last_codon_opt=start_index
                    #add 1 to the counter of this codon index
                    counter[start_index]=counter[start_index]+1
                        
                    #check if start codon is in frame of the orf
                    remainder=(start_index-first_start_index)%3
                    
                    if remainder == 0:
                        #find the spacing between rRNA binding site and start codon, and adjust for frame
                        space_adjusted=(this_dG_tot[6]-(this_dG_tot[6]%3))
                        
                        #recode an amino acid upstream, chosen by codon_change_dist 
                        rand_codon=random.choice(codon_change_dist)
                        recoding_target=Seq(mRNA_sequence[start_index-space_adjusted+rand_codon:start_index-space_adjusted+rand_codon+3])
                        recoding_target_tr=recoding_target.translate()
                        rec_out=""
                        for c in recoding_target_tr:
                            codon_dict={'A':random.choice(A_heg),'C':random.choice(C_heg),'D':random.choice(D_heg),'E':random.choice(E_heg),'F':random.choice(F_heg),'G':random.choice(G_heg),'H':random.choice(H_heg),'I':random.choice(I_heg),'K':random.choice(K_heg),'L':random.choice(L_heg),'M':random.choice(M_heg),'N':random.choice(N_heg),'P':random.choice(P_heg),'Q':random.choice(Q_heg),'R':random.choice(R_heg),'S':random.choice(S_heg),'T':random.choice(T_heg),'V':random.choice(V_heg),'W':random.choice(W_heg),'Y':random.choice(Y_heg),'*':'UAA'}
                            rec_out+=codon_dict[c]
                        
                        temp_mRNA_sequence=mRNA_sequence[0:start_index-space_adjusted+rand_codon]+rec_out+mRNA_sequence[start_index-space_adjusted+rand_codon+3:]
                        
                        #see if codon change improves the situation
                        temp_rel_mRNA=relevant_mRNA(temp_mRNA_sequence,start_index)                 
                        temp_this_dG_tot=dG_tot_calc(temp_rel_mRNA[0],temp_rel_mRNA[1],rRNA_sequence,0)
                        
                        if temp_this_dG_tot[0] >= this_dG_tot[0]:
                        
                            #make sure this change does not mess up the initial start codon, if in range of initial start codon
                            if start_index-rand_codon < first_start_index+70:
                                first_rel_mRNA=relevant_mRNA(temp_mRNA_sequence,first_start_index)
                                first_this_dG_tot=dG_tot_calc(first_rel_mRNA[0],first_rel_mRNA[1],rRNA_sequence,0)                            
                                if first_this_dG_tot[8] > first_start_TIR*first_codon_thresh:
                                    mRNA_sequence=temp_mRNA_sequence
                            else:
                                mRNA_sequence=temp_mRNA_sequence
                    else:
                        #try to just recode the start codon itself if start codon is out of frame
                        recoding_target=Seq(mRNA_sequence[start_index-remainder:start_index-remainder+6])
                        recoding_target_tr=recoding_target.translate()
                        rec_out=""
                        for c in recoding_target_tr:
                            codon_dict={'A':random.choice(A_heg),'C':random.choice(C_heg),'D':random.choice(D_heg),'E':random.choice(E_heg),'F':random.choice(F_heg),'G':random.choice(G_heg),'H':random.choice(H_heg),'I':random.choice(I_heg),'K':random.choice(K_heg),'L':random.choice(L_heg),'M':random.choice(M_heg),'N':random.choice(N_heg),'P':random.choice(P_heg),'Q':random.choice(Q_heg),'R':random.choice(R_heg),'S':random.choice(S_heg),'T':random.choice(T_heg),'V':random.choice(V_heg),'W':random.choice(W_heg),'Y':random.choice(Y_heg),'*':'UAA'}
                            rec_out+=codon_dict[c]
                            
                        mRNA_sequence=mRNA_sequence[0:start_index-remainder]+rec_out+mRNA_sequence[start_index-remainder+6:]
                    
                    break           
            i=i+1
        if i==len(start_indexes):
            break
    #return updated RBS and ORF sequence 
    return(mRNA_sequence[0:len(leader)],mRNA_sequence[len(leader):])
    
    
'''
test = dG_tot_calc('TCTCTGCTCCTGATGAGGGCGCTAA','ATGCCGCGTGGCCTGGAATTATTGATTGCTCAAAC','ACCUCCUUA',0)
print(test)
'''