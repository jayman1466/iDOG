# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 14:15:15 2017

@author: root
"""
import subprocess
import tempfile

def transterm(seq):

    #os.remove('ttermseq.fasta')
    
    #write the sequence temporary file
    termseq = tempfile.NamedTemporaryFile(suffix='.fa')
    termseq.write(b">seq\n")
    termseq.write(seq.encode('utf-8'))
    termseq.seek(0)

    #write the annotations temporary file
    annotation = tempfile.NamedTemporaryFile(suffix='.coords')
    size = len(seq)
    line = f"fakegene1\t1\t2\tseq\nfakegene2\t{str(size)}\t{str(size-1)}\tseq"
    annotation.write(line.encode('utf-8'))
    annotation.seek(0)
    
    
    #run the command line
    transterm_out = subprocess.check_output(['./transterm', '-p', 'transtermhp/expterm.dat', termseq.name, annotation.name])
    termseq.close()
    annotation.close()
    
    #decode output
    transterm_out = transterm_out.decode("utf-8")
    terminators = transterm_out[transterm_out.find("SEQUENCE seq"):]
    outlen = (len(terminators))
    if outlen > 125:
        return "terminators alert"
    else:
        return "clear of terminators"
'''    
a = transterm("AtggcaatgtttaaaagagaagaaatcattgaaatggccaataaggactttgaaaaagcatggatcgaaactaaagaccttataaaagctaaaaagataaacgaaagttacccaagaataaaaccagtttttggaaaaacacaccctgtaaatgacactattgaaaatttaagacaggcatatcttagaatgggttttgaagaatatataaacccagtaattgtcgatgaaagagatatttataaacaattcggcccagaagctatggcagttttggatagatgcttttatttagcgggacttccaagacctgacgttggtttgagcgatgaaaaaatttcacagattgaaaaacttggaattaaagtttctgagcacaaagaaagtttacaaaaaatacttcacggatacaaaaaaggaactcttgatggtgacgatttagttttagaaatttcaaatgcacttgaaatttcaagcgagatgggtttaaaaattttagaagatgttttcccagaatttaaggatttaaccgcagtttcttcaaaattaactttaagaagccacatgacttcaggatggttccttactgtttcagacctcatgaacaaaaaacccttgccatttaaactcttttcaatcgatagatgttttagaagagaacaaaaagaagataaaagccacttaatgacataccactctgcatcctgtgcaattgcaggtgaaggcgtggatattaatgatggaaaagcaattgcagaaggattattatcccaatttggctttacaaactttaaattcattcctgatgaaaagaaaagtaaatactacacccctgaaacacagactgaagtttacgcataccacccaaaattaaaagaatggctcgaagttgctacatttggagtatattcgccagttgcattaagcaaatacggaatagatGTACCTGTAATGAATTTGgcgtatGGTGTTGAAAGACTTGCaatgatttctggaaatttcgcagatgttcgagaaatggtatatcctcagttttacgaacacaaacttaatgaccggaatgtcgcttcaatggtaaaactcgataaagttccagtaatggatgaaatttacgatttaacaaaagaattaattgagtcatgtgttaaaaacaaagatttaaaatccccttgtgaattagctattgaaaaaacgttttcatttggaaaaaccaagaaaaatgtaaaaataaacatttttCCGAAATTCGAAGGTAAAAATTTACTCGGACCTTCAATTTTAAACGAAATCTACGTTTACGATGGAAATGTAATTGGAATTCCTGAAAGCTTTGACGGAGTAAAAGAAGAATTTAAAGACTTCTTAGAAAAAGGAAAATCAGAAGGGGTAGCAACAGGCATTCGATATATCGATGCGCTTTGCTTTAAAATTACTTCAAAATTAGAAGAAGCATTTGTGTCAAACACTACTGAATTCAAAGTTAAAGTTATGTGGGTCAGAAGTTTAAGCGACATTAACTTAAAAATCGATGATATCGCATTAAAACAGATCATGAGCAAAAATAAAGTAATCGACGTTAGAGGCCCAGTCAGCTTAAATGTCGAAGTAAAAATTGAATAAgagCCcggtacctaaagtatattagttaagtataagaaggagatatacatatgggatccatgtctaaagaaaagtttgaacgtacaaaaccgcacgttAACGTCGGTACTATCGGCCACGTTGACCATGGTAAAACAACGCTGACCGCTGCAATCACTACCGTACTGGCTAAAACCTACGGCGGTGCTGCTCGCGCATTCGACCAGATCGATAACGCGCCGGAAGAAAAAGCTCGTGGTATCACCATCAACACTTCTCGGGTTGAATACGACACCCCGACCCGTCACTACGCACACGTAGACTGCCCGGGGCACGCCGACTATGTTAAAAACATGATCACCGGTGCTGCGCAGATGGACGGCGCGATCCTGGTAGTTGCTGCGACTGACGGCCCGATGCCGCAGACTCGTGAGCACATCCTGCTGGGTCGTCAGGTAGGCGTTCCGTACATCATCGTGTTCCTGAACAAATGCGACATGGTTGATGACGAAGAGCTGCTGGAACTGGTTGAAATGGAAGTTCGTGAACTTCTGTCTCAGTACGACTTCCCGGGCGACGACACTCCGATCGTTCGTGGTTCTGCTCTGAAAGCGCTGGAAGGCGACGCAGAGTGGGAAGCGAAAATCCTGGAACTGGCTGGCTTCCTGGATTCTTACATTCCGGAACCAGAGCGTGCGATTGACAAGCCGTTCCTGCTGCCGATCgtCGGGGTATACTCCATCTCCGGTCGTGGTACCGTTGTTTCGGGTCGTGTAGAACGCGGTATCATCAAAGTTGGTGAAGAAGTTGAAATCGTTGGTATCAAAGAGACTCAGAAGTCTACCTGTACTGGCGTTGAAATGTTCCGCAAACTGCTGGACGAAGGCCGTGCTGGTGAGTGGGTAGGTGTTCTGCTGCGTGGTATCAAACGTGAAGAAATCGAACGTGGTCAGGTACTGGCTAAGCCGGGCACCATCAAGCCGCACACCAAGTTCGAATCTGAAGTGTACATTCTGTCCAAAGATGAAGGCGGCCGTCATACTCCGTTCTTCAAAGGCTACCGTCCGCAGTTCTACTTCCGTACTACTGACGTGACTGGTACCATCGAACTGCCGGAAGGCGTAGAGATGGTAATGCCGGGCGACAACATCAAAATGGTTGTTACCCTGATCCACCCGATCGCGATGGACGACGGTCTGCGTTTCGCAATCCGTGAAGGCGGCCGTACCGTTGGCGCGGGCGTTGTAGCAAAAGTTCTGAGCTAATGGGCTAACAGGAGGAAttaCAT")
'''