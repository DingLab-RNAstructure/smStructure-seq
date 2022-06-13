"""

 We parse the m5 formatted File and read the Ref length from the fasta file
 - we check if the bases match wild type  0
 - if observed base 


singularity exec forgi-centroid-mea-final.simg  python3 m5_to_bitvectors.py  --input_file  R1_5p.m5  --transcript COOLAIR3  --reference_file  reference/cool6.fasta  --output_file  R1_5p.bit

# md5sum: 5f3f2812fd3a9e0bfe44fcbfa8e4d4b0  R1_5p.bit

"""

import sys  
import os
import argparse
import jitu
from dotmap import DotMap

def handler():
    """
    Get command line inputs
    """ 
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference_file",   type=str, required=True)
    parser.add_argument("-t", "--transcript",   type=str, required=True)
    parser.add_argument("-i", "--input_file",   type=str, required=True)
    parser.add_argument("-o", "--output_file",  type=str, required=True)

    args = parser.parse_args()

    return args


def get_state_vector(Ref, d):
    """
    Find the obserevd mutation and set the bit on    

    """  

    bitD = { i: '.' for i in range(len(Ref))} # initilize NA ==> '.'


    beg = int(d.tStart)
    BLOCK = list(zip(d.tAlignedSeq, d.matchPattern, d.qAlignedSeq))

    if d.tStrand == '-': # only flip 
       BLOCK=BLOCK[::-1]  # whole orientation is flipped

    for t, tick, q in BLOCK: # we could only iterate and check tick == '|'
        if t != '-': # skip the indel in Reference position 
           if t == q: 
              bitD[beg] = '0'
           else:
              bitD[beg] = '1'  # either base difference or indel in read

           beg += 1 # skipping

    state = '\t'.join([bitD[i] for i in range(len(Ref))])

    return state 

if __name__ == "__main__":


   args = handler()
   if not os.path.exists(args.input_file):
      print(f'Input file not found: {args.input_file}')
      exit()

   if not os.path.exists(args.reference_file):
      print(f'Reference file not found: {args.reference_file}')
      exit()

   # inputs 
   ref_path  = args.reference_file  # 'reference/cool6.fasta'
   m5file    =  args.input_file      # 'R1_5p.m5'
   out_bitFile = args.output_file  # 'R1_5p.bit'

   # load ref...
   tube, seqD  = jitu.getTubeD(ref_path)

   assert args.transcript in tube, "check the transcript name in reference file: " +  args.transcript

   head = ['qName', 'qLength', 'qStart', 'qEnd', 'qStrand', 'tName', 'tLength', 'tStart', 'tEnd', 'tStrand',  'score', 'numMatch', 'numMismatch', 'numIns', 'numDel', 'mapQV', 'qAlignedSeq', 'matchPattern', 'tAlignedSeq']
   
  
   with open(m5file) as inp, open(out_bitFile, 'w') as outf: 

        for line in inp:
          A = line.strip().split()
          d = DotMap( dict(zip(head,A) ))
            
          this_transcript = d.tName

          if args.transcript != this_transcript:
             continue
 
          Ref = seqD[this_transcript] # COOLAIR3 : original ref seq

          state  = get_state_vector(Ref, d)
          state_na = state.replace('.', 'NA')  
          outf.write('\t'.join([d.qName] + state_na.split('\t') ) + '\n')


print ('DONE')
    
"""
COOLAIR4
CGCGCAGAGAGAGAGAGAGAGCCACGTCCCTGTTGCAAAATAAGCCGTAGGCTTCTTCACTGTGAAGCAAACACAAGTTTTTGACAGAAGTGAAGAACACATACACTCAAGATCTCGATGCAATTCTCACACGAATAAGAAAAGTAAAAGAGCACAAAACAGAAGATAAAAGGGGGAACAAATGAAAACCCAGGTAAGGAAAAGGCGTACTTATCGCCGGAGGAGAAGCTGTAGAGCTTGCCGGAGGCGGAGACGACGAGAAGAGCGACGGATGCGTCACAGAGAACAGAAAGCTGACGAGCTTTCTCGATGAGACCGTTGCGACGTTTGGAGAAGGTGACTTGTCGGCTACTTTTGTTCTCAATTCGCTTGATTTCTAGTTTTTTTCTTCCCATGGCTTCTCTCCGAGAGGGCTTTGTGCCCTAATTTGATCCTCAGGTTTGGGTTCAAGTCGCCGGAGATACTAAGCGTTTTCTCTTTCTATTTTTTTTTTTCCTTTTCTCGCTTTATTTCTTTCTATTTTTTGTGCCTATCTACTTTTTCTTCGTCGGGCCAGATATTTTTTTTTTTTTGGGGGTAAACGAGAGTGATGCAAAAAAACCAAATATGTGAATAAAAACGTTGTGTTTTGAAGACAAGATTGCCACGTGTACCGC
                                    AAAATAAGCCGTAGGCTTCTTCACTGTGAAGCAAACACAAGTTTTTGACAGAAGTGAAGAACACATACACTCAAGATCTCGATGCAATTCTCACACGAATAAGGTGGCTAATTAAGTAGTGGGAGAGT--CAC-CGGAAGATT-GTCGGAGATTTGTCC-AGCAGGTGACA-TCTCCA-TC-TCAGCTTCTG-CTCCCACATGATGATTATTCTCCATCTGTACGATAATC-A-T-A-GAAAAGTAAAAGAGCACAAAACAGAAGATAAAAGGGGGAACAAATGAAAACCCAGGTAAGGAAAAGGCGTACTTATCGCCGGAGGAGAAGCTGTAGAGCTTGCCGGAGGCGGAGACGACGAGAAGAGCGACGGATGCGTCACAGAGAACAGAAAGCTGACGAGCTTTCTCGATGAGACCGTTGCGACGTTTGGAGAAGGTGACTTGTCGGCTACTTTTGTTCTCAATTCGCTTGATTTCTAGTTTTTTTCTTCCCATGGCTTCTCTCCGAGAGGGCTTTGTGCCCTAATTTGATCCTCAGGTTTGGGTTCAAGTCGCCGGAGATACTAAGCGTTTTCTCTTTCTATTTTTTTTTTTCCTTTTCTCGCTTTATTTCTTTCTATTTTTTGTGCCTATCTACTTTTTCTTCGTCGGGCCAGATATTTTTTTTTTTTT-GGGGGTAAACGAGAGTGATGCAAAAAAACCAAATATGTGAATAAAAACGTTGTGTTTTGAAGACAA
656
DotMap(mapQV='0', numDel='34', numMatch='667', qEnd='709', score='-2951', 


qName='m54270_210831_092625/5111880/ccs/0_709', qStrand='+', tEnd='761', 

tStrand='-', qLength='709', tStart='36', tLength='779', tName='COOLAIR4', numIns='14', 

tAlignedSeq  = 'AAAATAAGCCGTAGGCTTCTTCACTGTGAAGCAAACACAAGTTTTTGACAGAAGTGAAGAACACATACACTCAAGATCTCGATGCAATTCTCACACGAATAAGGTGGCTAATTAAGTAGTGGGAGAGT--CAC-CGGAAGATT-GTCGGAGATTTGTCC-AGCAGGTGACA-TCTCCA-TC-TCAGCTTCTG-CTCCCACATGATGATTATTCTCCATCTGTACGATAATC-A-T-A-GAAAAGTAAAAGAGCACAAAACAGAAGATAAAAGGGGGAACAAATGAAAACCCAGGTAAGGAAAAGGCGTACTTATCGCCGGAGGAGAAGCTGTAGAGCTTGCCGGAGGCGGAGACGACGAGAAGAGCGACGGATGCGTCACAGAGAACAGAAAGCTGACGAGCTTTCTCGATGAGACCGTTGCGACGTTTGGAGAAGGTGACTTGTCGGCTACTTTTGTTCTCAATTCGCTTGATTTCTAGTTTTTTTCTTCCCATGGCTTCTCTCCGAGAGGGCTTTGTGCCCTAATTTGATCCTCAGGTTTGGGTTCAAGTCGCCGGAGATACTAAGCGTTTTCTCTTTCTATTTTTTTTTTTCCTTTTCTCGCTTTATTTCTTTCTATTTTTTGTGCCTATCTACTTTTTCTTCGTCGGGCCAGATATTTTTTTTTTTTT-GGGGGTAAACGAGAGTGATGCAAAAAAACCAAATATGTGAATAAAAACGTTGTGTTTTGAAGACAA', 
matchPattern = '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*|*|*|*|*||*|*|*||***|**||**||||||**|******||*|||**|||*|**|*|**|*|*|*|*|*|||*|**|*|***|*|**|||***||||**|*|*|**|*|*|****|*|*||**||**|*|*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 
qAlignedSeq  = 'TTGTCTTCAAAACACAACGTTTTTATTCACATATTTGGTTTTTTTGCATCACTCTCGTTTACCCCCAAAAAAAAAAAAAATATCTGGCCCGACGAAGAAAAAGTAGATAGGCACAAAAAATAGAAAGAAATAAAGCGAGAAAAGG-AAAAAAAAAATAGAAAGAGAAAACGCTTAGTATCTCCGGCGACTTGAACCCAAACCTGAGGATCAAATTAGGGCACAAAGCCCTCTCGGAGAGAAGCCATGGGAAGAAAAAAACTAGAAATCAAGCGAATTGAGAACAAAAGTAGCCGACAAGTCACCTTCTCCAAACGTCGCAACGGTCTCATCGAGAAAGCTCGTCAGCTTTCTGTTCTCTGTGACGCATCCGTCGCTCTTCTCGTCGTCTCCGCCTCCGGCAAGCTCTACAGCTTCTCCTCCGGCGATAAGTACGCCTTTTCCTTACCTGGGTTTTCATTTGTTCCCCCTTTTATCTTCTGTTTTGTGCTCTTTTACTTTTCTTAACTAG-TT-T-G-AC--TTTAAG-TTAATCA--A------AGCCAG-CGCT-ATCACT--AAACT-TTATCTG-TATG-CCTTTGTATGACTTTTCTT-TGAGGGAAAATGT----C-A-TT--TT--CAATCTTATTCGTGTGAGAATTGCATCGAGATCTTGAGTGTATGTGTTCTTCACTTCTGTCAAAAACTTGTGTTTGCTTCACAGTGAAGAAGCCTACGGCTTATTTT', 

qStart='4', numMismatch='24')




>18S
TACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTGTAAGTATGAACGAATTCAGACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAACCCCGACTTATGGAAGGGACGCATTTATTAGATAAAAGGTCGACGCGGGCTCTGCCCGTTGCTCTGATGATTCATGATAACTCGACGGATCGCATGGCCTCTGTGCTGGCGACGCATCATTCAAATTTCTGCCCTATCAACTTTCGATGGTAGGATAGTGGCCTACCATGGTGGTAACGGGTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACACGGGGAGGTAGTGACAATAAATAACAATACTGGGCTCTTTCGAGTCTGGTAATTGGAATGAGTACAATCTAAATCCCTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAACCTTGGGATGGGTCGGCCGGTCCGCCTTTGGTGTGCATTGGTCGGCTTGTCCCTTCGGTCGGCGATACGCTCCTGGTCTTAATTGGCCGGGTCGTGCCTCCGGCGCTGTTACTTTGAAGAAATTAGAGTGCTCAAAGCAAGCCTACGCTCTGGATACATTAGCATGGGATAACATCATAGGATTTCGATCCTATTGTGTTGGCCTTCGGGATCGGAGTAATGATTAACAGGGACAGTCGGGGGCATTCGTATTTCATAGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACAACTGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGAAAGTTGGGGGCTCGAAGACGATCAGATACCGTCCTAGTCTCAACCATAAACGATGCCGACCAGGGATCAGCGGATGTTGCTTATAGGACTCCGCTGGCACCTTATGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTTACCAGGTCCAGACATAGTAAGGATTGACAGACTGAGAGCTCTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGTTAATGAACGAGACCTCAGCCTGCTAACTAGCTACGTGGAGGCATCCCTTCACGGCCGGCTTCTTAGAGGGACTATGGCCGTTTAGGCCAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGTATTCAACGAGTTCACACCTTGGCCGACAGGCCCGGGTAATCTTTGAAATTTCATCGTGATGGGGATAGATCATTGCAATTGTTGGTCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGCGTTGACTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAATGATCCGGTGAAGTGTTCGGATCGCGGCGACGTGGGTGGTTCGCCGCCCGCGACGTCGCGAGAAGTCCACTAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTG
...........................................................................0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000........
............................................................................0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.............................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................        m54270_210831_092625/4915585/ccs/0_113
                                                                           TTGTTGCACGTATTAGCTCTAGAATTACTACGGTTATCCGAGTAGTAGTTACCATCAAACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT
                                                                           TTGTTGCACGTATTAGCTCTAGAATTACTACGGTTATCCGAGTAGTAGTTACCATCAAACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTct 

m54270_210831_092625/4915585/ccs/0_113 113 9 113 +  18S 1808 76 179 - -510 103 0 1 0 254 TTGTTGCACGTATTAGCTCTAGAATTTACTACGGTTATCCGAGTAGTAGTTACCATCAAACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT ||||||||||||||||||||||||*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| TTGTTGCACGTATTAGCTCTAGAA-TTACTACGGTTATCCGAGTAGTAGTTACCATCAAACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT

m54270_210831_092625/4915585/ccs/0_113 113 9 113 +  18S 1808 76 179 - -510 103 0 1 0 254 

TTGTTGCACGTATTAGCTCTAGAATTTACTACGGTTATCCGAGTAGTAGTTACCATCAAACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT 
||||||||||||||||||||||||*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
TTGTTGCACGTATTAGCTCTAGAA-TTACTACGGTTATCCGAGTAGTAGTTACCATCAAACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT

## answer


# good read  Easy  + + Approved 
m54270_210831_092625/5439813/ccs/0_108 108 8 107 +  18S 1808 84 186 + -480 99 0 0 3 254 AAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCT-ATACGTGCAAC-AA-CCCGA |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*|||||||||||*||*||||| AAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAACCCCGA
>18S
TACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTGTAAGTATGAACGAATTCAGACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAACCCCGACTTATGGAAGGGACGCATTTATTAGATAAAAGGTCGACGCGGGCTCTGCCCGTTGCTCTGATGATTCATGATAACTCGACGGATCGCATGGCCTCTGTGCTGGCGACGCATCATTCAAATTTCTGCCCTATCAACTTTCGATGGTAGGATAGTGGCCTACCATGGTGGTAACGGGTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACACGGGGAGGTAGTGACAATAAATAACAATACTGGGCTCTTTCGAGTCTGGTAATTGGAATGAGTACAATCTAAATCCCTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAACCTTGGGATGGGTCGGCCGGTCCGCCTTTGGTGTGCATTGGTCGGCTTGTCCCTTCGGTCGGCGATACGCTCCTGGTCTTAATTGGCCGGGTCGTGCCTCCGGCGCTGTTACTTTGAAGAAATTAGAGTGCTCAAAGCAAGCCTACGCTCTGGATACATTAGCATGGGATAACATCATAGGATTTCGATCCTATTGTGTTGGCCTTCGGGATCGGAGTAATGATTAACAGGGACAGTCGGGGGCATTCGTATTTCATAGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACAACTGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGAAAGTTGGGGGCTCGAAGACGATCAGATACCGTCCTAGTCTCAACCATAAACGATGCCGACCAGGGATCAGCGGATGTTGCTTATAGGACTCCGCTGGCACCTTATGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTTACCAGGTCCAGACATAGTAAGGATTGACAGACTGAGAGCTCTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGTTAATGAACGAGACCTCAGCCTGCTAACTAGCTACGTGGAGGCATCCCTTCACGGCCGGCTTCTTAGAGGGACTATGGCCGTTTAGGCCAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGTATTCAACGAGTTCACACCTTGGCCGACAGGCCCGGGTAATCTTTGAAATTTCATCGTGATGGGGATAGATCATTGCAATTGTTGGTCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGCGTTGACTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAATGATCCGGTGAAGTGTTCGGATCGCGGCGACGTGGGTGGTTCGCCGCCCGCGACGTCGCGAGAAGTCCACTAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTG
                                                                                    AAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAACCCCGA
....................................................................................000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000100100000......................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................        m54270_210831_092625/5439813/ccs/0_108

                                                                                    AAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCT-ATACGTGCAAC-AA-CCCga read
                                                                                    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*|||||||||||*||*||||| 
                                                                                    AAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAACCCCGA ref


AAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCT-ATACGTGCAAC-AA-CCCga read
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*|||||||||||*||*||||| 
AAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAACCCCGA ref

....................................................................................000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000100100000......................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................        m54270_210831_092625/5439813/ccs/0_108

# hard read m54270_210831_092625/7013299/ccs/0_109
>18S
TACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTGTAAGTATGAACGAATTCAGACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAACCCCGACTTATGGAAGGGACGCATTTATTAGATAAAAGGTCGACGCGGGCTCTGCCCGTTGCTCTGATGATTCATGATAACTCGACGGATCGCATGGCCTCTGTGCTGGCGACGCATCATTCAAATTTCTGCCCTATCAACTTTCGATGGTAGGATAGTGGCCTACCATGGTGGTAACGGGTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACACGGGGAGGTAGTGACAATAAATAACAATACTGGGCTCTTTCGAGTCTGGTAATTGGAATGAGTACAATCTAAATCCCTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAACCTTGGGATGGGTCGGCCGGTCCGCCTTTGGTGTGCATTGGTCGGCTTGTCCCTTCGGTCGGCGATACGCTCCTGGTCTTAATTGGCCGGGTCGTGCCTCCGGCGCTGTTACTTTGAAGAAATTAGAGTGCTCAAAGCAAGCCTACGCTCTGGATACATTAGCATGGGATAACATCATAGGATTTCGATCCTATTGTGTTGGCCTTCGGGATCGGAGTAATGATTAACAGGGACAGTCGGGGGCATTCGTATTTCATAGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACAACTGCGAAAGCATTTGCCAAGGATGTTTTCATTAATCAAGAACGAAAGTTGGGGGCTCGAAGACGATCAGATACCGTCCTAGTCTCAACCATAAACGATGCCGACCAGGGATCAGCGGATGTTGCTTATAGGACTCCGCTGGCACCTTATGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTTACCAGGTCCAGACATAGTAAGGATTGACAGACTGAGAGCTCTTTCTTGATTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGTTAATGAACGAGACCTCAGCCTGCTAACTAGCTACGTGGAGGCATCCCTTCACGGCCGGCTTCTTAGAGGGACTATGGCCGTTTAGGCCAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGTATTCAACGAGTTCACACCTTGGCCGACAGGCCCGGGTAATCTTTGAAATTTCATCGTGATGGGGATAGATCATTGCAATTGTTGGTCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGCGTTGACTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAATGATCCGGTGAAGTGTTCGGATCGCGGCGACGTGGGTGGTTCGCCGCCCGCGACGTCGCGAGAAGTCCACTAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTG

mapping                                                                     TTGTTGCACGTATTAGCTCTAGAATTACTACGGTTATCC--G-AGTAGTTACCATC-AACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT 
                                                                            |||||||||||||||||||||||||||||||||||||||**|*|||||||||||||*|||||||||||||||||||||||||||||||||||||||||||||| 
                                                                            TTGTTGCACGTATTAGCTCTAGAATTACTACGGTTATCCGAGTAGTAGTTACCATCAAACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT

flipped                                                                     TCTGACACTTTGACGCTTACCGAGTAATTTAGTCAATATCAAACAAACTACCATTGATGATGAGCCTATTGGCATCATTAAGATCTCGATTATGCACGTTGTT                                                                            
                                                                            ||||||||||||||||||||||||||||||||||||||||||||||*|||||||||||||*|**|||||||||||||||||||||||||||||||||||||||                                                                            
                                                                            TCTGACACTTTGACGCTTACCGAGTAATTTAGTCAATATCAAACAA-CTACCATTGATGA-G--CCTATTGGCATCATTAAGATCTCGATTATGCACGTTGTT


............................................................................0000000000000000000000000000000000000000000000100000000000001011000000000000000000000000000000000000000.............................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................        m54270_210831_092625/7013299/ccs/0_109





m54270_210831_092625/7013299/ccs/0_109 109 10 109 +  18S 1808 76 179 - -475 99 0 0 4 254 TTGTTGCACGTATTAGCTCTAGAATTACTACGGTTATCC--G-AGTAGTTACCATC-AACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT |||||||||||||||||||||||||||||||||||||||**|*|||||||||||||*|||||||||||||||||||||||||||||||||||||||||||||| TTGTTGCACGTATTAGCTCTAGAATTACTACGGTTATCCGAGTAGTAGTTACCATCAAACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT
#m54270_210831_092625/7013299/ccs/0_109 109 10 109 +  18S 1808 76 179 - -475 99 0 0 4 254 

TTGTTGCACGTATTAGCTCTAGAATTACTACGGTTATCC--G-AGTAGTTACCATC-AACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT 
|||||||||||||||||||||||||||||||||||||||**|*|||||||||||||*|||||||||||||||||||||||||||||||||||||||||||||| 
TTGTTGCACGTATTAGCTCTAGAATTACTACGGTTATCCGAGTAGTAGTTACCATCAAACAAACTATAACTGATTTAATGAGCCATTCGCAGTTTCACAGTCT




"""
