"""
Generate unique Folding candidates given the bit vectors

singularity exec forgi-centroid-mea-final.simg  python3 fold-contrafold-uniq-bits-vectors.py   --bit_file merged_R1_R2.bit --reference_file  cool6.fasta --transcript COOLAIR3 --size_file  sizer.tab 


"""

import sys
import os
import argparse
import jitu
from collections import defaultdict
import subprocess


def generate_constraints(ref, profile):
    """
        2. A nucleotide mapping to 0 is constrained to be unpaired.
        3. A nucleotide mapping to -1 is unconstrained.
        For example, given the following input BPSEQ file:
        1 A -1
        2 C -1
        3 G -1
        4 U 7
        5 U 0
        6 C 0
        7 G 4
        8 C -1
    """
    state = []
    for pos, (base, bit) in enumerate(zip(ref, profile), 1):
        if bit == '1':  # if bit vector mutation has one: single stranded, we signal contrafold to be unpaired 0
            state.append([str(pos), base, '0'])
        else:
            state.append([str(pos), base, '-1'])  # unconstrained
    return state


def handler():
    """
    Get command line inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bit_file", type=str, )  # required=True)
    parser.add_argument("-r", "--reference_file", type=str, )  # , required=True)
    parser.add_argument("-t", "--transcript", type=str, )  # , required=True)
    parser.add_argument("-s", "--size_file", type=str, )  # , required=True)

    parser.set_defaults(bit_file='merged_R1_R2.bit',
                        reference_file='cool6.fasta',
                        transcript='COOLAIR3',
                        size_file='sizer.tab')

    return parser.parse_args()


if __name__ == "__main__":

    args = handler()

    os.makedirs('constr', exist_ok=True)
    os.makedirs('folded', exist_ok=True)
    os.makedirs('posteriors_output', exist_ok=True)

    tube, seqD = jitu.getTubeD(fastaFiler=args.reference_file)
    reference_sequence = seqD[args.transcript]

    tot = 0
    unique_bits = {}  # unique bit vectors
    profileD = defaultdict(int)
    with open(args.bit_file) as inp:
        for line in inp:
            A = line.strip().split('\t')
            if A:
                bits = ''.join(['1' if b in '1' else '.' for b in A[1:]])
                unique_bits[bits] = bits.count('1')  # number of observed mutations or ones
                profileD[bits] += 1
                tot += 1
    print(f'unique bits vectors: {len(unique_bits)}  from total: {tot}')
    
    seen = defaultdict(int)

    with open('sizer.tab', 'w') as outS, open('forgi-vect-ser.txt', 'w') as serF:

        for t, (profile, num_mutations) in enumerate(sorted(unique_bits.items(), key=lambda x: x[1], reverse=True),
                                                     1):  # sorted by decreasing num mutants
            assert len(reference_sequence) == len(profile), 'length mismatch bad reference !'
            state = generate_constraints(reference_sequence, profile)
            print(f'{profile[:10]}  {num_mutations} ')

            ticks_string = ''.join([x[-1] for x in state])  # 1-1-1-1-1-1-1-1001-1-1-1-1-1-1-1
            seen[ticks_string] += 1

            if seen[ticks_string] == 1:
                constraint_file = 'constr/' + 'bit_' + str(t) + '.bpseq'
                bit_prefix = 'bit_' + str(t)
                with open(constraint_file, 'w') as outf:
                    for i, base, tick in state:
                        outf.write('\t'.join([i, base, tick]) + '\n')
                
                pipe = ["contrafold", "predict", "--constraints", constraint_file, '--parens',
                        'folded/' + bit_prefix + '.fold',
                        '--posteriors', '0.0', 'posteriors_output/' + bit_prefix + '.txt']

                #print(' '.join(pipe))
                result = subprocess.run(
                    pipe, capture_output=True, text=True
                )
                print("stdout:", result.stdout)
                print("stderr:", result.stderr)
               
                # dotbracket                 
                pipe = [  "fold2dotbracketFasta.py", "--input_file", 'folded/' + bit_prefix + '.fold',
                          "--tag", bit_prefix,
                          "--output_file", 'folded/' +  bit_prefix + '.db']

                print(' '.join(pipe))
             
                result = subprocess.run(
                    pipe, capture_output=True, text=True
                )
                print("stdout:", result.stdout)
                print("stderr:", result.stderr)

                # rnaconvert        
                pipe = ["rnaConvert.py", 'folded/' +  bit_prefix + '.db', '-T', 'element_string', '--force', '--to-file', '--filename', 'folded/' +  bit_prefix + '.txt']

                print(' '.join(pipe))

                result = subprocess.run(
                    pipe, capture_output=True, text=True
                )
                print("stdout:", result.stdout)
                print("stderr:", result.stderr)
               
                # adding to forgi strings
                filer = 'folded/' +  bit_prefix + '.txt001' + '.element_string'
                
                with open(filer) as inp:
                    for t, line in enumerate(inp,1):
                        if t == 3:
                            serF.write('\t'.join([bit_prefix, line.strip()]) + '\n')  

                
                #------- sizer information this bit profile is represented by 
                outS.write('\t'.join([bit_prefix, str(profileD[profile]), profile] ) + '\n')

    print('DONE')            
                
"""

unique bits vectors: 802  from total: 3005

"""