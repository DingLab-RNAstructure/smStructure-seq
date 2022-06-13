"""

 We parse the each bit vector file and merge

 - we append a serial name  


singularity exec forgi-centroid-mea-final.simg  python3 merge_bitvectors.py   --bit_file R1_5p.bit   R2_5p.bit    --output_file   merged_R1_R2.bit


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
    parser.add_argument("-b", "--bit_file",    type=argparse.FileType('r'), nargs='+', required=True)
    parser.add_argument("-o", "--output_file",  type=str, required=True)

    args = parser.parse_args()

    return args


if __name__ == "__main__":


   args = handler()
   ser = 1
   with open(args.output_file, 'w') as outf:
        for f in  args.bit_file: # iterate overs files
            for line in f:
                read_name, *vect = line.strip().split('\t')
                new_read_name = str(ser) + '_' + read_name

                pipe = [new_read_name] + vect
                outf.write('\t'.join(pipe) + '\n')  

                ser += 1

   print (f'Reads written: {ser-1}')
 
