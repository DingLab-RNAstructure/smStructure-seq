#!/opt/software/conda/envs/myenv/bin/python3

"""

Write the fourth line as fasta with a given tag: 

python3 fold2dotbracketFasta.py   --input_file  bit_1.fold  --tag bit_1  --output_file  bit_1.db

or 

singularity exec ~/BUILD/FORGI-SING/forgi-centroid-mea.simg python3 fold2dotbracketFasta.py   --input_file  bit_1.fold  --tag bit_1  --output_file  bit_1.db

"""

import sys  
import os
from pathlib import Path
import argparse

def handler():
    """
    Get command line inputs
    """ 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file",   type=str, required=True)
    parser.add_argument("-t", "--tag",          type=str, required=True)
    parser.add_argument("-o", "--output_file",  type=str, required=True)

#    parser.set_defaults(
#        input_file   = 'bit_1.fold',
#        output_file  = 'bit_1.db',
#        tag  = 'bit_1',
#    )
#
    args = parser.parse_args()

    return args


def run():
    """ 
    Main function 
    """
    args = handler()

    if not os.path.exists(args.input_file):
       print(f'Input file not found: {args.input_file}')
       exit()

    
    if not args.tag:
       print(f'Tag not defined : {args.tag}  ')
       exit()
  
    with open(args.input_file) as inp, open(args.output_file,'w') as  db_outf:
       for t, line in enumerate(inp,1):
          if t == 4:
             db_outf.write('>'+ args.tag + '\n') 
             db_outf.write(line.strip() + '\n') 


if __name__ == "__main__":
 
    run()

