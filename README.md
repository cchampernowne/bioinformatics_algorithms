This code runs a local alignment using the Smith-Waterman algorithm, and a global alignment using the Needleman-Wunsch algorithm. 
Sequences must have length between 5 - 50

The scoring function used on both algorithms is:
s(i,j) = {+2 if i and j match, -1 if mismatch}, gamma(g) = -2

A sample input file, 'sequences.in', has been included 

In order to run this code: 
- python 3.7.6 or higher should be enabled
- the input file must be located in the same folder as the executable
- the command prompt location should be in the same folder as the executable

to run local alignment on the input file, use the following command:
python local-alignment.py <input_file_name>

to run global alignment on the input file, use the following command:
python global-alignment.py <input_file_name>