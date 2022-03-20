### **Local web application**

In order to run this code: 
* python 3.7.6 or higher should be enabled
* required libraries:
  * flask
  * numpy
  * re
  * pandas

Steps to run:
1. via command-line, navigate to your repository
2. use command: <python3 server.py>
3. open a web browser
4. navigate to: http://127.0.0.1:5000/ or the url specified by your command-line


### **Smith-Waterman via command-line:**
This code runs a local alignment using the Smith-Waterman algorithm. 
Sequences must have length between 5 - 50

The scoring function used is:
s(i,j) = {+2 if i and j match, -1 if mismatch}, gamma(g) = -2

In order to run this code: 
* python 3.7.6 or higher should be enabled
* the input file must be located in the same folder as the executable
* the command prompt location should be in the same folder as the executable
* required libraries:
  * numpy
  * re

Steps to run:
1. via command-line, navigate to /algorithms within your repository
2. use command: <python3 algorithms/local-alignment.py>
3. enter two sequences to be aligned, using the following formatting requirements:
  * sequence lengths must be between 5 - 50
  * sequences must not have spaces, or any other characters between nucleotides
  * nucleotides must be represented using only the characters A, T, C, G, U, a, t, c, g, and u


### **Needleman-Wunsch via command-line:**

This code runs a global alignment using the Needleman-Wunsch algorithm. 
Sequences must have length between 5 - 50

The scoring function used is:
s(i,j) = {+2 if i and j match, -1 if mismatch}, gamma(g) = -2
  
In order to run this code: 
* python 3.7.6 or higher should be enabled
* the input file must be located in the same folder as the executable
* the command prompt location should be in the same folder as the executable
* required libraries:
  * numpy
  * re

Steps to run:
1. via command-line, navigate to /algorithms within your repository
2. use command: <python3 algorithms/global-alignment.py>
3. enter two sequences to be aligned, using the following formatting requirements:
  * sequence lengths must be between 5 - 50
  * sequences must not have spaces, or any other characters between nucleotides
  * nucleotides must be represented using only the characters A, T, C, G, U, a, t, c, g, and u


**UPGMA via the command-line:**


**Hidden Markov machine via command-line:**


**Nussinov via command-line:**



