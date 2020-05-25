# PicXAA-R

Official implementation of **P**robabilist**I****C** ma**X**imum **A**ccuracy **A**lignment of **R**NA sequences
            
-----------------------------------------------------------------
PicXAA-R is a probabilistic non-progressive alignment algorithm
that finds multiple sequence alignments with maximum
expected accuracy. PicXAA greedily builds up the multiple
alignment from sequence regions with high local similarities,
thereby yielding an accurate global alignment that 
effectively grasps the local similarities among sequences. 

PicXAA-R is an extension to PicXAA for greedy structural 
alignment of ncRNAs. PicXAA-R efficiently grasps both 
folding information within each sequence and local similarities 
between sequences. Using a graph-based scheme, it greedily build
up the structural alignment from sequence regions with high 
base-pairing and base alignment probabilities

PicXAA-R is one of the fastest algorithms for structural alignment 
of multiple RNAs and it consistently yields accurate alignment
results, especially for datasets with locally similar sequences

PicXAA-R Version 1.0 (Jan 2011) is developed by Sayed Mohammad 
Ebrahim Sahraeian and Byung-Jun Yoon in GSP lab, Department 
of Electrical and Computer Engineering, Texas A&M University. 

PicXAA-R uses codes from PicXAA version 1.02, PROBCONS version 
1.12 (written by Chuong Do), MXSCARNA version 2.1 (written by 
Yasuo Tabei), and the energy parameters of Vienna RNA package version 1.5
(written by Ivo Hofacker in Institute for Theoretical Chemistry 
of the  University of Vienna).

For more information on the algorithms, please see

Sahraeian, S.M.E., Yoon, B.J. (2011), PicXAA-R: efficient 
structural alignment of multiple RNA sequences using a greedy 
approach, BMC Bioinformatics, In Press.


Sahraeian, S.M.E., Yoon, B.J. (2010), PicXAA: Greedy 
Probabilistic Construction of Maximum Expected Accuracy 
Alignment of Multiple Sequences, Nucleic Acids Research,
38(15): 4917-4928, 2010.

-----------------------------------------------------------------
PicXAA-R has been made freely available as PUBLIC DOMAIN
software and hence is not subject to copyright in the  United
States.  This system and/or any portion of  the  source  code
may be used, modified, or redistributed without restrictions.  
PicXAA is distributed WITHOUT WARRANTY, express or implied.
The authors accept NO LEGAL LIABILITY OR  RESPONSIBILITY  for
loss due to reliance on the program.
-----------------------------------------------------------------

## Install:
       make clean
       make
## Usage:
       ./picxaa-r [options]  MFAFILE
## Example:
       ./picxaa-r test/sample1.fasta > test/output.fasta

## Description:
       Align sequences in MFAFILE(s) and print result to standard output

       -stockholm
              use STOCKHOLM output format instead of MFA

       -clustalw
              use CLUSTALW output format instead of MFA

       -noss
              do not show the consensus secondary structure

       -columnscore
              write annotation for multiple alignment

       -v, --verbose
              report progress while aligning (default: off)

       -a, --alignment-order
              print sequences in alignment order rather than
              input order (default: off)

       -c, --consistency REPS
              use 0 <= REPS <= 5 (default: 1) passes of
              inter-sequence consistency transformation

       -ic, --intraconsistency REPS
              use 0 <= REPS <= 5 (default: 1) passes of
              intra-sequence consistency transformation

       -al, --alpha
              Sets the weight parameter (alpha) for 
              intra-sequence consistency transformation 
              (0<alpha<1 , default: 0.4)
       -bt, --beta
              Sets the weight parameter (beta) for 
              Four-Way consistency transformation 
              (0<beta<1 , default: 0.1)
       -Tb
              Sets the threshold value for identifying
              the most probale base-pairing probabilities 
              (0<Tb<1 , default: 0.5)
       -r, --refinement REPS
              use 0 <= REPS <= 0 (default: 100) passes of 
              refinement realignments

       -p, --paramfile FILENAME
             read parameters for Pai-HMM from FILENAME 
              (default: )

## Version History
1.0, 1/01/2011
   -- source code released