# hwetest
Tests if a population is in the Hardy Weinberg Proportion (HWP) for the observed locus. <br/>
The conditional probability, under the Hardy-Weinberg
equilibrium, to obtain the sample X is compute as described by: <br/>
Howard Levene - "On a matching problem arising in genetics". <br/>
Annals of Mathematical Statistics. 1949; 20:91-94.<br/>

If the locus is biallelic, the function perform an Exact test, computing
the p-value of all possible tables and the summing all p-value<=p(observed
table). If you have downloaded ternplot (ID: 2299) the routine plots a De
Finetti's Diagram; if not the routine will try to download it from FEX.
If the locus is m-allelic (m>2), the function uses a Monte Carlo conventional
method to evaluate the p-value. 

Syntax: pvalue=hwetest(x,verbose,delta,alpha)

Input: X - Genotype matrix. If the locus is biallelic, X is a vector
          x=[AA AB BB]; else if the locus is m-allelic X is a lower
          triangular matrix of size=[m m]. If X is not a lower
          triangular matrix it will be triangularized.
      VERBOSE (optional)- a logical variable to display more results and comments:
             0=does not display (default)
             1=display 
      DELTA and ALPHA (optional)- If Monte Carlo method is used (if locus is more
          than bi-allelic), it is necessary to evaluate how many times 
          the process must be reiterated to ensure that p-value is 
          within DELTA units of the true one with (1-ALPHA)*100% confidence. 
          (Default DELTA=ALPHA=0.01).
Output: the probability that the population is in HWP
       the De Finetti's Diagram if the locus is biallelic
       if VERBOSE:
          Polymorphism Information Content (PIC)
          Matching probability
          Power of discrimination
          Power of exclusion
          Typical Paternity Index

Example: 
         Run hwedemo

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2007) HWtest: a routine to test if a locus is in Hardy
Weinberg equilibrium (exact test). 
http://www.mathworks.com/matlabcentral/fileexchange/14425
