%HWEDEMO
%This is a demo on the use of HWETEST function
%Example
%Run hwedemo
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2007) HWtest: a routine to test if a locus is in Hardy
% Weinberg equilibrium (exact test). 
% http://www.mathworks.com/matlabcentral/fileexchange/14425

clc; home
disp('If the locus has only 2 alleles, you can input a vector x=[AA AB BB].')
x=[35 250 235];
disp(' ')
disp(x)
disp('Press a key to continue'); pause; clc; home
hwetest(x,1)
disp('Press a key to continue'); pause; close; clc; home
x=[12 10 4 3 3 3 0; 3 5 5 4 0 2 0; 10 5 1 0 0 3 0; 4 0 3 3 1 1 0; 0 4 4 0 0 0 0; 3 0 0 0 0 1 0; 0 0 0 0 0 0 0];
disp('This matrix show the genotypes of a population for a locus with 7 possible alleles')
disp(' ')
disp(x)
disp(' ')
disp('Note that X is a [7 7] square matrix')
disp('Press a key to continue'); pause; clc; home
hwetest(x,1)