# HBV_model
This is the simulation code used to model the infection of single cell with Hepatitis B virus, corresponding to the paper submitted to Viruses.
This stochastic model uses a Gillespie algorithm.
The code is written in Fortran 90 and can be compiled using gfortran.  The -O4 flag should be used to optimise the executable output.  
This code as currently configured runs 100 repeats of the treatment-free control case, using integer seeds drawn from the seed.list file.
The seed list should be refreshed for each change in the variables to obtain statistically separable datapoints.
The results are written to the file "fort.9"

On a desktop computer (~2Ghz processor) the code should take around 25hrs to run and uses <1Gb of RAM.  
The output file fort.9 should be ~ 6Gb.

The code is available for use by anyone, however if in a academic context, we would appreciate a citation; the details of which will be added following publication of the paper.  
