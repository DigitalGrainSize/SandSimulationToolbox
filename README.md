 sand_simulation toolbox

MATLAB and Fortran software to simulate realistic models of packed sediment/granular beds consisting of individual discrete grains. The grains possess a realistic size and shape distribution. User is able to control these attributes and others such as packing and structure according to the underlying stochastic process which governs the placement of sand centres.

This program implements the algorithms of:
Buscombe, D. and Rubin, D.M., 2012, Advances in the Simulation and Automated Measurement of Well-Sorted Granular Material, Part 1: Simulations. Journal of Geophysical Research - Earth Surface 117, F02001.

To run the program, download and unzip, open MATLAB, cd to the master directory, then in the command window type:

simgrains_demo

to run the demo. This will implement the main program (simgrains.m) using each of the supplied config files which control inputs and outputs. Please note the config files have the UNIX/MAC file separator '/'. Windows users will need to replace / with \
These config files are:
sim.config
sim_from_input_coords.config
sim_from_input_image.config
sim_using_model1.config
sim_using_model4.config
sim_using_model5.config

Optimized performance:
The toolbox calculations are considerably speeded up if the particle centres, when a model is requested to generate the points, are generated using the Fortran code within ./cvt_f. The software will run without it, but the matlab implementation is considerably slower. You will need a fortran compiler to compile these programs. Please see the shell scripts ./cvt_f/cvt.sh and ./cvt_f/cvt_allmodels.sh for an example of how to compile with gfortran compiler. These models are invoked using the config file variables:
use_compiled=1; 
compiled_file='./cvt_f/cvt'
or 
compiled_file='./cvt_f/cvt_allmodels'

Written by Daniel Buscombe, various times in 2011-2013
while at
School of Marine Science and Engineering, University of Plymouth, UK
and now:
Grand Canyon Monitoring and Research Center, U.G. Geological Survey, Flagstaff, AZ 

Please contact:
dbuscombe@usgs.gov

to report bugs and discuss the code, algorithm, collaborations

For the latest code version please visit:
https://github.com/dbuscombe-usgs

See also the project blog: 
http://dbuscombe-usgs.github.com/

Please download, try, report bugs, fork, modify, evaluate, discuss. Thanks for stopping by!
