
% EXAMPLE USAGE: simgrains_demo
% 
% Written by Daniel Buscombe, various times in 2011 - 2013
% while at
% School of Marine Science and Engineering, University of Plymouth, UK
% then
% Grand Canyon Monitoring and Research Center, U.G. Geological Survey, Flagstaff, AZ 
% please contact:
% dbuscombe@usgs.gov
% for lastest code version please visit:
% https://github.com/dbuscombe-usgs
% see also (project blog):
% http://dbuscombe-usgs.github.com/
% Buscombe, D. and Rubin, D.M., 2012, Advances in the Simulation and Automated Measurement of Well-Sorted Granular Material, Part 1: Simulations. 
% Journal of Geophysical Research - Earth Surface 117, F02001.
%====================================
%   This function is part of 'sand simulation toolbox' software
%   This software is in the public domain because it contains materials that originally came 
%   from the United States Geological Survey, an agency of the United States Department of Interior. 
%   For more information, see the official USGS copyright policy at 
%   http://www.usgs.gov/visual-id/credit_usgs.html#copyright
%====================================
clc
disp('=============================================')
disp('======= SAND SIMULATION TOOLBOX =============')
disp('========= By Daniel Buscombe ================')
disp('Some demos illustrating some example uses of this toolbox')
disp('see the supplied config files for settings')
disp('see the following paper for more details:')
disp('Buscombe, D. and Rubin, D.M., 2012,')
disp('Advances in the Simulation and Automated Measurement of Well-Sorted Granular Material, Part 1: Simulations.')
disp('Journal of Geophysical Research - Earth Surface 117, F02001')
disp('=============================================')


disp('=========DEMO 1==============')
disp('200 grains using a supplied image, grain concentration of 0.7, no modified throats')
disp('saving slices and polytopes, printing surface')
disp('=============================================')
simgrains('sim.config')
clc

disp('===========DEMO 2============')
disp('1000 grains using a supplied file, grain concentration of 0.7, no modified throats')
disp('saving slices and polytopes, printing 3d model')
disp('=============================================')
simgrains('sim_from_input_coords.config')
clc

disp('=========DEMO 3==========')
disp('200 grains using a supplied image, grain concentration of 0.5, modified throats')
disp('saving slices and polytopes, printing surface and 3d model')
disp('=============================================')
simgrains('sim_from_input_image.config')
clc

disp('==========DEMO 4==============')
disp('1000 grains using a pvt model using a compiled fortran program, grain concentration of 0.7, modified throats')
disp('saving slices and polytopes, printing 3d-slice composite')
disp('=============================================')
simgrains('sim_using_model1.config')
clc

disp('=========DEMO 5=================')
disp('200 grains using a CP model in matlab, grain concentration of 0.7, no modified throats')
disp('saving slices and polytopes, printing 3d slice composite')
disp('=============================================')
simgrains('sim_using_model4.config')
clc

disp('========DEMO 6================')
disp('200 grains using a Strauss model in matlab, grain concentration of 0.7, no modified throats')
disp('saving slices and polytopes, printing 3d model')
disp('=============================================')
simgrains('sim_using_model5.config')
clc




