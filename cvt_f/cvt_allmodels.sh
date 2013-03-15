#Written by Daniel Buscombe, various times in 2011-2013
#while at
#School of Marine Science and Engineering, University of Plymouth, UK
#and now:
#Grand Canyon Monitoring and Research Center, U.G. Geological Survey, Flagstaff, AZ 

#Please contact:
#dbuscombe@usgs.gov

#to report bugs and discuss the code, algorithm, collaborations

#For the latest code version please visit:
#https://github.com/dbuscombe-usgs

#See also the project blog: 
#http://dbuscombe-usgs.github.com/

#rm *.o #cvt
gfortran -c -g cvt.f90 

gfortran -c -g cvt_allmodels.f90 
gfortran cvt_allmodels.o cvt.o -o cvt_allmodels
rm *.o
# inputs
# 1 = num of grains (e.g. 1000)
# 2 = model (e.g. 1) pvt=1, cvt=2, hal=3
# 3 = batch (e.g. 10000) suggest: num_grains*10
# 4 = it_max (e.g. 20)
# 5 = outfile (e.g. 'n1000_out.txt')

./cvt_allmodels 100 1 10000 20 'n1000_out.txt' 
