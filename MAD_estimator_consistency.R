## MAD estimator and SD

rm(list=ls())
load("alpha_change.Rdata")


#https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/
# setting the consistency constant equal to the reciprocal of the 75th percentile of the standard distribution will achieve the nice estimation property mentioned above. 
# By "standard" I mean the distribution that results from shift-and-scale tranforming 
# the distribution in question to a distribution 
# with a mean of zero and a standard deviation of 1.

for (i in 1:7039){
  i=1
  z = ( data_alpha[[i]] - mean(data_alpha[[i]]) )/sd( data_alpha[[i]] ) 
  
  

}
