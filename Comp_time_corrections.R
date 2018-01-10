# Case 1 - Differing valus of K
rm(list=ls())
library(compiler)
enableJIT(3)

library(tbart)
library(ecp)

# Simulate data
source("sim_data.R")

###################
# load in functions - for segmentation
###################
# DOUBLE CUSUM 
source("bs_DC.R")
# INDEP
source("PELT.R")
source("ind.R")
# MRC
source("mrc_mean.R")
source("multiple_mrc.R")
# MV
source("PELT_MV.R")

n = 500
N = 100

############################################
#set.seed(1)
############################################

# K = 1
true.chpts = vector("list",10)
ind.chpts = vector("list",10)
dcbs.chpts = vector("list",10)
mrc.chpts = vector("list",10)
agg.chpts = vector("list",10)
ecp.chpts = vector("list",10)

comp_time = matrix(nrow=10 , ncol = 5)

for (i in 1:10){

  set.seed(i)
  sim = simulate.data( K=5 , eps=1.0 )
  data = sim$data
  
  # true MRC's each series
  true.chpts[[i]] = sim$mrc[ sim$series.mrc ]
  
  # AGG
  start_time <- Sys.time()
  agg = apply( data , 2 , sum )
  pagg = PELT( agg , 200*log(dim(data)[2]) )
  agg.chpts[[i]] = rep( rev( pagg$cpts )[1] , N )
  end_time <- Sys.time()
  comp_time[i,1] = end_time - start_time
  
  # ECP
  start_time <- Sys.time()
  ecp.div = e.divisive( t(data) , alpha = 1)
  ecp.chpts[[i]] = rep( rev( ecp.div$estimates )[2] , N )
  end_time <- Sys.time()
  comp_time[i,2] = end_time - start_time
  
  # IND
  start_time <- Sys.time()
  ind.chpts[[i]] = ind( data )
  end_time <- Sys.time()
  comp_time[i,3] = end_time - start_time
  
  # DCBS
  start_time <- Sys.time()
  dcbs.chpts[[i]] = Bin_seg(data,10)
  end_time <- Sys.time()
  comp_time[i,4] = end_time - start_time
  
  # MRC
  start_time <- Sys.time()
  mrc = mrc.mean( data , beta = 1.5*log(n) )
  c = multiple.mrc( mrc , pmax=10 )
  p.hat = c$MDL
  mrc.chpts[[i]] = c$locs[[p.hat]][ c$affected[[p.hat]] ]
  end_time <- Sys.time()
  comp_time[i,5] = end_time - start_time
  
  print(i)
}

# agg , mv ecp
times = apply( comp_time , 2 , mean)

names(times) = c("AGG","ECP","IND","DCBS","MRC")
t( data.frame(times) )
library(xtable)
xtable(t( times ) )

## comp time -- MRC
comp_time = numeric(15)

for (i in 1:15){
print(i)  
  set.seed(1)
  sim = simulate.data( K=5 , eps=1.0 )
  data = sim$data
  
  # MRC
  start_time <- Sys.time()
  mrc = mrc.mean( data , beta = 1.5*log(n) )
  c = multiple.mrc( mrc , pmax=i )
  p.hat = c$MDL
  end_time <- Sys.time()
  comp_time[i] = end_time - start_time
  
}

pdf("Kmax_comp_cosr.pdf")
plot( 2:15 ,  comp_time[2:15]  , ylab = "Computational cost (seconds)" , xlab = expression(K[max]) )
dev.off()
