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
true.chpts = vector("list",100)
ind.chpts = vector("list",100)
dcbs.chpts = vector("list",100)
mrc.chpts = vector("list",100)
agg.chpts = vector("list",100)

for (i in 1:100){
    set.seed(i)
  sim = simulate.data( K=1 , eps=1.0 )
  data = sim$data
  
  # true MRC's each series
  true.chpts[[i]] = sim$mrc[ sim$series.mrc ]
  
  # AGG
  agg = apply( data , 2 , sum )
  pagg = PELT( agg , 200*log(dim(data)[2]) )
  agg.chpts[[i]] = rep( rev( pagg$cpts )[1] , N )

  # MV
  pmv = PELT.MV( data , 101*log(dim(data)[2]) )  
  mv.chpts =  rep( rev( pmv$cpts )[1] , N )

  # ECP
  ecp.div = e.divisive( t(data.train) , alpha = 1)
  ecp.chpts = rep( rev( ecp.div$estimates )[2] , N )
  ecp.mse[i] = pred.mse( ecp.chpts , data.train , data.test )
  
  # IND
  ind.chpts[[i]] = ind( data )
  
  # DCBS
  dcbs.chpts[[i]] = Bin_seg(data,10)
  
  # MRC
  mrc = mrc.mean( data , beta = 1.5*log(n) )
  c = multiple.mrc( mrc , pmax=10 )
  p.hat = c$MDL
  mrc.chpts[[i]] = c$locs[[p.hat]][ c$affected[[p.hat]] ]

  print(i)
  }

save( true.chpts , ind.chpts , dcbs.chpts , mrc.chpts , file = "K1_resultsb.Rdata" )

############################################

# K = 2
true.chpts = vector("list",100)
ind.chpts = vector("list",100)
dcbs.chpts = vector("list",100)
mrc.chpts = vector("list",100)


for (i in 1:100){
     set.seed(i) 
  sim = simulate.data( K=2 , eps=1.0 )
  data = sim$data
  
  # true MRC's each series
  true.chpts[[i]] = sim$mrc[ sim$series.mrc ]
  
  # IND
  ind.chpts[[i]] = ind( data )
  
  # DCBS
  dcbs.chpts[[i]] = Bin_seg(data,10)
  
  # MRC
  mrc = mrc.mean( data , beta =  1.5*log(n) )
  c = multiple.mrc( mrc , pmax=10 )
  p.hat = c$MDL
  mrc.chpts[[i]] = c$locs[[p.hat]][ c$affected[[p.hat]] ]
  
  print(i)
}

save( true.chpts , ind.chpts , dcbs.chpts , mrc.chpts , file = "K2_resultsb.Rdata" )

############################################

# K = 3
true.chpts = vector("list",100)
ind.chpts = vector("list",100)
dcbs.chpts = vector("list",100)
mrc.chpts = vector("list",100)


for (i in 1:100){
     set.seed(i) 
  sim = simulate.data( K=3 , eps=1.0 )
  data = sim$data
  
  # true MRC's each series
  true.chpts[[i]] = sim$mrc[ sim$series.mrc ]
  
  # IND
  ind.chpts[[i]] = ind( data )
  
  # DCBS
  dcbs.chpts[[i]] = Bin_seg(data,10)
  
  # MRC
  mrc = mrc.mean( data , beta = 1.5*log(n)  )
  c = multiple.mrc( mrc , pmax=10 )
  p.hat = c$MDL
  mrc.chpts[[i]] = c$locs[[p.hat]][ c$affected[[p.hat]] ]
  
  print(i)
}

save( true.chpts , ind.chpts , dcbs.chpts , mrc.chpts , file = "K3_resultsb.Rdata" )

############################################

# K = 4
true.chpts = vector("list",100)
ind.chpts = vector("list",100)
dcbs.chpts = vector("list",100)
mrc.chpts = vector("list",100)


for (i in 1:100){
     set.seed(i) 
  sim = simulate.data( K=4 , eps=1.0 )
  data = sim$data
  
  # true MRC's each series
  true.chpts[[i]] = sim$mrc[ sim$series.mrc ]
  
  # IND
  ind.chpts[[i]] = ind( data )
  
  # DCBS
  dcbs.chpts[[i]] = Bin_seg(data,10)
  
  # MRC
  mrc = mrc.mean( data , beta = 1.5*log(n)  )
  c = multiple.mrc( mrc , pmax=10 )
  p.hat = c$MDL
  mrc.chpts[[i]] = c$locs[[p.hat]][ c$affected[[p.hat]] ]

}

save( true.chpts , ind.chpts , dcbs.chpts , mrc.chpts , file = "K4_resultsb.Rdata" )

############################################

# K = 5
true.chpts = vector("list",100)
ind.chpts = vector("list",100)
dcbs.chpts = vector("list",100)
mrc.chpts = vector("list",100)


for (i in 1:100){
  set.seed(i)
  sim = simulate.data( K=5 , eps=1.0 )
  data = sim$data
  
  # true MRC's each series
  true.chpts[[i]] = sim$mrc[ sim$series.mrc ]
  
  # IND
  ind.chpts[[i]] = ind( data )
  
  # DCBS
  dcbs.chpts[[i]] = Bin_seg(data,10)
  
  # MRC
  mrc = mrc.mean( data ,  beta = 1.5*log(n)  )
  c = multiple.mrc( mrc , pmax=10 )
  p.hat = c$MDL
  mrc.chpts[[i]] = c$locs[[p.hat]][ c$affected[[p.hat]] ]
  print(i)
}

save( true.chpts , ind.chpts , dcbs.chpts , mrc.chpts , file = "K5_resultsb.Rdata" )

############################################

# K = 10
true.chpts = vector("list",100)
ind.chpts = vector("list",100)
dcbs.chpts = vector("list",100)
mrc.chpts = vector("list",100)


for (i in 1:100){
      set.seed(i)
  sim = simulate.data( K=10 , eps=1.0 )
  data = sim$data
  
  # true MRC's each series
  true.chpts[[i]] = sim$mrc[ sim$series.mrc ]
  
  # IND
  ind.chpts[[i]] = ind( data )
  
  # DCBS
  dcbs.chpts[[i]] = Bin_seg(data,10)
  
  # MRC
  mrc = mrc.mean( data , beta = 1.5*log(n)  )
  c = multiple.mrc( mrc , pmax=20 )
  p.hat = c$MDL
  mrc.chpts[[i]] = c$locs[[p.hat]][ c$affected[[p.hat]] ]
  
  print(i)
}

save( true.chpts , ind.chpts , dcbs.chpts , mrc.chpts , file = "K10_resultsb.Rdata" )

