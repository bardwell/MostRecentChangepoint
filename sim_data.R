## function to sim a single series with given chpts and eps
## new chpt sim function
sim_series_chpts = function( n , chpts , eps ){
  
  data = rnorm(n)
  mu = rnorm(1,0,2)
  data[1:chpts[1]] = data[1:chpts[1]] + mu
  if ( length(chpts) > 2 ){
    for (i in 2:( length(chpts) - 1) ){
      mu = rnorm(1,0,2)
      data[ (chpts[i] + 1):( chpts[(i+1)] ) ] =  data[ (chpts[i] + 1):( chpts[(i+1)] ) ] + mu
    }
  }
  data[ ( tail(chpts,1) + 1 ):n] = data[ ( tail(chpts,1) + 1 ):n] + mu + eps 
  return(data)
}

# ## simulate
# # length of time series
# n = 500
# # dimension
# N = 100
# # number of MRC's
# K = 10
# # mu + eps - mean of last seg 
# eps = 10

simulate.data = function( n=500 , N=100 , K=2 , eps=5 ){
  
  ### alternative K<=10###
  true.mrc.chpts = n-sample(20*(1:10) , K , replace = FALSE)
   
  # which series carry MRC's
  f = floor( N/K )
  # reorder series
  tsr = sample(1:N,N)
  
  # locations of ordinary chpts
  chpt.pot.locs = rbinom( min(true.mrc.chpts) , 1 , prob = 0.02) 
  chpt.locs = which( chpt.pot.locs == 1 )
  # prop of series each chpt affects
  alpha = runif(length(chpt.locs))
  
  chpts.each.series = vector("list",N)
  series.which.mrc = numeric(N)
  data = matrix(nrow=N,ncol=n)
  for (i in 1:N){
    
    # which of the chpts are in this series
    probs = runif(length(chpt.locs))
    wc = which( probs < alpha )
    
    # which most recent chpt is series affected by
    w = which(tsr == i)
    m = ceiling(w/f)
    if (m >K){
      m <- K
    }
    # which MRC affects ith series 
    series.which.mrc[i] = m 
    # changepoints in each series
    chpts.each.series[[i]] = c( chpt.locs[wc] , true.mrc.chpts[m] )
    
    data[i,] = sim_series_chpts( n , chpts.each.series[[i]] , eps )
    
  }
  
  newlist = list("data" = data , "mrc" = true.mrc.chpts , "series.mrc" = series.which.mrc , "series.chpts" =  chpts.each.series )
  return(newlist)
  
}





