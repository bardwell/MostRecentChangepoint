multiple.mrc = function( G , pmax = 5 , alpha = 2 , elbow.thresh = 0.5 ){

  # no of series
  N = dim(G)[1]
  location.vec = 0:(n-1)

  cost = numeric( pmax )
  mmrc = vector( "list" , pmax )
  affected = vector( "list" , pmax )
  # p=1 separate as simpler
  mmrc[[1]] = tb.raw( G , c(1) ) 
  # say all series are affected by this 1 change
  index = rep(1,times=N)
  # find which affected (more evidence above threshold)
  #index[ G[ , mmrc[[1]] ] > - alpha ] <- 0
  affected[[1]] = index
  # objective cost
  cost[1] =  sum( G[ , mmrc[[1]] ] )

  if (pmax>1){

    for (p in 2:pmax){

      # mmrc[[i]] gives locations of the i best locations
      mmrc[[p]] = tb.raw( G , c(1:p) )
      # affected[[i]], gives each dimension a label from 1:i
      # depending on which change it is associated with
      affected[[p]] = apply( G[ , mmrc[[p]] ] , 1 , which.min )
      # cost[[i]] gives the objective cost for solving with i different changes/sets
      csum = 0
      for (i in 1:p){
        wa = which( affected[[p]] == i)
        csum = csum + sum( G[ wa , mmrc[[p]][i] ] )  
      }
      cost[p] = csum

    }

  }
  
  locations = vector("list",pmax)
  for ( i in 1:pmax ){
    locations[[i]] = location.vec[ mmrc[[i]] ]
  }

  BIC = which.min( cost + ( (1:pmax) )*log(N)*log(N*n) )  ###
  MDL = which.min(cost + (N*log(1:pmax)+ (1:pmax)*log(n))/log(2)) ##MDL

  # # selecting best p - Lavielle
  # J = numeric(pmax)
  # for (k in 1:pmax){
  #  J[k] = ( ( cost[pmax] - cost[k] )/( cost[pmax] - cost[1] ) ) * ( pmax - 1 ) + 1
  # }
  # D = numeric(pmax-1)
  # D[1] = Inf
  # for ( k in 2:(pmax-1) ){
  #   D[k] = J[k-1] - 2*J[k] + J[k+1]
  # }
  # elbow = max( which( D > elbow.thresh ) )
#"elbow" = elbow,
  newlist = list( "locs" = locations , "affected" = affected , "cost" = cost , "bic" = BIC,  "MDL"=MDL )
  return(newlist)

}

