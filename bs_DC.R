## DOUBLE CUSUM ##
CUSUM.dims = function( data ){
  
  n = dim(data)[2]
  i = 1:(n-1)
  N = dim(data)[1]  
  X = matrix( nrow=N , ncol=n-1 )
  
  for (j in 1:N){
    cs = cumsum( data[j,] )
    X[j,] = abs( sqrt( ( i*(n-i) )/n ) * ( cs[i]/i  - ( cs[n] - cs[i] )/(n-i) ) )      
  }

  # need to order by time pt b 
  Xord = apply( X , 2 , order )
  X = apply(X,2,sort)
  
  newlist = list("X"= X , "Xord" = Xord)
  return(newlist)  
  
}

# inout matrix from CUSUM.dims
DC = function( newlist ){
  
  X = newlist$X
  Xord = newlist$Xord
  
  # X is ordered
  N = dim(X)[1]
  c = numeric(N)
  c_max = numeric(N)
  for (m in 1:N){

    if ( m==1 ){
      y = sqrt( ( m * ( 2*N - m )/(2*N) ) ) * (  X[ N:N , ]/m - apply( X[ 1:(N-1) , ] , 2 , sum )/(2*N - m) ) ##added brackets
      c_max[(m)] = which.max(y)-1
      c[(m)] = max(y)
    }
    
    else if ( m==N ){
      y = sqrt( ( m * ( 2*N - m )/(2*N) ) ) * apply( X[ (N-m+1):N , ] , 2 , sum )/m # - X[1,]/(2*N - m) #removed this -- empty sum ==0 I assume
      c_max[(m)] = which.max(y)-1
      c[(m)] = max(y)
    }
  
    else if (m==(N-1)){
      y = sqrt( ( m * ( 2*N - m )/(2*N) ) ) *( apply( X[ (N-m+1):N , ] , 2 , sum )/m - X[1:(N-m),]/(2*N - m) ) ##added brackets
      c_max[(m)] = which.max(y)-1
      c[(m)] = max(y)
    }
    
    else{
      y = sqrt( ( m * ( 2*N - m )/(2*N) ) ) *  (apply( X[ (N-m+1):N , ] , 2 , sum )/m - apply( X[ 1:(N-m) , ] , 2 , sum )/(2*N - m) ) ##added brackets
      c_max[(m)] = which.max(y)-1
      c[(m)] = max(y)
    }
    
  }
  
  
  # no of series change 
  mb = which.max( c )
  # location of change
  tau = c_max[which.max(c)]
  # which series change
  series.changes = Xord[ N:(N-mb+1) , tau ]   ##I think this should have N-mb+1
  
  # c_max is location of changes
  newlist = list( "max" = max(c) , "tau" = tau , "series" = series.changes )
  return(newlist)
  
}


# DCBS algo

Bin_seg = function(data,beta){

  n = dim(data)[2]  
  N = dim(data)[1]
  S = c(1,n) 
  CP = c()
  affected.series = list()
  
  while ( length(S) > 0 ){
    
    # between start point and end point
    start = S[1]
    end = S[2]
  
    if ( (end-start)>1 ){
      
      X = CUSUM.dims( data[ , start:end ] )
      test.stat = DC(X)
      
      # if changepoint is significant add segments to end of S
      chpt = test.stat$tau + start
      
      if ( (test.stat$max > beta) & !(chpt %in% CP) ){
        
        CP = c( CP , chpt )
        affected.series[[ length(affected.series) + 1]] = test.stat$series
        
        if (chpt != start + 1){
          S = c( S , c( start , chpt ) )       
        }
        if (chpt != end - 1){
          S = c( S , c( chpt + 1 , end ) )
        }
      }
      
    }
    
    S = S[-c(1,2)]
   
  }

  # find mrc chpt for each series
  mrc.chpt = matrix( 0 , nrow = max(1,length(CP)) , ncol = N )
  mrc = rep(0,N)
  # check whether no of affected series is more than one
  if (length(affected.series)>0){
    for (i in 1:length(affected.series)){
      mrc.chpt[ i , affected.series[[i]] ] = CP[i]
    }
    mrc = apply( mrc.chpt , 2 , max )
  }
  #newlist = list( "CP" = CP , "series" = affected.series  )
  return(mrc)
  
}



