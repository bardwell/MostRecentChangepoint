mrc.mean = function( hts , beta = 2*log(dim(hts)[2])  ){

  n = dim(hts)[2]

  # first find the dimensions which arent all zero
  non_null_dims = which( apply( hts , 1 , sum ) != 0 )

  # mean change in individual series
  PELT_hts = matrix( nrow = length(non_null_dims) , ncol = n )
  for ( i in 1:length(non_null_dims) ){
    d = non_null_dims[i]
    PELT_hts[i,] = PELT( hts[d,] , beta )$F[1:(n)] + vec_C( hts[d,] )
  }

  return(PELT_hts)

}

# vector of costs for (t+1):n FOR MEAN
vec_C = function(data){

  n = length(data)
  cd = cumsum( c(0,data) )
  cd_2 = cumsum( c(0,data^2) )
  t = 0:(n-1)
  seg_cost = cd_2[(n+1)] - cd_2[(t+1)] - ( ( cd[(n+1)] - cd[(t+1)] )^2 )/(n-t)
  return(seg_cost)

}

