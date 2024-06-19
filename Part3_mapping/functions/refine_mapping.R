
refine_map_mt = function(map_mt, source, target){
  # source: a column vector of length q
  # target: a column vector of length p
  # map_mt: q-by-p map_mtine matrix 
  
  p = length(target)
  q = length(source)
  
  A = matrix(c(sum(target^2), sum(target),sum(target),p),nrow=2, ncol=2)
  A_inv = solve(A)
  newMap_mt = map_mt
  
  for(i in 1:q){
    u = matrix(map_mt[i,],ncol=1)
    x = matrix(c(source[i] - t(target)%*%u,1- sum(u)),ncol=1)
    newMap_mt[i,] = as.numeric(u + cbind(target,1)%*%A_inv%*%x)
  }
  
  return(newMap_mt)
  
}