#########################################################
# This snippet of code will produce a random binary symmetric
# matrix
#
# input
# positive integer K - dimension of square matrix
# probabilility of elements being equal to 1
# SEED value in order to reproduce random matrix



SparseRandom_Mat<- function(K, prob = 0.5, SEED = 123){
  
  set.seed(SEED)
  
  mat.rand = matrix(0, K, K)
  
  for(i in 1:(K-1)){   
    for(j in (i+1):K ){
      mat.rand[i,j]<-mat.rand[j,i] <- rbinom(1,1,prob = prob)
    }
  }
  
  mat.rand
  return(mat.rand)  
}
