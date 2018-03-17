#######################################################
# This snippet of code is to generate a Kth dimensional 
# Symmetric first or second order matrix

# input
# integer K for dimension of matrix
# order

OrderMat_1stOR2nd<- function(K, order = c('first','second') ){
  
  order <- match.arg(order)  # <-- this means it'll accept 'fir' or 'sec' arguments
  mat.2nd.ord <- mat.1st.ord<- matrix(0, K, K)
  
  ## order == 'first' create first order matrix
  if(order == 'first'){
    for(i in 1:(K-1)){   
      for(j in (i+1):K ){
        
        if(abs(i-j) == 1){  # <-- 1st order neighbours
          mat.1st.ord[i,j]<- mat.1st.ord[j,i]<- 1
        }
      }
    }
    mat.op<- mat.1st.ord
    
  }else{  # otherwise order == 'second'
    for(i in 1:(K-1)){   
      for(j in (i+1):K ){
        if(abs(i-j) <= 2){  # <-- 2nd order neighbours (?)
          mat.2nd.ord[i,j]<- mat.2nd.ord[j,i]<- 1
        }
      }
    }
    mat.op<- mat.2nd.ord
  }

  mat.op
  return(mat.op)  
}
