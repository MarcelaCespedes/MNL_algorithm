#########################################################
# This is the same as MNL_v1.r only here we try to use a faster 
# dmnorm() implementation via C++

# using dnorm() loping thru all 70 ROI off diagonals takes 4 minutes!!!!
# using the C++ version this takes 0.0886 minutes ==> 5 seconds!!!!

MNL_v2<- function(init.s2s, init.W, K=10, spatDat, rho.val = 0.9, no.runs = 10){
  
  library(mvtnorm)
  library(stats) # <== for optim()
  
#  library(RcppArmadillo)
#  library(Rcpp)
#  library(installr) 
#  install.Rtools()  # <== test for Rcpp
#  evalCpp("1+1") 
  
  Rcpp::sourceCpp('dmvnrm_arma.cpp')
  
  I = dim(spatDat)[1]
  I
  p <- rep(0,I)
  p_star <- p
  w.keep<- array(0, dim = c(K, K, no.runs))
  w.keep[,,1]<- init.W
  s2s.vec<- rep(0, no.runs)
  s2s.vec[1]<- init.s2s
  
  # this is the function to optimise and solve for s2s; given a particular inv.Q
  loglike <- function(s2s){
    p = rep(0, I)
    for (l in 1:I){
      p[l]<- dmvnorm(x=spatDat[l,], mean=rep(0,K), sigma=s2s*inv.Q, log=TRUE)
    }
    return(-sum(p))  # <==== why is this negative? because it's is the negative of the log-density
    #return(sum(p))  
  }
  
  
  for (r in 2:no.runs){
    
    w = w.keep[,,r-1]
    
    #tiMe<- proc.time()   
    for(i in 1:(K-1)){
      for(j in (i+1):K){
        
        Q.curr <- rho.val*(diag(K)*rowSums(w) - w) + (1-rho.val)*diag(K)
        inv.Q.curr <- solve(Q.curr)
        
        sigma.curr<- s2s.vec[r-1]*inv.Q.curr
        #for (l in 1:I){
        #  p1[l]<- dmvnorm(spatDat[l,], rep(0,K), s2s.vec[r-1]*inv.Q.curr,log=TRUE)
        #}
        p<- dmvnrm_arma(spatDat, mean = rep(0, K), sigma.curr, TRUE)
        
        mi <- sum(p)
        
        # permute W
        w_star <- w
        if(w_star[i,j] == 0){
          w_star[i,j] <- w_star[j,i] <-  1
        }else{
          w_star[i,j] <- w_star[j,i] <- 0
        }
        
        Q.star <- rho.val*(diag(K)*rowSums(w_star) - w_star) + (1-rho.val)*diag(K)
        inv.Q.star <- solve(Q.star)
        sigma.star<- s2s.vec[r-1]*inv.Q.star
        
        #for (l in 1:I){
        #  p_star[l] <- dmvnorm(spatDat[l,], rep(0,K), s2s.vec[r-1]*inv.Q.star, log=TRUE)
        #}
        p_star = dmvnrm_arma(spatDat, mean = rep(0, K), sigma.star, TRUE)
        
        mi_star <- sum(p_star)
        
        if(mi_star>mi){  # if there is an improvement in the log-likelihood by permuting 
          w <- w_star    # the elements of W, then keep it
        }
        #if(runif(1)< exp(mi_star - mi)){w <- w_star}
        #w.keep[,,r] <- w
        
      }
      
      #print(i)
    }# end nested for-loop
    #(proc.time() - tiMe)/60 
    
    w.keep[,,r]<- w
    
    #Q <- rho.val*(diag(K)*rowSums(w) - w) + (1-rho.val)*diag(K)
    Q <- rho.val*(diag(K)*rowSums(w.keep[,,r]) - w.keep[,,r]) + (1-rho.val)*diag(K)
    inv.Q <- solve(Q)
    
    #out <- optim(s2s.vec[r-1], loglike, 
    #             method="BFGS", #lower=0.00001, upper=2,
    #             control=list(fnscale = -1, trace=TRUE) ) #optimisation algorithm
    
    invisible(capture.output(out<-optim(s2s.vec[r-1], loglike, 
                                        method="L-BFGS-B", lower=0.01, upper=100,
                                        control=list(trace=TRUE) ))) #optimisation algorithm
    s2s.vec[r] <- out$par
    
  }
  w.keep[,,r]
  s2s.vec
  
  #chk.W<- raster(xmn = 0, xmx = K, ymn = 0, ymx = K, nrows = K, ncols = K)
  #chk.W[]<- as.vector(w.keep[,,r])
  #chk.W[]<- as.vector(Q)
  #x11()
  #plot(chk.W, legend = F)
  
  
  ##############################################
  final.op<- c()
  Q.final <- rho.val*(diag(K)*rowSums(w.keep[,,r]) - w.keep[,,r]) + (1-rho.val)*diag(K)
  
  for (i in 1:I){
    final.op[i]<- dmvnorm(spatDat[i,], rep(0,K), s2s.vec[r]*solve(Q.final),log=TRUE)
  }
  final.log.lik = sum(final.op)
  ###############################
  return(list(w.keep = w.keep, s2s = s2s.vec, final.log.lik = final.log.lik))
}
