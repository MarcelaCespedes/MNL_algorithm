######################################################
# MNL Simulation study for S_2 matrix/ data
#
#


library(MASS)
library(raster)
library(ggplot2)
library(glasso)

rm(list=ls())

#################################################
## To generate new binary connectivity matrix  ##
## uncomment the code below                    ##
#################################################

#SEED = 456
#source("OrderMat_1stOR2nd.R")
#source("SparseRandom_Mat.R")
#set.seed(SEED)
#K=70 # dimension of matrix
#no.reps<- 10  # each simulatin study will be for no.reps replicates

###
### Generate first or second order neighbours
#mat.1<- OrderMat_1stOR2nd(K=K, order='second')

###
### Generate random symmetric matrix
#prob = 0.05
#mat.2<- SparseRandom_Mat(K=K, prob = prob, SEED = SEED)

###
### combine mat.2nd.ord and mat.rand
#comb.mat<- mat.1 + mat.2
#comb.mat[1:10, 1:10]

#comb.mat<- ifelse(comb.mat == 2 | comb.mat == 1, 1, 0) # return to binary
#comb.mat[1:10, 1:10]

#tot.no.links<- sum(rowSums(comb.mat))/2
#tot.no.links

# on paper 1) they multiply these binary matrices by a constant, 0.04, and make the diagonal 
# a positive constant. Here I fiddle with different values

#off.diag.cov = 1
#f.mat<- off.diag.cov*comb.mat  # when this constant is small --> poor MNL results, better with high values 
#diag.cov <- 6
#diag(f.mat)<- diag.cov      # likewise when the above is large, then this has to be larger
#f.mat[1:10, 1:10]

#eigen(f.mat)$values

# If all eigenvalues are positive, then generate MVN dat:  100, 250, 500, 1000 
# Need to generate 10 sets of random data for each matrix config
#no.samp = 1000

#Sim.1<- list()
#for(i in 1:no.reps){
#  Sim.1[[i]] <- mvrnorm(n = no.samp, mu =rep(0, K), Sigma = f.mat) 
#}

#Solution.mat<- comb.mat
#Sim.1[[no.reps + 1]]<- Solution.mat
#notes = paste("SEED: ", SEED, " Samp.size: ", no.samp, 
#              " Rand.mat.prob: ", prob , " Tot.no.links: ",tot.no.links, 
#              " Off.diag.cov: ", off.diag.cov, " diag.cov:", diag.cov, sep = "")
#Sim.1[[no.reps + 2]]<- notes

#str(Sim.1)
#length(Sim.1)

##save(Sim.1, file = paste("SimData.Size_Cov1Var6_", no.samp,".Rdata", sep = "") )

# view solution network to recover
#sol.W<- raster(xmn = 0, xmx = K, ymn = 0, ymx = K, nrows = K, ncols = K)
#sol.W[]<- as.vector(comb.mat)

#x11()
#plot(sol.W, legend = F, main = paste("Solution to recover (",prob,")", sep = ""))

#######################
# view solution network to recover
#load("SimData.Size_Cov1Var6_1000.Rdata")
#ss<- 100
#t1<- Sim.1[[sample(1:10, 1)]]  # randomly pick 1/10 reps to show covariance
#cov.t1<- cov(t1)
#K=70
#p.cov<- raster(xmn = 0, xmx = K, ymn = 0, ymx = K, nrows = K, ncols = K)
#p.cov[]<- as.vector(cov(t1))
#x11()
#plot(p.cov, main = paste("Sample covariance for sample size ", ss,sep = ""))
#x11()
#p.cov[]<- as.vector(Sim.1[[11]])
#plot(p.cov, main = "Solution to recover (0.05)", legend = FALSE)

## **********************************************************************************************
## **********************************************************************************************
## **********************************************************************************************


            ###################################################################
            ### Begin simulation study        #################################
            ###################################################################

rm(list = ls())

library(raster)
library(ggplot2)
library(reshape2)
source("MNL_v2.r")
source("Specificity.Sensitivity.r")
source("OrderMat_1stOR2nd.R")
source("rand_W.R")
load("SimData.Size_Cov1Var6_100.Rdata")

comb.mat<- Sim.1[[11]]
comb.mat[1:10, 1:10]  # <-- solution matrix

tot.no.links = sum(rowSums(comb.mat))/2
tot.no.links

samp.size = dim(Sim.1[[1]])[1]
samp.size

K=70
no.reps = 10   # <-- for now just o get a feel of the results, run this for only a few runs

SEED = 456
set.seed(SEED)

no.runs = 10 

init.W<-  OrderMat_1stOR2nd(K=K, order='first') 
init.W[1:10, 1:10]

# Alternatively  - for different W starting value
# init.W<- rand_W(N=K)

spec.sens.mnl <- W.keep<- list()
no.links.mnl<- log.lik<- s2s.est<- mnl.spec<- mnl.sens<- init.s2s<- c()


        ## **********************
        ## Begin simulation
        ## **********************

tiMe<- proc.time()
pb <- txtProgressBar(min = 0, max = no.reps, style = 3)
for(s in 1:no.reps){
  
  init.s2s[s]<- runif(1, min = 5, max = 20)
  op.mnl<- MNL_v2(init.s2s = init.s2s[s], init.W = init.W, 
                  spatDat = Sim.1[[s]], rho.val = 0.9, K=K, no.runs = no.runs)
  #str(op.mnl)
  
  log.lik[s]<- op.mnl$final.log.lik 
  s2s.est[s]<- op.mnl$s2s[no.runs]
  W.keep[[s]]<- op.mnl$w.keep[1:K, 1:K, no.runs]
  
  spec.sens.mnl[[s]]<- Specificity.Sensitivity(W.sol = 
                                                 comb.mat, W.final = op.mnl$w.keep[,,no.runs])$twoW.tbl
  mnl.sens[s] <- Specificity.Sensitivity(W.sol = 
                                           comb.mat, W.final = op.mnl$w.keep[,,no.runs])$twoW.tbl[1,1]
  mnl.spec[s] <- Specificity.Sensitivity(W.sol = 
                                           comb.mat, W.final = op.mnl$w.keep[,,no.runs])$twoW.tbl[2,2]   
  
  no.links.mnl[s]<- sum(rowSums(op.mnl$w.keep[,,no.runs]))/2
  
  setTxtProgressBar(pb, s) 
}
close(pb)
(proc.time() - tiMe)/60

        ## **********************
        ## End simulation
        ## **********************

no.links.mnl
mnl.sens
mnl.spec
s2s.est
spec.sens.mnl[[s]]


mnl.sim<- list(samp.size = samp.size, sim.data = paste("SimData.Size_Cov1Var6_", samp.size,".Rdata", sep = ""),
               mnl.spec = mnl.spec, mnl.sens = mnl.sens, no.links.mnl = no.links.mnl, init.s2s = rep(10, no.reps),
               W.keep = W.keep, s2s.est = s2s.est, no.runs = no.runs, spec.sens.mnl = spec.sens.mnl,
               no.reps = no.reps)

str(mnl.sim)


#save(mnl.sim, file = paste("redoSimStudy_S2_SampSize", samp.size,".Rdata",sep = ""))

## *********************************
## view and compare results
## 
#sol.mnl<- diff.mnl.sol<- raster(xmn = 0, xmx = K, ymn = 0, ymx = K, nrows = K, ncols = K)
#sol.mnl[]<- as.vector(op.mnl$w.keep[1:K, 1:K, no.runs])
#diff.mnl.sol[]<- as.vector(op.mnl$w.keep[1:K, 1:K, no.runs] - comb.mat)

#x11()
#plot(sol.mnl, legend = FALSE, main = paste("MNL result", sep = ""))

#x11()
#plot(diff.mnl.sol, legend=TRUE, main = "Difference MNL and solution")

## *************************************
## Visualise results as per manuscript
## *************************************

dat.m<- data.frame(no.links = no.links.mnl, reps = 1:no.reps)
dat.m

x11()
ggplot(dat.m, aes(x = reps, y = no.links)) + 
  geom_point(size = 2) + geom_line(size = 2) + 
  theme_bw() + ggtitle(paste("MNL algorithm tot.no.link (sample size ", samp.size, ")", sep = "")) + 
  geom_hline(yintercept = tot.no.links, colour = "red") +
  ylab("Total number of links") + xlab("Replicate number") 


## **************************************
## Visualise specificity and sensitivity
## **************************************

mnl.perf<- data.frame(reps = c(1:no.reps), 
                      value = c(mnl.sens, mnl.spec),
                      Performance = c(rep("Sensitivity", no.reps), rep("Specificity", no.reps)))

x11()
ggplot(mnl.perf, aes(x = reps, y = value, group = Performance, 
                     colour = Performance)) + 
  geom_point(size = 2) + geom_line(size = 2) +
  theme_bw() + ggtitle(paste("Sensitivity & Specificity: MNL algorithm (sample size ", samp.size, ")", sep = "")) + 
  ylab("Algorithm performance") + xlab("Replicate.no") +
  ylim(c(0,1))












