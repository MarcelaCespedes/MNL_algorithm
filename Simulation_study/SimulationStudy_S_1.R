######################################################
# MNL Simulation study for S_1 matrix/ data
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
#prob = 0.1
#mat.2<- SparseRandom_Mat(K=K, prob = prob, SEED = 123)

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

#off.diag.cov = 0.5
#f.mat<- off.diag.cov*comb.mat  # when this constant is small --> poor MNL results, better with high values 
#diag.cov <- 3
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

##save(Sim.1, file = paste("SimData.Size_Cov0.5Var3_", no.samp,".Rdata", sep = "") )

# view solution network to recover
#sol.W<- raster(xmn = 0, xmx = K, ymn = 0, ymx = K, nrows = K, ncols = K)
#sol.W[]<- as.vector(comb.mat)
#x11()
#plot(sol.W, legend = F, main = paste("Solution to recover (prob = 0.1)", sep = ""))

#ss<- 100
#load(paste("SimData.Size_Cov0.5Var3_", ss, ".Rdata", sep = "") )
#t1<- Sim.1[[sample(1:10, 1)]]  # randomly pick 1/10 reps to show covariance
#cov.t1<- cov(t1)
#p.cov<- raster(xmn = 0, xmx = K, ymn = 0, ymx = K, nrows = K, ncols = K)
#p.cov[]<- as.vector(cov(t1))
#x11()
#plot(p.cov, main = paste("Sample covariance for sample size ", ss,sep = ""))


## **********************************************************************************************
## **********************************************************************************************
## **********************************************************************************************


                ###################################################################
                ### Begin simulation study        #################################
                ###################################################################


rm(list = ls())

## interestingly it takes 6,9,15,xx minutes to run the algorithm for the different sample sizes

library(raster)
library(ggplot2)
library(reshape2)
source("MNL_v2.r")
source("Specificity.Sensitivity.r")
source("OrderMat_1stOR2nd.R")
load("SimData.Size_Cov0.5Var3_100.Rdata")

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

no.runs = 10 # noticeable improvement when this was doubled to 20

other.init<-  OrderMat_1stOR2nd(K=K, order='first') #matrix(0, K, K)
other.init[1:10, 1:10]
spec.sens.mnl <- W.keep<- list()
no.links.mnl<- log.lik<- s2s.est<- mnl.spec<- mnl.sens<- init.s2s<- c()


## Is there a difference with different starting values for s2s and W??
## should I keep track of the covariance matrix and compare with the sample covariance
## of the data?

##
tiMe<- proc.time()
pb <- txtProgressBar(min = 0, max = no.reps, style = 3)
for(s in 1:no.reps){
  
  init.s2s[s]<- runif(1, min = 5, max = 20)
  op.mnl<- MNL_v2(init.s2s = init.s2s[s], init.W = other.init,#mat.1st.ord, 
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


mnl.sim<- list(samp.size = samp.size, sim.data = paste("SimData.Size_Cov0.5Var3_", samp.size,".Rdata", sep = ""),
               mnl.spec = mnl.spec, mnl.sens = mnl.sens, no.links.mnl = no.links.mnl, init.s2s = init.s2s,
               W.keep = W.keep, s2s.est = s2s.est, no.runs = no.runs, spec.sens.mnl = spec.sens.mnl,
               no.reps = no.reps)

str(mnl.sim)


save(mnl.sim, file = paste("SimStudy_S1_SampSize", samp.size,".Rdata",sep = ""))


# Ugh! this takes about half an hour to run on a single core!!!
## quick look at a single matrix

# view solution network to recover
sol.mnl<- diff.mnl.sol<- raster(xmn = 0, xmx = K, ymn = 0, ymx = K, nrows = K, ncols = K)
sol.mnl[]<- as.vector(op.mnl$w.keep[1:K, 1:K, no.runs])
diff.mnl.sol[]<- as.vector(op.mnl$w.keep[1:K, 1:K, no.runs] - comb.mat)

x11()
plot(sol.mnl, legend = FALSE, main = paste("MNL result", sep = ""))

x11()
plot(diff.mnl.sol, legend=TRUE, main = "Difference MNL and solution")
##

# init W at zero matrix, rho.val = 0.9, no.runs=20
#          Pred.pos Pred.neg
#Sol.pos     0.28     0.72
#Sol.neg     0.01     0.99

##
## I tried to make it run in parallel and wasted bout 1.5 hours on this - turns out there is a problem
## with running MNL_v2 because it calls an C++ function dmvnrm_arma() and the workers have a problem
## compiling this. Running the code from here
## https://stackoverflow.com/questions/18245193/doparallel-issue-with-inline-function-on-windows-7-works-on-linux?rq=1
## shows a replication of the vague error I had above "Error in f1() : task 1 failed - "NULL value passed as symbol address" "
## and their work around was to include the C++ in the function MNL_v2
##


no.links.mnl
mnl.sens
mnl.spec
s2s.est
# it seems that this is not performing its best - and one potential reason why
# is because I have to allow the estimate of s2s to be "unconstrained"
# for values of s2s too low (eg max(s2s)~ 2) then it underestimates the number of links

spec.sens.mnl[[s]]


#> no.links.mnl on a not-sparse network
#[1] 203 201 204 199 199 199 200 197 198 203
#> s2s.est
#[1] 20 20 20 20 20 20 20 20 20 20 # <--- limit for this was set to between 0.01 to 100


### make plot for paper
dat.m<- data.frame(no.links = no.links.mnl, reps = 1:no.reps)
dat.m

x11()
ggplot(dat.m, aes(x = reps, y = no.links)) + 
  geom_point(size = 2) + geom_line(size = 2) + 
  theme_bw() + ggtitle(paste("MNL algorithm tot.no.link (sample size ", samp.size, ")", sep = "")) + 
  geom_hline(yintercept = tot.no.links, colour = "red") +
  ylab("Total number of links") + xlab("Replicate number") 


### do similar plots for specificity and sensitivity
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