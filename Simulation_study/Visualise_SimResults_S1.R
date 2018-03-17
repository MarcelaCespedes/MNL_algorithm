################################################################
## Post simulation study - here I redo the plots for the paper
## which combines and compares the results for the MNL vs gLASSO and
# MNL vs PPC
#
# Tuesday 13.3.2018
#
# this is for S1 simulation results ONLY!!!


library(ggplot2)
library(reshape2)

source("ViewBinary_Mat.R")
source("multiplot.r")
source("ViewCov_Mat.R")

####
####
####  Binary matrices to recover
load("SimData.Size_Cov1Var6_100.Rdata")  # this is sparser network, high covariance values
p.0.05<- Sim.1[[11]]

rm(Sim.1) 

load("SimData.Size_Cov0.5Var3_100.Rdata")  # this is sparser network, high covariance values
p.0.1<- Sim.1[[11]]

K=70

plot0.05<- plot0.1<- raster(xmn = 0, xmx = K, ymn = 0, ymx = K, nrows = K, ncols = K)

###
### raster
###

#x11()
tiff("SolutionMatrices.tiff", units = "in", width = 5, height = 3.1, res = 300)
par(mfrow = c(1,2))
plot0.1[]<- as.vector(p.0.1)
plot0.05[]<- as.vector(p.0.05)

plot(plot0.1, main = "S_1 (prob = 0.1)", legend = FALSE, xlab = "", ylab = "")
plot(plot0.05, main = "S_2 (prob = 0.05)", legend = FALSE,xlab = "", ylab = "")
dev.off()

###
### ggplot function
###

op.0.1<- ViewBinary_Mat(mat = p.0.1, R=70, tiTle = expression(paste(S[1]," matrix (off-diagonal prob = 0.1)", sep = "")))
op.0.05<- ViewBinary_Mat(mat = p.0.05, R=70, tiTle = expression(paste(S[2]," matrix (off-diagonal prob = 0.05)", sep = "")))

#x11()
tiff("SolutionMatrices.tiff", units = "in", width = 8, height = 5, res = 300)
#op.0.1$hm.2
#op.0.05$hm.2
multiplot(op.0.1$hm.2, op.0.05$hm.2, cols = 2)
dev.off()

##
ggsave(filename = "SolutionMatrices.png", plot = multiplot(op.0.1$hm.2, op.0.05$hm.2, cols = 2),
       width = 8, height = 5, units = "in", dpi = 300)

######################################################################################
######################################################################################

####
####  Less sparse matrix S_1
####

            ###############################
            ### Pearson's correlation   ###
            ###############################

rm(list=ls())

source("multiplot.r")
source("Specificity.Sensitivity.r")
library(reshape2)
library(ggplot2)

K=70
no.reps = 10
SEED = 456
set.seed(SEED)
thresh.range<- seq(from = 0.01, to= 0.2, by = 0.01)
thresh.range
plot.list<- list()

samp.size<- c(100, 250, 500, 1000)


for(i in 1:4){
    
    load(paste("SimData.Size_Cov0.5Var3_",samp.size[i],".Rdata", sep = ""))
    comb.mat<- Sim.1[[11]]

    pearson.spec.sens<- list()
    count = 1
    pearson.links<- pearson.spec<- pearson.sens<- matrix(0, nrow = no.reps, ncol = length(thresh.range))
    
    ##
    #pb <- txtProgressBar(min = 0, max = no.reps*length(thresh.range), style = 3)
    for(s in 1:no.reps){
      for(t in 1:length(thresh.range)){
        
        dat.temp<- Sim.1[[s]]
        corr.mat<- matrix(0, K, K)
        
        for(r in 1:K){
          for(c in 1:K){
            # consider the absolute value of the correlations
            corr.mat[r, c]<- abs(cor(dat.temp[, r], dat.temp[, c])) 
          }
        }
        
        p.mat<- ifelse(corr.mat > thresh.range[t], 1, 0)
        diag(p.mat) <- 0
        p.mat[1:10, 1:10]
        
        spec.sens.pear<-Specificity.Sensitivity(W.sol = comb.mat, W.final = p.mat)
        pearson.spec.sens[[count]] <- spec.sens.pear$twoW.tbl
        
        pearson.sens[s,t]<- spec.sens.pear$twoW.tbl[1,1]
        pearson.spec[s,t]<- spec.sens.pear$twoW.tbl[2,2]
        
        #setTxtProgressBar(pb, count)
        count = count + 1
      }
    }
    #close(pb)
    ##
    
    ### do similar plots for specificity and sensitivity
    p.tot<- rbind(data.frame(pearson.sens), data.frame(pearson.spec))
    p.tot$Performance<- c(rep("Sensitivity", no.reps), rep("Specificity", no.reps))
    colnames(p.tot)<- c(thresh.range, "Performance")
    
    p.tot.1<- melt(p.tot, id.vars = "Performance")
    colnames(p.tot.1)<- c("Performance", "Threshold.Range", "value")
    p.tot.1$Performance<- factor(p.tot.1$Performance)
    
    ##########################
    #x11()  # ----- make plot
    p1<-ggplot(p.tot.1, aes(x = Threshold.Range, y = value, group = Threshold.Range, 
                        colour = Performance)) + 
      geom_point() +
      theme_bw() + 
      #ggtitle(paste("Sensitivity & Specificity: Pearson pairwise correlation (samp. size ", samp.size[i],")", sep = "")) + 
      ggtitle(paste("Sample size: ", samp.size[i], sep = "")) +
      ylab("Algorithm performance") + xlab("Threshold range") + theme(legend.position="none") 
    #p1
    
    load(paste("SimStudy_S1_SampSize",samp.size[i],".Rdata", sep = ""))
    #str(mnl.sim)
    mnl.dat<- data.frame(Threshold.Range = thresh.range,
                         Sensitivity = rep(mean(mnl.sim$mnl.sens), length(thresh.range)),
                         Specificity = rep(mean(mnl.sim$mnl.spec, length(thresh.range))))
    mnl.dat<- melt(mnl.dat, id.vars = "Threshold.Range")
    mnl.dat$Threshold.Range<- factor(mnl.dat$Threshold.Range)
    colnames(mnl.dat)[2]<- "Performance"
    
    p1<- p1 + geom_line(data=mnl.dat, aes(x = Threshold.Range, y = value, group = Performance,colour = Performance), size = 1) +
      theme(legend.position="none") 
    #p1
   plot.list[[i]]<- p1
   
   ##
   print(i)
}

## combine ggplot objects
length(plot.list)

x11()
multiplot(plotlist = plot.list, cols = 2)

##
#save(plot.list, file = "SimStudy_S1_allPearsonCovar.Rdata")


## *********************************************************

          ###########################
          ###    glasso           ###
          ###########################

rm(list = ls())

source("multiplot.r")
library(reshape2)
library(glasso)
source("Specificity.Sensitivity.r")

K=70
no.reps = 10
SEED = 456
set.seed(SEED)
sparsity<- seq(from = 0.1, to = 1, by = 0.05)  # this range varies 
sparsity
samp.size<- c(100,250,500,1000)
plot.list<- list()

for(i in 1:4){
    load(paste("SimData.Size_Cov0.5Var3_", samp.size[i],".Rdata", sep = ""))
    comb.mat<- Sim.1[[11]]
    count = 1
    spec.sens.glass<- list()
    glass.links<- g.spec<- g.sens<- matrix(0, nrow = no.reps, ncol = length(sparsity))
    
    ##
    for(s in 1:no.reps){
      temp.covar<- cov(Sim.1[[s]])
      
      for(r in 1:length(sparsity)){
        op.glass<- glasso(temp.covar, rho = sparsity[r])
        #str(op.glass)
        
        mat1 = ifelse(round(op.glass$w,9) == 0, 0, 1)  # <-- !! rounding error
        diag(mat1) = 0
        #mat1[1:10, 1:10]
        
        spec.sens.glass[[count]] <- Specificity.Sensitivity(W.sol = 
                                                              comb.mat, W.final = mat1)$twoW.tbl
        count = count + 1
        
        glass.links[s, r]<- sum(rowSums(mat1))/2
        g.sens[s,r]<- Specificity.Sensitivity(W.sol = 
                                                comb.mat, W.final = mat1)$twoW.tbl[1,1]
        g.spec[s,r]<- Specificity.Sensitivity(W.sol = 
                                                comb.mat, W.final = mat1)$twoW.tbl[2,2]
      }
    }
    ##
    
    ### do similar plots for specificity and sensitivity
    g.tot<- rbind(data.frame(g.sens), data.frame(g.spec))
    g.tot$Performance<- c(rep("Sensitivity", no.reps), rep("Specificity", no.reps))
    colnames(g.tot)<- c(sparsity, "Performance")
    
    g.tot.1<- melt(g.tot, id.vars = "Performance")
    colnames(g.tot.1)<- c("Performance", "Sparsity", "value")
    
    g.tot.1$Performance<- factor(g.tot.1$Performance)
    #head(g.tot.1)
    
    p1<- ggplot(g.tot.1, aes(x = Sparsity, y = value, group = Sparsity, 
                        colour = Performance)) + 
      geom_point() +
      theme_bw() + 
      ggtitle(paste("Sample size: ", samp.size[i], sep = "")) +
      #ggtitle(paste("Sensitivity & Specificity: gLASSO (sample size ", samp.size, ")", sep = "")) + 
      ylab("Algorithm performance") + xlab("Sparsity") + theme(legend.position="none") 
    
    ##
    load(paste("SimStudy_S1_SampSize",samp.size[i],".Rdata", sep = ""))
    #str(mnl.sim)
    mnl.dat<- data.frame(Sparsity= sparsity,
                         Sensitivity = rep(mean(mnl.sim$mnl.sens), length(sparsity)),
                         Specificity = rep(mean(mnl.sim$mnl.spec, length(sparsity))))
    mnl.dat<- melt(mnl.dat, id.vars = "Sparsity")
    mnl.dat$Sparsity<- factor(mnl.dat$Sparsity)
    colnames(mnl.dat)[2]<- "Performance"
    
    p1<- p1 + geom_line(data=mnl.dat, aes(x = Sparsity, y = value, group = Performance,colour = Performance), size = 1) +
      theme(legend.position="none") 
    #p1
    plot.list[[i]]<- p1
    
    ##
    print(i)

}

## combine ggplot objects
length(plot.list)

x11()
multiplot(plotlist = plot.list, cols = 2)

##
save(plot.list, file = "SimStudy_S1_allgLASSO.Rdata")

###################################################################################
##
##  Combine results above
##

rm(list = ls())

library(ggpubr) # see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
source("multiplot.r")

load("SimStudy_S1_allPearsonCovar.Rdata")
p.plots<- plot.list

rm(plot.list)

load("SimStudy_S1_allgLASSO.Rdata")
g.plots<- plot.list

tot<- c(p.plots, g.plots) #append(p.plots, g.plots)
length(tot)

x11()
multiplot(plotlist = tot, cols = 2)  # <--- NO labels


tot2<- list()
tot2[[1]]<- p.plots[[1]]
tot2[[2]]<- g.plots[[1]]
tot2[[3]]<- p.plots[[2]]
tot2[[4]]<- g.plots[[2]]

tot2[[5]]<- p.plots[[3]]
tot2[[6]]<- g.plots[[3]]
tot2[[7]]<- p.plots[[4]]
tot2[[8]]<- g.plots[[4]]


x11()
p2<- ggarrange(plotlist = tot2, ncol =2, nrow = 4, labels = c("A","E", "B", "F", "C", "G", "D", "H"))
p2


#tiff("S1_allSimResults.tiff", units = "in", width = 13, height =15, res = 300)
#p2
#dev.off()

#ggsave(filename = "S1_allSimResults.png", plot = p2,
#       width = 13, height = 15, units = "in", dpi = 300)

#######################################################################################
#######################################################################################

####
#### covariance matrices


##
load(paste("SimData.Size_Cov0.5Var3_100.Rdata", sep = "") )
t1<- Sim.1[[sample(1:10, 1)]]  # randomly pick 1/10 reps to show covariance
cov.t1<- cov(t1)

s100<- ViewCov_Mat(mat = cov.t1, R = 70, tiTle = "Sample size 100", range = c(-1,4), legend = FALSE)

#x11()
#s100$hm.1

##
load(paste("SimData.Size_Cov0.5Var3_250.Rdata", sep = "") )
t1<- Sim.1[[sample(1:10, 1)]]  # randomly pick 1/10 reps to show covariance
cov.t1<- cov(t1)

s250<- ViewCov_Mat(mat = cov.t1, R = 70, tiTle = "Sample size 250", range = c(-1,4), legend=FALSE)

#x11()
#s250$hm.1

##
load(paste("SimData.Size_Cov0.5Var3_500.Rdata", sep = "") )
t1<- Sim.1[[sample(1:10, 1)]]  # randomly pick 1/10 reps to show covariance
cov.t1<- cov(t1)

s500<- ViewCov_Mat(mat = cov.t1, R = 70, tiTle = "Sample size 500", range = c(-1,4), legend = TRUE)

#x11()
#s500$hm.1


##
load(paste("SimData.Size_Cov0.5Var3_1000.Rdata", sep = "") )
t1<- Sim.1[[sample(1:10, 1)]]  # randomly pick 1/10 reps to show covariance
cov.t1<- cov(t1)

s1000<- ViewCov_Mat(mat = cov.t1, R = 70, tiTle = "Sample size 1000", range = c(-1,4), legend=FALSE)

#x11()
#s1000$hm.1

x11()
multiplot(plotlist = list(s100$hm.1, s250$hm.1,s500$hm.1,s1000$hm.1), cols=2)

#tiff("S1_sampleCovriances.tiff", units = "in", width = 8, height =8, res = 300)
#  multiplot(plotlist = list(s100$hm.1, s250$hm.1,s500$hm.1,s1000$hm.1), cols=2)
#dev.off()

#ggsave(filename = "S1_sampleCovriances.png", plot = multiplot(plotlist = list(s100$hm.1, s250$hm.1,s500$hm.1,s1000$hm.1), cols=2),
#       width = 8, height = 8, units = "in", dpi = 300)
