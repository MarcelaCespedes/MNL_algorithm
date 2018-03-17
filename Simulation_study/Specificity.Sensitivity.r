##########################################################
# Function to define in matrices - which are the true positive and true negative
# given the solution and the final matrix

# Maybe if needed - I can add the 

# Sensitivity = true positive
# specificity = true negative

# to return the % of true positive, and % or true negative

Specificity.Sensitivity<- function(W.sol, W.final){
  
  K = dim(W.sol)[1]
  no.off.diag = K*(K-1)/2

  c1 <- c2<- c3<- c4<- 1
  sol.rn<- est.rn<- c()  # locate links
  sol.nl<- est.nl<- c()  # locate where tha absence of links are
  
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      
      if(W.final[i,j] == 1){
        est.rn[c1]<- paste("R", i,".", j, sep = "")  #<== list of estimated links
        c1 = c1 + 1
      }
      
      if(W.final[i,j] == 0){
        est.nl[c2]<- paste("R", i,".", j, sep = "")  #<== list of estimated no links
        c2 = c2 + 1
      }
  #-------------------------------------      
      
      if(W.sol[i,j] == 1){
        sol.rn[c3]<- paste("R", i,".", j, sep = "")  #<== list of correct links
        c3 = c3 + 1
      }
      
      if(W.sol[i,j] == 0){
        sol.nl[c4]<- paste("R", i,".", j, sep = "")  #<== list of correct no links
        c4 = c4 + 1
      }
      
    }
  }  
  
  no.sol.links<- length(sol.rn)
  no.sol.links
  
  absent.sol.links<- length(sol.nl)
  absent.sol.links  
  
  
  ## True positive
  TP<- length(intersect(sol.rn, est.rn))
  
  ## True negative
  TN<- length(intersect(sol.nl, est.nl))
  
  ## False positive - regions links identified as neighbours when they really aren't
  FP <- length(intersect(est.rn, sol.nl) )
  
  ## False negative - links which are thought to NOT be neighbours when they really are
  FN <- length(intersect(est.nl, sol.rn) )
  
  
  #### -- off wikipedia
  Sensitivity = TP/(TP + FN)
  Sensitivity
  
  Specificity = TN/(TN + FP)
  Specificity
  
  ### two way confusion table
  res<- matrix(c(TP, FP, FN, TN), 2,2)
  colnames(res)<- c("Pred.pos", "Pred.neg")
  rownames(res)<- c("Sol.pos", "Sol.neg")
  
  res<- as.table(res)
  res
  
  res1<- res
  res1[1,] = res1[1,]/no.sol.links
  res1[2,] = res1[2,]/absent.sol.links
  round(res1,2)
  ## OR get it by directly comparing it with the solution
  ## True positive
  #Sens<- length(intersect(sol.rn, est.rn))/ no.sol.links
  #Sens
  ## True negative
  #Spec<- length(intersect(sol.nl, est.nl))/ absent.sol.links
  #Spec
  ########################################
  return(list(Specifi = Specificity, Sentiv = Sensitivity, TP=TP, TN=TN, FP=FP, FN=FN, twoW.tbl = round(res1,2)))
}