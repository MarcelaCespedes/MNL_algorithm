#############################################
# This function creates a random adjacency matrix
# of weights 1.
# That is a binary symmetric matrix 

rand_W<- function(N, seEd = FALSE, seEd.val = 1){
  matr<- matrix(0, N, N)
  
  lowtri <- function(matr) {
    if(seEd == TRUE){
      set.seed(seEd.val)
    }
    lt<- apply(matr, c(1, 2), function(x) sample(c(0, 1), size = 1, prob = c(0.9,0.1) ))
    for (i in 1:(nrow(lt))) {
      for (j in 1:i) {
        lt[i,j] <- abs(matr[i,j])
      }
    }
    #print(lt)
    return(lt)
  }
  
  m1<- lowtri(matr)
  #m1
  M <- m1
  for(i in 1:nrow(M)) {
    for(j in 1:i) {
      M[i,j]=M[j,i] }
  }
  #M
  #isSymmetric(M)
  rand.adj<-M
  # make off diagonals 0
  diag(rand.adj)<- 0
  return(rand.adj)  #
}