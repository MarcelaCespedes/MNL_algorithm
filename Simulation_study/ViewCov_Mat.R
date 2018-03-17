########################################################################
# Modified from ViewMatrices.r
# this is for non-binary covariance matrices
#
# For MNL paper
# Tuesday 13.3.2018



ViewCov_Mat<- function(mat, R = 35, tiTle = "A title", range = c(-5, 5), legend = TRUE){
  
  temp1<- melt(mat)
  #head(temp1)
  
  colnames(temp1)[1]<- "Var1"; colnames(temp1)[2]<- "Var2"
  colnames(temp1)[3]<- "Covariance"
  #class(temp1$Var1)
  
  temp1$Var1<- factor(temp1$Var1)  
  temp1$Var2<- factor(temp1$Var2)  
  
  # ************************************ first heap map plot  ********************
  # tweek it a little more
  base_size <- 9
  
  # in ggplot opts() is depricated
  #b = c(-1, 0.4, 0.6, 1)
  
  if(legend == TRUE){
    hm.1<- ggplot(temp1, aes(x= Var1, y = Var2)) + geom_tile(aes(fill = Covariance), colour = "white") + 
      scale_fill_gradientn(limits = range,colours = c("white","yellow", "red","green","blue", "black")) + 
      labs(x = "", y = "") + 
      #scale_x_discrete(expand = c(0, 0)) +
      ylim(rev(levels(temp1$Var2))) + # <== allows for left to right diagonal plots
      theme(axis.ticks = element_blank(),
            axis.text = element_blank()) +
      ggtitle(tiTle)
    
  }else{
    hm.1<- ggplot(temp1, aes(x= Var1, y = Var2)) + geom_tile(aes(fill = Covariance), colour = "white") + 
      scale_fill_gradientn(limits = range,colours = c("white","yellow", "red","green","blue", "black")) + 
      labs(x = "", y = "") + 
      ylim(rev(levels(temp1$Var2))) + # <== allows for left to right diagonal plots
      theme(axis.ticks = element_blank(),
            axis.text = element_blank()) +
      ggtitle(tiTle) +
      theme(legend.position="none") 
  }
 
  return(list(hm.1 = hm.1))
  
}