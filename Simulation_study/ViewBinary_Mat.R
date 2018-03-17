########################################################################
# View binary solution matrices



ViewBinary_Mat<- function(mat, R = 35, tiTle = "A title"){
  
  temp1<- melt(mat)
  #head(temp1)
  
  colnames(temp1)[1]<- "Var1"; colnames(temp1)[2]<- "Var2"
  #class(temp1$Var1)
  
  temp1$Var1<- factor(temp1$Var1)  
  temp1$Var2<- factor(temp1$Var2)  
  

  #### suppose for inference or discussion - we wish to see those regions with p-values > 0.8
  #head(temp1)
  
  hm.2 = NULL
 
  # Create plot
  hm.2 <- ggplot(temp1, aes(x= Var1, y = Var2)) + geom_tile(aes(fill = value), colour = "black") + 
    scale_fill_gradient(low = "white", high = "blue4", na.value = "black") +
    labs(x = "", y = "") + 
    ylim(rev(levels(temp1$Var2))) + # <== allows for left to right diagonal plots
    ggtitle(tiTle) +
    theme(legend.position="none") +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank())
    
  #x11()
  #hm.2

  
  return(list(hm.2 = hm.2))
  
}