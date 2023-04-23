#' Plot Optimal Contributions
#' 
#' Plot Optimal Contributions
#' 
#' @param data matrix of optimal contributions
#' @param min.c minimum contribution for plotting
#' 
#' @return ggplot2 object
#' 
#' @export
#' @importFrom tidyr pivot_longer
#' @import ggplot2

plot_opt <- function(data, min.c) {
  Fvec <- as.numeric(colnames(data))
  ix <- which(!is.na(Fvec) & !duplicated(Fvec))
  if (length(ix)==0) {
    stop("No optimal solutions")
  }
  data <- data[,ix,drop=FALSE]
  Fvec <- as.numeric(colnames(data))
  m <- ncol(data)
  colnames(data) <- 1:m
  
  max.c <- apply(data,1,max)
  ix <- which(max.c >= min.c)
  data2 <- data[ix,,drop=FALSE]
  data2 <- data.frame(id=rownames(data2),data2,check.names=F,row.names = NULL)
  
  plot.data <- pivot_longer(data2,cols=2:ncol(data2),names_to="F",
                            values_to="contrib",names_prefix = "X",
                            names_transform = list(F=as.numeric))
  plot.data$F <- Fvec[match(plot.data$F,1:m)]

  if (m > 1) {
    ggplot(data=plot.data,aes(x=F,y=contrib,fill=id)) + 
      scale_fill_viridis_d(name="") +
      ylab("Contribution") + xlab("Inbreeding") + 
      scale_y_continuous(breaks=seq(0,1,by=0.1)) + geom_area(colour="white") 
      
  } else {
    ggplot(data=plot.data,aes(x=1,y=contrib,fill=id)) + 
      scale_fill_viridis_d(name="") +
      ylab("Contribution") + xlab("Inbreeding") + 
      scale_y_continuous(breaks=seq(0,1,by=0.1)) + geom_col(colour="white") + 
      scale_x_continuous(breaks=1,labels=Fvec[1])
  }
}