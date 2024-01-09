#' Ribbon Plot of Optimal Contributions
#' 
#' Ribbon Plot of Optimal Contributions
#' 
#' The inbreeding rates are taken from the names of \code{oc}. Each element of \code{oc} should be a data frame with columns "id" and "value".
#' 
#' @param oc named list of optimal contributions/allocations
#' @param min.c minimum contribution/allocation for plotting
#' 
#' @return ggplot2 object
#' 
#' @export
#' @importFrom tidyr pivot_wider pivot_longer
#' @import ggplot2

plot_ribbon <- function(oc, min.c=0.001) {
  stopifnot(inherits(oc,"list"))
  stopifnot(!is.null(names(oc)))
  m <- length(oc)
  dF <- as.numeric(names(oc))
  oc2 <- NULL
  for (i in 1:m) {
    oc2 <- rbind(oc2,data.frame(dF=dF[i],oc[[i]]))
  }
  ocw <- pivot_wider(oc2[oc2$value >= min.c,], names_from="dF", 
                     values_from="value")
  ocw[is.na(ocw)] <- 0
  for (i in 1:m) {
    ocw[,i+1] <- ocw[,i+1]/sum(ocw[,i+1])
  }
  plot.data <- pivot_longer(ocw,cols=2:ncol(ocw),names_to="dF",
                            values_to="value",names_prefix = "X",
                            names_transform = list(dF=as.numeric))

  if (m > 1) {
    ggplot(data=plot.data,aes(x=dF,y=value,fill=id)) + 
      scale_fill_viridis_d(name="") +
      ylab("Contribution") + xlab("Inbreeding Rate") + 
      geom_area(colour="white") +
      scale_y_continuous(breaks=seq(0,1,by=0.1)) 
      
  } else {
    ggplot(data=plot.data,aes(x=1,y=value,fill=id)) + 
      scale_fill_viridis_d(name="") +
      ylab("Contribution") + xlab("Inbreeding Rate") + 
      scale_y_continuous(breaks=seq(0,1,by=0.1)) + 
      geom_col(colour="white") + 
      scale_x_continuous(breaks=1,labels=dF[1])
  }
}