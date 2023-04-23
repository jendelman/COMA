#' Optimal Contribution Selection
#'
#' Optimize the contribution of each individual to the next generation
#' 
#' The data.frame \code{data} has four mandatory columns: id, merit, min, max. An optional fifth column specifies the sex of the individual as "M" (male) or "F" (female).
#' 
#' @param data Input data frame (see Details)
#' @param K kinship matrix to control inbreeding
#' @param F.max maximum inbreeding coefficient (vector allowed)
#' @param solver solver for CVXR (default is "ECOS")
#' 
#' @return list containing
#' \describe{
#' \item{response}{data frame with F and merit}
#' \item{oc}{matrix of optimal contributions for each F.max}
#' }
#' @import CVXR
#' @export
#' 
ocs <- function(data, K, F.max, solver="ECOS") {
  
  m <- length(F.max)
  n <- nrow(data)
  H <- matrix(data$merit,nrow=1)
  oc <- matrix(as.numeric(NA),nrow=n,ncol=m,dimnames=list(data$id,F.max))
  response <- data.frame(F=numeric(m)*NA, merit=numeric(m)*NA)
  
  for (i in 1:m) {
    x <- Variable(n)
    objective <- Maximize(H%*%x)
    constraints <- list(x >= data$min, 
                        x <= data$max,
                        sum(x)==1,
                        quad_form(x,K) <= F.max[i])
    problem <- Problem(objective,constraints)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      x.opt <- result$getValue(x)
      oc[,i] <- x.opt
      response$F[i] <- as.numeric(t(x.opt)%*%K%*%x.opt)
      response$merit[i] <- as.numeric(H%*%x.opt)
    }
  }
  colnames(oc) <- round(response$F,3)
  
  Fvec <- as.numeric(colnames(oc))
  ix <- which(!is.na(Fvec) & !duplicated(Fvec))
  if (length(ix)==0) {
    stop("No optimal solutions")
  }
  oc <- oc[,ix,drop=FALSE]
  response <- response[ix,]
  
  return(list(response=response, oc=oc))
}