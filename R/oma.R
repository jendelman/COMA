#' Optimal Mate Allocation
#'
#' Optimize the allocation for each mating
#' 
#' The data.frame \code{data} has five columns: mother, father, min, max, merit
#' 
#' @param data Input data frame (see Details)
#' @param K kinship matrix to control inbreeding
#' @param F.max maximum inbreeding coefficient (vector allowed)
#' @param solver solver for CVXR (default is "ECOS")
#' 
#' @return list containing
#' \describe{
#' \item{response}{data frame with F and merit}
#' \item{oc}{matrix of optimal contributions for each individual}
#' \item{om}{matrix of optimal allocations for each mating}
#' }
#' @import CVXR
#' @export
#' 
oma <- function(data, K, F.max, solver="ECOS") {
  
  m <- length(F.max)
  n <- nrow(data)
  H <- matrix(data$merit,nrow=1)
  mating <- paste(data$mother,data$father,sep="/")
  om <- matrix(as.numeric(NA),nrow=n,ncol=m,dimnames=list(mating,F.max))
  response <- data.frame(F=numeric(m)*NA, merit=numeric(m)*NA)
  Kvec <- matrix(apply(data[,c("mother","father")],1,function(z,K){K[z[1],z[2]]},K=K),nrow=1)
  
  for (i in 1:m) {
    x <- Variable(n)
    objective <- Maximize(H%*%x)
    constraints <- list(x >= data$min, 
                        x <= data$max,
                        sum(x)==1,
                        Kvec%*%x <= F.max[i])
    problem <- Problem(objective,constraints)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      x.opt <- result$getValue(x)
      om[,i] <- x.opt
      response$F[i] <- as.numeric(Kvec%*%x.opt)
      response$merit[i] <- as.numeric(H%*%x.opt)
    }
  }
  colnames(om) <- round(response$F,3)
  
  Fvec <- as.numeric(colnames(om))
  ix <- which(!is.na(Fvec) & !duplicated(Fvec))
  if (length(ix)==0) {
    stop("No optimal solutions")
  }
  om <- om[,ix,drop=FALSE]
  response <- response[ix,]
  
  id <- sort(union(data$mother,data$father))
  data$mother <- factor(data$mother,levels=id)
  data$father <- factor(data$father,levels=id)
  oc <- apply(om,2,function(x){
    maternal <- tapply(x,data$mother,sum,na.rm=T)
    maternal[is.na(maternal)] <- 0
    paternal <- tapply(x,data$father,sum,na.rm=T)
    paternal[is.na(paternal)] <- 0
    (maternal+paternal)/2
  })
  return(list(response=response, om=om, oc=oc))
}