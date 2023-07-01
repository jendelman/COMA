#' Optimal Contribution Selection
#'
#' Optimize the contribution of each individual to the next generation
#' 
#' The data.frame \code{parents} can have up to five columns: id, merit, min, max, female. Min and max are real numbers between 0 and 1 specifying the minimum and maximum contribution for each parent. The "female" column is logical TRUE/FALSE. For hermaphroditic species, the last column is omitted.
#' 
#' The average inbreeding coefficient of the current generation is based on all individuals in \code{K}, which may exceed the list of individuals in \code{parents}.
#' 
#' After optimization, contributions less than \code{min.c} are set equal to zero, and the remainder are renormalized.
#' 
#' @param parents input data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param dF inbreeding rate (can be vector of values)
#' @param min.c minimum contribution
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
ocs <- function(parents, ploidy, K, dF, min.c=0, solver="ECOS") {
  
  stopifnot(c("id","merit","min","max") %in% colnames(parents))
  parents$id <- as.character(parents$id)
  if ("female" %in% colnames(parents)) {
    sexed <- TRUE
    stopifnot(!duplicated(parents$id))
    stopifnot(class(parents$female)=="logical")
    stopifnot(sum(parents$female) > 0 & sum(!parents$female) > 0)
    parents$female <- as.integer(parents$female)
  } else {
    sexed <- FALSE
  }
  stopifnot(parents$max <= 1)
  stopifnot(parents$min >= 0)
  
  m <- length(dF)
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  F.avg <- mean(as.numeric(Fi))
  n <- nrow(parents)
  tmp <- matrix(as.numeric(NA),nrow=n,ncol=m)
  colnames(tmp) <- dF
  oc <- data.frame(parents,tmp,check.names=F)
  response <- data.frame(dF.target=dF, dF1=numeric(m)*NA, 
                         merit=numeric(m)*NA)
  h.t <- matrix(parents$merit,nrow=1)
  K <- K[parents$id,parents$id]
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  
  for (i in 1:m) {
    RHS <- (ploidy-1)*(dF[i]+(1-dF[i])*F.avg)
    y <- Variable(n)
    objective <- Maximize(h.t%*%y)
    constraints <- 
      list(y <= parents$max, y >= parents$min, sum(y)==1, 
           (ploidy/2)*quad_form(y,K) + (ploidy/2-1)*Fi%*%y <= RHS)
    if (sexed)
      constraints <- c(constraints, sum(parents$female*y)==0.5)
    problem <- Problem(objective,constraints)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      y.opt <- result$getValue(y)
      y.opt[y.opt < min.c] <- 0
      y.opt <- y.opt/sum(y.opt)
      oc[,4+i] <- y.opt
      if (sexed) {
        y.f <- parents$female*y.opt
        y.m <- (1-parents$female)*y.opt
      } else {
        y.f <- y.m <- y.opt
      }
      F.t1 <- as.numeric((ploidy/2)*crossprod(y.f,K%*%y.m) + (ploidy/2-1)*Fi%*%y.opt)/(ploidy-1)
      response$dF1[i] <- round((F.t1-F.avg)/(1-F.avg),3)
      response$merit[i] <- result$value
    }
  }
  colnames(oc)[4+1:m] <- response$dF1

  return(list(response=response, oc=oc))
}