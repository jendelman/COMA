#' Optimal Contribution Selection
#'
#' Optimize the contribution of each individual to the next generation
#' 
#' The data.frame \code{parents} can have up to five columns: id, merit, min, max, female. Min and max are real numbers between 0 and 1 specifying the minimum and maximum contribution for each parent. The "female" column is logical TRUE/FALSE. For hermaphroditic species, the last column is omitted.
#' 
#' The average inbreeding coefficient of the current generation is based on all individuals in \code{K}, which may exceed the list of individuals in \code{parents}.
#' 
#' After optimization, contributions less than \code{min.c} are set equal to zero, and the remainder are renormalized. The realized inbreeding rate can exceed the specified limit after applying this threshold. It is also possible that no feasible solution exists for the specified \code{dF}. In this case, argument \code{dF.adapt} can be used to find other solutions. It is a list with two named components: step, max. The software increases the dF limit by dF.adapt$step up to the smaller of dF.adapt$max or the realized value under the original dF, in an attempt to find a solution with less inbreeding. 
#' 
#' @param parents input data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param dF maximum inbreeding rate 
#' @param min.c minimum contribution
#' @param dF.adapt see Details
#' @param solver solver for CVXR (default is "ECOS")
#' 
#' @return list containing
#' \describe{
#' \item{response}{named vector with realized dF and merit}
#' \item{oc}{data frame of optimal contributions}
#' }
#' @import CVXR
#' @export
#' 
ocs <- function(parents, ploidy, K, dF, min.c=0, dF.adapt=NULL, solver="ECOS") {
  
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
  
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  F.avg <- mean(as.numeric(Fi))
  n <- nrow(parents)
  
  oc <- data.frame(parents,y=numeric(n),check.names=F)
  
  h.t <- matrix(parents$merit,nrow=1)
  K <- K[parents$id,parents$id]
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  
  y <- Variable(n)
  objective <- Maximize(h.t%*%y)
  
  dF1 <- dF
  dFr <- merit <- as.numeric(NA)
  dF.best <- Inf
  flag <- TRUE
  
  while (flag) { 
    RHS <- (ploidy-1)*(dF1+(1-dF1)*F.avg)
    constraints <- 
        list(y <= parents$max, y >= parents$min, sum(y)==1, 
             (ploidy/2)*quad_form(y,K) + (ploidy/2-1)*Fi%*%y <= RHS)
    if (sexed)
      constraints <- c(constraints, sum(parents$female*y)==0.5)
      
    problem <- Problem(objective,constraints)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      y.opt <- result$getValue(y)
      ix <- which(y.opt < min.c)
      if (length(ix) > 0)
        y.opt[ix] <- 0
      if (sum(y.opt) > 0) {
        y.opt <- y.opt/sum(y.opt)
        
        if (sexed) {
          y.f <- parents$female*y.opt
          y.m <- (1-parents$female)*y.opt
        } else {
          y.f <- y.m <- y.opt
        }
        
        F.t1 <- as.numeric((ploidy/2)*crossprod(y.f,K%*%y.m) + 
                             (ploidy/2-1)*Fi%*%y.opt)/(ploidy-1)
        dFr <- (F.t1-F.avg)/(1-F.avg)
        
        if (dFr < dF.best) {
          oc[,5] <- y.opt
          merit <- result$value
          dF.best <- dFr
        }
      } 
    } 
    
    if (is.null(dF.adapt)) {
      flag <- FALSE
    } else {
      if ((!is.na(merit) & dFr <= dF1) | (dF1 > min(dF.adapt$max,dF.best+dF.adapt$step))) {
        flag <- FALSE
      } else {
        dF1 <- dF1 + dF.adapt$step
      }
    }
  }

  if (is.na(merit)) {
    return(list(response=c(dF=as.numeric(NA), merit=as.numeric(NA)), 
                oc=oc[integer(0),]))
  } else {
    return(list(response=c(dF=dF.best, merit=merit), 
                oc=oc[which(oc[,5] > 0),,drop=FALSE]))
  }
}