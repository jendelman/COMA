#' Optimal Mate Allocation
#'
#' Optimize the allocation for each mating
#' 
#' The data frame \code{parents} needs columns id, min, max, with optional column female (logical) for separate sexes. The data.frame \code{matings} has five columns: female, male, merit, min, max. 
#'
#' After optimization, allocations less than \code{min.c} are set equal to zero, and the remainder are renormalized. The realized inbreeding rate can exceed the specified limit after applying this threshold. It is also possible that no feasible solution exists for the specified \code{dF}. In this case, argument \code{dF.adapt} can be used to find other solutions. It is a list with two named components: step, max. The software increases the dF limit by dF.adapt$step up to the smaller of dF.adapt$max or the realized value under the original dF, in an attempt to find a solution with less inbreeding. 
#'
#' @param parents parents data frame (see Details)
#' @param matings matings data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param dF max inbreeding rate
#' @param min.c minimum allocation
#' @param dF.adapt see Details
#' @param solver solver for CVXR (default is "ECOS")
#' 
#' @return list containing
#' \describe{
#' \item{response}{named vector with realized dF and merit}
#' \item{oc}{data frame of optimal contributions for each individual}
#' \item{om}{data frame of optimal allocations for each mating}
#' }
#' @import CVXR
#' @importFrom stats model.matrix
#' @export
#' 
oma <- function(parents, matings, ploidy, K, dF, min.c=0, 
                dF.adapt=NULL, solver="ECOS") {
  
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
  stopifnot(c("female","male","merit","min","max") %in% colnames(matings))
  stopifnot(matings$max <= 1)
  stopifnot(matings$min >= 0)
  
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  F.avg <- mean(as.numeric(Fi))
  
  n <- nrow(parents)
  p <- nrow(matings)
  
  h.t <- matrix(matings$merit,nrow=1)
  parent.id <- sort(union(matings$female,matings$male))
  ix <- match(parent.id,parents$id,nomatch=0)
  if (any(ix==0))
    stop("parents data.frame missing individuals in matings")
  parents <- parents[ix,]
  
  matings$female <- factor(matings$female,levels=parent.id)
  matings$male <- factor(matings$male,levels=parent.id)

  M1 <- model.matrix(~female-1,matings)
  M2 <- model.matrix(~male-1,matings)
  
  matings$female <- as.character(matings$female)
  matings$male <- as.character(matings$male)
  
  M <- t((M1+M2)/2)
  rownames(M) <- parent.id
  
  oc <- data.frame(parents,y=numeric(n),check.names=F)
  om <- data.frame(matings,x=numeric(p),check.names=F)
  #response <- data.frame(dF.target=dF, dF1=numeric(m)*NA, dF2=numeric(m)*NA, merit=numeric(m)*NA)
  
  K <- K[parents$id,parents$id]
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  Kvec <- matrix(apply(matings[,c("female","male")],1,
                       function(z,K){K[z[1],z[2]]},K=K),nrow=1)
  
  dF1 <- dF
  dFr <- merit <- as.numeric(NA)
  dF.best <- Inf
  flag <- TRUE
  y <- Variable(n)
  x <- Variable(p)
  objective <- Maximize(h.t%*%x)
  theta1 <- ploidy*(3*ploidy/4-1)
  
  while (flag) {
    theta2 <- (ploidy/2-1)*(dF1*(ploidy-1)-ploidy/2)
    theta3 <- ploidy/2*(ploidy-1)*(1-dF1)
    RHS1 <- (ploidy-1)*(dF1+(1-dF1)*F.avg)
#    if (!is.null(dF.min))
#      LHS1 <- (ploidy-1)*(dF.min+(1-dF.min)*F.avg)
    RHS2 <- (ploidy-1)^2*dF1
    constraints <- 
      list(y <= parents$max, y >= parents$min, y==M%*%x, 
           x <= matings$max, x >= matings$min, sum(x)==1,
           (ploidy/2)*Kvec%*%x + (ploidy/2-1)*Fi%*%y <= RHS1,
           theta1*quad_form(y,K) + theta2*Fi%*%y - theta3*Kvec%*%x <= RHS2)
#    if (!is.null(dF.min)) {
#      constraints <- c(constraints, LHS1 <= (ploidy/2)*Kvec%*%x + (ploidy/2-1)*Fi%*%y)
#    }
    if (sexed)
      constraints <- c(constraints, sum(parents$female*y)==0.5)
    problem <- Problem(objective,constraints)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      x.opt <- result$getValue(x)
      ix <- which(x.opt < min.c)
      if (length(ix) > 0)
        x.opt[ix] <- 0
      if (sum(x.opt) > 0) {
        x.opt <- x.opt/sum(x.opt)
        y.opt <- M%*%x.opt
        
        F.t1 <- as.numeric((ploidy/2)*Kvec%*%x.opt + (ploidy/2-1)*Fi%*%y.opt)/(ploidy-1)
        dFr <- (F.t1-F.avg)/(1-F.avg)
        
#        F.t2 <- as.numeric(ploidy*(3*ploidy/4-1)*crossprod(y.opt,K%*%y.opt) + (ploidy/2-1)^2*Fi%*%y.opt)/(ploidy-1)^2
#        response$dF2[i] <- round((F.t2-F.t1)/(1-F.t1),3)
        if (dFr < dF.best) {
          oc[,5] <- y.opt
          om[,6] <- x.opt
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
                oc=oc[integer(0),], om=om[integer(0),]))
  } else {
    return(list(response=c(dF=dF.best, merit=merit), 
                oc=oc[which(oc[,5] > 0),,drop=FALSE],
                om=om[which(om[,6] > 0),,drop=FALSE]))
  }
}
