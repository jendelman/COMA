#' Optimal Mate Allocation
#'
#' Optimize the allocation for each mating
#' 
#' The data frame \code{parents} needs columns id, min, max, with optional column female (logical) for separate sexes. The data.frame \code{matings} has five columns: female, male, merit, min, max. 
#'
#' After optimization, allocations less than \code{min.c} are set equal to zero, and the remainder are renormalized.
#'
#' @param parents parents data frame (see Details)
#' @param matings matings data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param dF inbreeding rate (can be vector of values)
#' @param min.c minimum allocation
#' @param solver solver for CVXR (default is "ECOS")
#' @param dF.min minimum inbreeding rate
#' 
#' @return list containing
#' \describe{
#' \item{response}{data frame with merit and dF}
#' \item{oc}{data frame of optimum contributions for each individual}
#' \item{om}{data frame of optimum allocations for each mating}
#' }
#' @import CVXR
#' @importFrom stats model.matrix
#' @export
#' 
oma <- function(parents, matings, ploidy, K, dF, min.c=0, solver="ECOS",
                dF.min=NULL) {
  
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
  
  m <- length(dF)
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
  tmp <- matrix(as.numeric(NA),nrow=p,ncol=m)
  colnames(tmp) <- dF
  om <- data.frame(matings,tmp,check.names=F)
  tmp <- matrix(as.numeric(NA),nrow=n,ncol=m)
  colnames(tmp) <- dF
  oc <- data.frame(parents,tmp,check.names=F)
  response <- data.frame(dF.target=dF, dF1=numeric(m)*NA, dF2=numeric(m)*NA, merit=numeric(m)*NA)
  
  K <- K[parents$id,parents$id]
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  Kvec <- matrix(apply(matings[,c("female","male")],1,
                       function(z,K){K[z[1],z[2]]},K=K),nrow=1)
  i=1
  for (i in 1:m) {
    theta1 <- ploidy*(3*ploidy/4-1)
    theta2 <- (ploidy/2-1)*(dF[i]*(ploidy-1)-ploidy/2)
    theta3 <- ploidy/2*(ploidy-1)*(1-dF[i])
    RHS1 <- (ploidy-1)*(dF[i]+(1-dF[i])*F.avg)
    if (!is.null(dF.min))
      LHS1 <- (ploidy-1)*(dF.min+(1-dF.min)*F.avg)
    RHS2 <- (ploidy-1)^2*dF[i]
    y <- Variable(n)
    x <- Variable(p)
    objective <- Maximize(h.t%*%x)
    constraints <- 
      list(y <= parents$max, y >= parents$min, y==M%*%x, 
           x <= matings$max, x >= matings$min, sum(x)==1,
           (ploidy/2)*Kvec%*%x + (ploidy/2-1)*Fi%*%y <= RHS1,
           theta1*quad_form(y,K) + theta2*Fi%*%y - theta3*Kvec%*%x <= RHS2)
    if (!is.null(dF.min)) {
      constraints <- c(constraints, LHS1 <= (ploidy/2)*Kvec%*%x + (ploidy/2-1)*Fi%*%y)
    }
    if (sexed)
      constraints <- c(constraints, sum(parents$female*y)==0.5)
    problem <- Problem(objective,constraints)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      x.opt <- result$getValue(x)
      x.opt[x.opt < min.c] <- 0
      x.opt <- x.opt/sum(x.opt)
      y.opt <- M%*%x.opt
      
      oc[,4+i] <- y.opt
      om[,5+i] <- x.opt
      
      F.t1 <- as.numeric((ploidy/2)*Kvec%*%x.opt + (ploidy/2-1)*Fi%*%y.opt)/(ploidy-1)
      response$dF1[i] <- round((F.t1-F.avg)/(1-F.avg),3)
      
      F.t2 <- as.numeric(ploidy*(3*ploidy/4-1)*crossprod(y.opt,K%*%y.opt) + (ploidy/2-1)^2*Fi%*%y.opt)/(ploidy-1)^2
      response$dF2[i] <- round((F.t2-F.t1)/(1-F.t1),3)
      
      response$merit[i] <- result$value
    }
  }
  colnames(oc)[4+1:m] <- round(response$dF1,3)
  colnames(om)[5+1:m] <- round(response$dF1,3)
  max.c <- apply(om[,5+1:m,drop=FALSE],1,max)
  ix <- which(max.c >= min.c)
  om <- om[ix,,drop=FALSE]
  max.c <- apply(oc[,4+1:m,drop=FALSE],1,max)
  ix <- which(max.c >= min.c)
  oc <- oc[ix,,drop=FALSE]
  
  return(list(response=response, om=om, oc=oc))
}
