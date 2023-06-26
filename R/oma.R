#' Optimal Mate Allocation
#'
#' Optimize the allocation for each mating
#' 
#' The data frame \code{parents} needs columns id, min, max, with optional column female (logical) for separate sexes. The data.frame \code{matings} has five columns: female, male, merit, min, max. 
#'
#' @param parents parents data frame (see Details)
#' @param matings matings data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param dF inbreeding rate (can be vector of values)
#' @param solver solver for CVXR (default is "ECOS")
#' 
#' @return list containing
#' \describe{
#' \item{response}{data frame with F and merit}
#' \item{oc}{matrix of optimal contributions for each individual}
#' \item{om}{matrix of optimal allocations for each mating}
#' }
#' @import CVXR
#' @importFrom stats model.matrix
#' @export
#' 
oma <- function(parents, matings, ploidy, K, dF, solver="ECOS") {
  
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
  theta1 <- (ploidy/2-1)/(ploidy-1)
  theta2 <- (ploidy-1)/(ploidy/2)
  F.avg <- (ploidy*mean(diag(K))-1)/(ploidy-1)
  n <- nrow(parents)
  p <- nrow(matings)
  
  h.t <- matrix(matings$merit,nrow=1)
  matings$id <- paste(matings$female,matings$male,sep="/")
  parent.id <- sort(union(matings$female,matings$male))
  ix <- match(parent.id,parents$id,nomatch=0)
  if (any(ix==0))
    stop("parents data.frame missing individuals in matings")
  parents <- parents[ix,]
  
  matings$female <- factor(matings$female,levels=parent.id)
  matings$male <- factor(matings$male,levels=parent.id)

  M1 <- model.matrix(~female-1,matings)
  M2 <- model.matrix(~male-1,matings)
  M <- t((M1+M2)/2)
  rownames(M) <- parent.id
#  K2 <- crossprod(Z,K%*%Z)
  
  om <- matrix(as.numeric(NA),nrow=p,ncol=m,dimnames=list(matings$id,dF))
  oc <- matrix(as.numeric(NA),nrow=n,ncol=m,dimnames=list(parents$id,dF))
  response <- data.frame(dF.target=dF, dF.realized=numeric(m)*NA, merit=numeric(m)*NA)
  
  K <- K[parents$id,parents$id]
  Kvec <- matrix(apply(matings[,c("female","male")],1,
                       function(z,K){K[z[1],z[2]]},K=K),nrow=1)
  
  for (i in 1:m) {
    b2 <- theta2*(2-dF[i]-theta1)*dF[i] + theta2*(1-dF[i])*(1-dF[i]-theta1)*F.avg
    b1 <- theta2*dF[i] + theta2*(1-dF[i])*F.avg + 2*theta1/ploidy
    y <- Variable(n)
    x <- Variable(p)
    objective <- Maximize(h.t%*%x)
    constraints <- 
      list(y <= parents$max, y >= parents$min, y==M%*%x, 
           x <= matings$max, x >= matings$min, sum(x)==1,
           Kvec%*%x + 2*theta1*t(y)%*%diag(K) <= b1,
           quad_form(y,K) <= b2)
    if (sexed)
      constraints <- c(constraints, sum(parents$female*y)==0.5)
    problem <- Problem(objective,constraints)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      y.opt <- result$getValue(y)
      x.opt <- result$getValue(x)
      oc[,i] <- y.opt
      om[,i] <- x.opt
      F1 <- as.numeric(Kvec%*%x.opt/theta2 + theta1*(ploidy*t(y.opt)%*%diag(K)-1)/(ploidy-1))
      response$dF.realized[i] <- (F1-F.avg)/(1-F.avg)
      response$merit[i] <- result$value
    }
  }
  colnames(oc) <- round(response$dF.realized,3)
  oc <- oc[,!is.na(colnames(oc)),drop=F]
  colnames(om) <- round(response$dF.realized,3)
  om <- om[,!is.na(colnames(om)),drop=F]
  return(list(response=response, om=om, oc=oc))
}