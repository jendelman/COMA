#' Optimal Mate Allocation
#'
#' Optimize the allocation for each mating
#' 
#' The data frame \code{parents} needs columns id, min, max, with optional column female (logical) for separate sexes. Additional columns can be used to specify constraints on linear combinations of the contributions. The values are the coefficients of the contribution variables, and the name of each column specifies the right-hand side of the constraint. Each name must begin with "lt","gt", or "eq", followed by a non-negative numeric value. For example, "lt0.5" means less than or equal to 0.5.
#' 
#' The data.frame \code{matings} has five columns: female, male, merit, min, max. 
#' 
#' The average inbreeding coefficient of the current generation is based on all individuals in \code{K}, which may exceed the list of individuals in \code{parents}.
#'
#' After optimization, allocations less than \code{min.a} are set equal to zero, and the remainder are renormalized. The realized inbreeding rate can exceed the specified limit after applying this threshold. It is also possible that no feasible solution exists for the specified \code{dF}. In either case, argument \code{dF.adapt} can be used to find other solutions. It is a list with two named components: step, max. The software increases the dF limit by dF.adapt$step up to the smaller of dF.adapt$max or the realized value under the original dF, in an attempt to find a solution with less inbreeding. 
#'
#' @param dF inbreeding rate
#' @param parents parents data frame (see Details)
#' @param matings matings data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param min.a minimum allocation
#' @param dF.adapt see Details
#' @param solver solver for CVXR (default is "ECOS")
#' 
#' @return list containing
#' \describe{
#' \item{response}{named vector with realized dF, merit, Shannon diversity for parents}
#' \item{oc}{data frame of optimal contributions for each individual}
#' \item{om}{data frame of optimal allocations for each mating}
#' }
#' @import CVXR
#' @importFrom stats model.matrix
#' @export
#' 
oma <- function(dF, parents, matings, ploidy, K, min.a=0, 
                dF.adapt=NULL, solver="ECOS") {
  
  stopifnot(c("id","merit","min","max") %in% colnames(parents))
  parents$id <- as.character(parents$id)
  stopifnot(parents$max <= 1)
  stopifnot(parents$min >= 0)
  stopifnot(c("female","male","merit","min","max") %in% colnames(matings))
  stopifnot(matings$max <= 1)
  stopifnot(matings$min >= 0)
  
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  F.avg <- mean(as.numeric(Fi))
  
  h.t <- matrix(matings$merit,nrow=1)
  parent.id <- sort(union(matings$female,matings$male))
  ix <- match(parent.id,parents$id,nomatch=0)
  if (any(ix==0))
    stop("parents data.frame missing individuals in matings")
  parents <- parents[ix,]
  n <- nrow(parents)
  p <- nrow(matings)
  
  matings$female <- factor(matings$female,levels=parent.id)
  matings$male <- factor(matings$male,levels=parent.id)

  M1 <- model.matrix(~female-1,matings)
  M2 <- model.matrix(~male-1,matings)
  
  matings$female <- as.character(matings$female)
  matings$male <- as.character(matings$male)
  
  M <- t((M1+M2)/2)
  rownames(M) <- parent.id
  
  oc <- data.frame(id=parents$id, value=numeric(n))
  om <- data.frame(matings[,1:2], value=numeric(p))
  
  K <- K[parents$id,parents$id]
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  Kvec <- matrix(apply(matings[,c("female","male")],1,
                       function(z,K){K[z[1],z[2]]},K=K),nrow=1)
  
  y <- Variable(n)
  x <- Variable(p)
  objective <- Maximize(h.t%*%x)
  
  constraints <- NULL
  if ("female" %in% colnames(parents)) {
    sexed <- TRUE
    stopifnot(!duplicated(parents$id))
    stopifnot(class(parents$female)=="logical")
    stopifnot(sum(parents$female) > 0 & sum(!parents$female) > 0)
    parents$female <- as.integer(parents$female)
    if (ncol(parents) > 5) 
      constraints <- parents[,6:ncol(parents),drop=FALSE]
  } else {
    sexed <- FALSE
    if (ncol(parents) > 4) 
      constraints <- parents[,5:ncol(parents),drop=FALSE]
  }
  con.list <- list(y <= parents$max, y >= parents$min, y==M%*%x, 
                   x <= matings$max, x >= matings$min, sum(x)==1)
  if (sexed)
    con.list <- c(con.list, sum(parents$female*y)==0.5)
  if (!is.null(constraints)) {
    tmp <- colnames(constraints)
    nc <- length(tmp)
    signs <- substr(tmp,1,2)
    stopifnot(signs %in% c("lt","gt","eq"))
    lt <- which(signs=="lt")
    gt <- which(signs=="gt")
    eq <- which(signs=="eq")
    values <- apply(array(tmp),1,
                    function(z){as.numeric(substr(z,3,nchar(z)))})
    if (length(gt) > 0) {
      values[gt] <- -values[gt]
      constraints[,gt] <- -constraints[,gt]
      signs[gt] <- "lt"
      lt <- which(signs=="lt")
    }
    if (length(lt) > 0) {
      con.list <- c(con.list, 
                    t(as.matrix(constraints[,lt,drop=FALSE])) %*% y <= values[lt])
    }
    if (length(eq) > 0) {
      con.list <- c(con.list, 
                    t(as.matrix(constraints[,eq,drop=FALSE])) %*% y == values[eq])
    }
  }
  
  dF1 <- dF
  dFr <- merit <- as.numeric(NA)
  dF.best <- Inf
  flag <- TRUE
  theta1 <- ploidy*(3*ploidy/4-1)
  while (flag) {
    theta2 <- (ploidy/2-1)*(dF1*(ploidy-1)-ploidy/2)
    theta3 <- ploidy/2*(ploidy-1)*(1-dF1)
    RHS1 <- (ploidy-1)*(dF1+(1-dF1)*F.avg)
    RHS2 <- (ploidy-1)^2*dF1
    con2 <- c(con.list, 
      list((ploidy/2)*Kvec%*%x + (ploidy/2-1)*Fi%*%y <= RHS1,
           theta1*quad_form(y,K) + theta2*Fi%*%y - theta3*Kvec%*%x <= RHS2))
    
    problem <- Problem(objective,con2)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      x.opt <- as.numeric(result$getValue(x))
      ix <- which(x.opt < min.a)
      if (length(ix) > 0)
        x.opt[ix] <- 0
      if (sum(x.opt) > 0) {
        x.opt <- x.opt/sum(x.opt)
        y.opt <- as.numeric(M%*%x.opt)
        
        F.t1 <- as.numeric((ploidy/2)*Kvec%*%x.opt + (ploidy/2-1)*Fi%*%y.opt)/(ploidy-1)
        dFr <- (F.t1-F.avg)/(1-F.avg)
        
        F.t2 <- as.numeric(theta1*crossprod(y.opt,K%*%y.opt) + (ploidy/2-1)^2*Fi%*%y.opt)/(ploidy-1)^2
        dF2r <- (F.t2-F.t1)/(1-F.t1)
        if (dFr < dF.best) {
          oc$value <- y.opt
          om$value <- x.opt
          merit <- result$value
          dF.best <- dFr
          dF2.best <- dF2r
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
    return(list(response=c(dF=as.numeric(NA), merit=as.numeric(NA), diversity=as.numeric(NA)), 
                oc=oc[integer(0),], om=om[integer(0),]))
  } else {
    oc <- oc[which(oc$value > 0),,drop=FALSE]
    om <- om[which(om$value > 0),,drop=FALSE]
    diversity <- -sum(oc$value*log(oc$value))
    return(list(response=c(dF1=dF.best, 
                           dF2=dF2.best,
                           merit=merit,
                           diversity=diversity),oc=oc,om=om)) 
  }
}
