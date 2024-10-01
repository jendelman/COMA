#' Optimal Contribution Selection
#'
#' Optimize the contribution of each individual to the next generation
#' 
#' The first four columns of \code{parents} should be named as follows: id, merit, min, max. Min and max are real numbers between 0 and 1 specifying the minimum and maximum contribution for each parent. An optional fifth column named "female" is a logical TRUE/FALSE variable for species with separate sexes. 
#' 
#' Additional columns can be used to specify constraints on linear combinations of the contributions. The values are the coefficients of the contribution variables, and the name of each column specifies the right-hand side of the constraint. Each name must begin with "lt","gt", or "eq", followed by a non-negative numeric value. For example, "lt0.5" means less than or equal to 0.5.
#' 
#' The average inbreeding coefficient of the current generation is based on all individuals in \code{K}, which may exceed the list of individuals in \code{parents}.
#' 
#' It is possible that no feasible solution exists for the specified \code{dF}. Argument \code{dF.adapt} can be used to automatically increase dF by dF.adapt$step up to dF.adapt$max.
#' 
#' @param dF inbreeding rate 
#' @param parents input data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param tol tolerance, values below this set to 0
#' @param dF.adapt see Details
#' @param solver solver for CVXR (default is "ECOS")
#' 
#' @return list containing
#' \describe{
#' \item{response}{data.frame with realized dF, merit, n.parent}
#' \item{oc}{data frame of optimal contributions}
#' }
#' @import CVXR
#' @export
#' 
ocs <- function(dF, parents, ploidy, K, tol=1e-6, 
                dF.adapt=NULL, solver="ECOS") {
  
  stopifnot(c("id","merit","min","max") == colnames(parents)[1:4])
  stopifnot(parents$max <= 1)
  stopifnot(parents$min >= 0)
  ploidy <- as.integer(ploidy)
  
  parents$id <- as.character(parents$id)
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
  
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  Ft0 <- mean(as.numeric(Fi))
  n <- nrow(parents)
  
  oc <- data.frame(id=parents$id, value=numeric(n),check.names=F)
  
  h.t <- matrix(parents$merit,nrow=1)
  K <- K[parents$id,parents$id]
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  
  y <- Variable(n)
  objective <- Maximize(h.t%*%y)
  con.list <- list(y <= parents$max, y >= parents$min, sum(y)==1)
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
  
  done <- FALSE
  while (!done) { 
    Ft1 <- dF+(1-dF)*Ft0
    
    con2 <- c(con.list,
              (ploidy/2)*quad_form(y,K) + (ploidy/2-1)*Fi%*%y <= (ploidy-1)*Ft1)

    problem <- Problem(objective,con2)
    result <- suppressWarnings(solve(problem,solver=solver))
    if (result$status %in% c("optimal","optimal_inaccurate")) {
      done <- TRUE
    } else {
      if (is.null(dF.adapt)) {
        done <- TRUE
      } else {
        dF <- dF + dF.adapt$step
        if (dF > dF.adapt$max)
          done <- TRUE
      }
    }
  }
  
  if (result$status %in% c("optimal","optimal_inaccurate")) {
    if (result$status=="optimal_inaccurate")
      warning("Optimal inaccurate solution")
    
    y.opt <- as.numeric(result$getValue(y))
    y.opt[y.opt < tol] <- 0
    y.opt <- y.opt/sum(y.opt)
    
    if (sexed) {
      y.f <- 2*parents$female*y.opt
      y.m <- 2*(1-parents$female)*y.opt
    } else {
      y.f <- y.m <- y.opt
    }
        
    Ft1 <- as.numeric((ploidy/2)*crossprod(y.f,K%*%y.m) + 
                         (ploidy/2-1)*Fi%*%y.opt)/(ploidy-1)
    dF1 <- (Ft1-Ft0)/(1-Ft0)

    oc$value <- y.opt
    oc <- oc[which(oc$value >= tol),,drop=FALSE]
    
    return(list(response=data.frame(dF=round(dF1,4), 
                           merit=result$value,
                           n.parent=nrow(oc)),oc=oc)) 
  } else {
    return(list(response=data.frame(dF=as.numeric(NA),
                           merit=as.numeric(NA),
                           n.parent=0L),
                oc=oc[integer(0),]))
  }
}