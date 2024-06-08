#' Optimal Mate Allocation
#'
#' Optimize the allocation for each mating
#' 
#' The first three columns of \code{parents} should be named "id", "min", "max", with an optional fourth column "female" to indicate sex in dioecious species. Additional columns can be used to specify constraints on linear combinations of the contributions. The values are the coefficients of the contribution variables, and the name of each column specifies the right-hand side of the constraint. Each name must begin with "lt","gt", or "eq", followed by a non-negative numeric value. For example, "lt0.5" means less than or equal to 0.5.
#' 
#' The data.frame \code{matings} has five columns: "female, male, merit, min, max" for dioecious species, or else "parent1, parent2, merit, min, max".
#' 
#' The average inbreeding coefficient of the current generation is based on all individuals in \code{K}, which may exceed the list of individuals in \code{parents}.
#'
#' It is possible that no feasible solution exists for the specified \code{dF}. Argument \code{dF.adapt} can be used to automatically progressively higher values.The software increases the dF limit by dF.adapt$step up to the smaller of dF.adapt$max or the realized value under the original dF, in an attempt to find a solution with less inbreeding. 
#'
#' @param dF inbreeding rate
#' @param parents parents data frame (see Details)
#' @param matings matings data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param tol tolerance, values below this set to 0
#' @param dF.adapt see Details
#' @param solver solver for CVXR (default is "ECOS")
#' 
#' @return list containing
#' \describe{
#' \item{response}{data.frame with realized dF, merit, n.parent, n.mate}
#' \item{oc}{data frame of optimal contributions for each individual}
#' \item{om}{data frame of optimal allocations for each mating}
#' }
#' @import CVXR
#' @importFrom stats model.matrix
#' @export
#' 
oma <- function(dF, parents, matings, ploidy, K,
                tol=1e-6, dF.adapt=NULL, solver="ECOS") {
  
  stopifnot(c("id","min","max") == colnames(parents)[1:3])
  parents$id <- as.character(parents$id)
  stopifnot(parents$max <= 1)
  stopifnot(parents$min >= 0)
  stopifnot(c("merit","min","max") == colnames(matings)[3:5])
  stopifnot(matings$max <= 1)
  stopifnot(matings$min >= 0)
  ploidy <- as.integer(ploidy)
                       
  colnames(matings) <- c("female","male","merit","min","max")
  
  Fi <- matrix((ploidy*diag(K)-1)/(ploidy-1),nrow=1)
  Ft0 <- mean(as.numeric(Fi))
  
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

  Lf <- model.matrix(~female-1,matings)
  Lm <- model.matrix(~male-1,matings)
  M <- t(Lf+Lm)/2
  rownames(M) <- parent.id
  
  matings$female <- as.character(matings$female)
  matings$male <- as.character(matings$male)
  
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
    if (ncol(parents) > 4) 
      constraints <- parents[,5:ncol(parents),drop=FALSE]
  } else {
    sexed <- FALSE
    if (ncol(parents) > 3) 
      constraints <- parents[,4:ncol(parents),drop=FALSE]
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
  
  done <- FALSE
  while (!done) {
    Ft1 <- dF+(1-dF)*Ft0
    Ft2 <- dF*(2-dF)+(1-dF)^2*Ft0
    
    con2 <- c(con.list, 
      list((ploidy/2)*Kvec%*%x + (ploidy/2-1)*Fi%*%y <= (ploidy-1)*Ft1,
           (ploidy/2)*(ploidy-1)*quad_form(y,K) + 
             (ploidy/2-1)*((ploidy/2-1)*Fi%*%y + ploidy/2*Kvec%*%x) <= (ploidy-1)^2*Ft2))
    
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
      
    x.opt <- as.numeric(result$getValue(x))
    x.opt[x.opt < tol] <- 0
    x.opt <- x.opt/sum(x.opt)
  
    y.opt <- as.numeric(M%*%x.opt)
    Ft1 <- as.numeric((ploidy/2)*Kvec%*%x.opt + (ploidy/2-1)*Fi%*%y.opt)/(ploidy-1)
    dF1 <- (Ft1-Ft0)/(1-Ft0)
    Ft2 <- ((ploidy/2)*(ploidy-1)*crossprod(y.opt,K%*%y.opt) + 
              (ploidy/2-1)*((ploidy/2-1)*Fi%*%y.opt + 
                              ploidy/2*Kvec%*%x.opt))/(ploidy-1)^2
    dF2 <- sqrt(1+(as.numeric(Ft2)-Ft0)/(1-Ft0)) - 1
    
    oc$value <- y.opt
    om$value <- x.opt
  
    om <- om[om$value >= tol,,drop=FALSE]
    oc <- oc[oc$id %in% union(om$female,om$male),,drop=FALSE]
    
    if (!sexed)
      colnames(om) <- replace(colnames(om),1:2,c("parent1","parent2"))
    
    #diversity <- -sum(oc$value*log(oc$value))
    return(list(response=data.frame(dF1=round(dF1,4),
                           dF2=round(dF2,4), 
                           merit=result$value,
                           n.parent=nrow(oc),
                           n.mate=nrow(om)),oc=oc,om=om)) 
  } else {
    return(list(response=data.frame(dF1=as.numeric(NA),dF2=as.numeric(NA),
                           merit=as.numeric(NA),
                           n.parent=0L, n.mate=0L),
                oc=oc[integer(0),],om=om[integer(0),]))
  }
}
