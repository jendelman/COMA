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
#' After optimization, contributions less than \code{min.c} are set equal to zero, and the remainder are renormalized. The realized inbreeding rate can exceed the specified limit after applying this threshold. It is also possible that no feasible solution exists for the specified \code{dF}. In either case, argument \code{dF.adapt} can be used to find other solutions. It is a list with two named components: step, max. The software increases the dF limit by dF.adapt$step up to the smaller of dF.adapt$max or the realized value under the original dF, in an attempt to find a solution with less inbreeding. 
#' 
#' @param dF inbreeding rate 
#' @param parents input data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param min.c minimum contribution
#' @param dF.adapt see Details
#' @param solver solver for CVXR (default is "ECOS")
#' 
#' @return list containing
#' \describe{
#' \item{response}{named vector with realized dF, merit, Shannon diversity for parents}
#' \item{oc}{data frame of optimal contributions}
#' }
#' @import CVXR
#' @export
#' 
ocs <- function(dF, parents, ploidy, K, min.c=0, 
                dF.adapt=NULL, solver="ECOS") {
  
  stopifnot(c("id","merit","min","max") == colnames(parents)[1:4])
  stopifnot(parents$max <= 1)
  stopifnot(parents$min >= 0)
  
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
  F.avg <- mean(as.numeric(Fi))
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
  
  dF1 <- dF
  dFr <- merit <- as.numeric(NA)
  dF.best <- Inf
  flag <- TRUE
  
  while (flag) { 
    RHS <- (ploidy-1)*(dF1+(1-dF1)*F.avg)
    con2 <- c(con.list,
             (ploidy/2)*quad_form(y,K) + (ploidy/2-1)*Fi%*%y <= RHS)

    problem <- Problem(objective,con2)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      y.opt <- as.numeric(result$getValue(y))
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
          oc$value <- y.opt
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
    return(list(response=c(dF=as.numeric(NA), 
                           merit=as.numeric(NA),
                           diversity=as.numeric(NA)), 
                oc=oc[integer(0),]))
  } else {
    oc <- oc[which(oc$value > 0),,drop=FALSE]
    diversity <- -sum(oc$value*log(oc$value))
    return(list(response=c(dF=dF.best, 
                           merit=merit,
                           diversity=diversity),oc=oc)) 
  }
}