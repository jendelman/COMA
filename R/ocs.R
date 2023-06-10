#' Optimal Contribution Selection
#'
#' Optimize the contribution of each individual to the next generation
#' 
#' The data.frame \code{parents} can have up to five columns: id, merit, min, max, female. Min and max are real numbers between 0 and 1 specifying the minimum and maximum contribution for each parent. The "female" column is logical TRUE/FALSE. For hermaphroditic species, the last column is omitted.
#' 
#' The average inbreeding coefficient of the current generation is based on all individuals in \code{K}.
#' 
#' @param parents input data frame (see Details)
#' @param ploidy ploidy
#' @param K kinship matrix 
#' @param dF inbreeding rate (can be vector of values)
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
ocs <- function(parents, ploidy, K, dF, solver="ECOS") {
  
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
  theta1 <- (ploidy/2-1)/(ploidy-1)
  theta2 <- (ploidy-1)/(ploidy/2)
  F.avg <- (ploidy*mean(diag(K))-1)/(ploidy-1)
  n <- nrow(parents)
  oc <- matrix(as.numeric(NA),nrow=n,ncol=m,dimnames=list(parents$id,dF))
  response <- data.frame(dF.target=dF, dF.realized=numeric(m)*NA, merit=numeric(m)*NA)
  h.t <- matrix(parents$merit,nrow=1)
  K <- K[parents$id,parents$id]
  
  for (i in 1:m) {
    b2 <- theta2*(2-dF[i]-theta1)*dF[i] + theta2*(1-dF[i])*(1-dF[i]-theta1)*F.avg
    y <- Variable(n)
    objective <- Maximize(h.t%*%y)
    constraints <- 
      list(y <= parents$max, y >= parents$min, sum(y)==1, quad_form(y,K) <= b2)
    if (sexed)
      constraints <- c(constraints, sum(parents$female*y)==0.5)
    problem <- Problem(objective,constraints)
    result <- solve(problem,solver="ECOS")
    if (result$status=="optimal") {
      y.opt <- result$getValue(y)
      oc[,i] <- y.opt
      if (sexed) {
        y.f <- parents$female*y.opt
        y.m <- (1-parents$female)*y.opt
      } else {
        y.f <- y.m <- y.opt
      }
      F1 <- as.numeric(t(y.f)%*%crossprod(K,y.m))
      response$dF.realized[i] <- (F1-F.avg)/(1-F.avg)
      response$merit[i] <- result$value
    }
  }
  colnames(oc) <- round(response$dF.realized,3)
  oc <- oc[,!is.na(colnames(oc))]
  
  #Fvec <- as.numeric(colnames(oc))
  # ix <- which(!is.na(Fvec) & !duplicated(Fvec))
  # if (length(ix)==0) {
  #   stop("No optimal solutions")
  # }
  # oc <- oc[,ix,drop=FALSE]
  # response <- response[ix,]
  # 
  return(list(response=response, oc=oc))
}