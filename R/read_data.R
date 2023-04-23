#' Read data for OCS/OMA
#' 
#' Constructs the G matrix and predicts merit for OCS/OMA
#' 
#' The first column of \code{geno.file} is marker name. The second column contains the additive effects for the breeding value parameterization. When \code{dominance=TRUE}, the third column contains (digenic) dominance effects. Subsequent columns contain the marker data for the population, coded as allele dosage, from 0 to ploidy. 
#' 
#' When \code{matings} is present, merit is predicted for all matings in the data frame. The columns "lower" and "upper" specify the limits (0-1) for the contribution of each mating.
#'
#' @param geno.file file with marker effects and genotypes
#' @param ploidy even integer
#' @param dominance TRUE/FALSE
#' @param sex optional, data frame with columns id and sex ("M" or "F")
#' @param matings data frame with four columns: mother, father, lower, upper
#' @param n.core multi-core evaluation
#'
#' @return list containing
#' \describe{
#' \item{K}{genomic kinship matrix}
#' \item{ocs.data}{data frame with individual merits and limits}
#' \item{oma.data}{data frame with mating merits and limits}
#'}
#' @export
#' @importFrom utils read.csv
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster


read_data <- function(geno.file, ploidy, dominance, sex=NULL, matings=NULL, n.core=1) {
  
  data <- read.csv(file = geno.file,check.names=F,row.names=1)
  if (dominance) {
    effects <- as.matrix(data[,1:2])
    geno <- as.matrix(data[,-(1:2)])
  } else {
    effects <- matrix(data[,1],ncol=1)
    geno <- as.matrix(data[,-1])
  }
  rownames(geno) <- rownames(data)
  n <- ncol(geno)
  id <- colnames(geno)
  p <- apply(geno,1,mean,na.rm=T)/ploidy
  polymorphic <- which(p >= 1/(ploidy*n))
  m <- length(polymorphic)
  stopifnot(m > 0)
  p <- p[polymorphic]
  geno <- geno[polymorphic,]
  effects <- effects[polymorphic,]
  markers <- rownames(geno)
  
  coeff <- scale(t(geno),scale=F)
  dimnames(coeff) <- list(id,markers)
  coeff[which(is.na(coeff))] <- 0
  scale <- ploidy*sum(p*(1-p))
  K <- tcrossprod(coeff)/scale/ploidy
  w <- 1e-5 
  K <- (1-w)*K + w*mean(diag(K))*diag(n)
  
  #predict merit
  #OCS
  ocs.data <- data.frame(id=id, add=as.numeric(coeff %*% effects[,1]))
  
  if (dominance) {
    Pmat <- kronecker(matrix(p,nrow=1,ncol=m),matrix(1,nrow=n,ncol=1))
    coeff.D <- -2*choose(ploidy,2)*Pmat^2 + 2*(ploidy-1)*Pmat*t(geno) - t(geno)*(t(geno)-1)
    coeff.D[is.na(coeff.D)] <- 0
    gamma <- (ploidy/2 - 1)/(ploidy - 1)
    ocs.data$merit <- ocs.data$add + gamma*as.numeric(coeff.D %*% effects[,2])
  } else {
    ocs.data$merit <- ocs.data$add
  }
  
  if (!is.null(sex)) {
    colnames(sex) <- c("id","sex")
    ocs.data <- merge(ocs.data,sex)
  }
  
  EQz <- function(parents,ploidy,geno,p) {
    p1 <- geno[,parents[1]]/ploidy
    q1 <- 1-p1
    p2 <- geno[,parents[2]]/ploidy
    q2 <- 1-p2
    -ploidy^2/4*((p1+p2)^2+(p1*q1+p2*q2)/(ploidy-1)) + 
      ploidy*(p*(ploidy-1)+1/2)*(p1+p2) - ploidy*(ploidy-1)*p^2
  }
  
  if (is.null(matings)) {
    return(list(K=K, ocs.data=ocs.data[,-2]))
  } else {
    #OMA
    matings$mother <- as.character(matings$mother)
    matings$father <- as.character(matings$father)
    matings$merit <- (ocs.data$add[match(matings$mother,ocs.data$id)] + 
                        ocs.data$add[match(matings$father,ocs.data$id)])/2
    if (dominance) {
      mate.list <- split(as.matrix(matings[,1:2]),f=1:nrow(matings))
      if (n.core > 1) {
        cl <- makeCluster(n.core)
        clusterExport(cl=cl,varlist=NULL)
        ans <- parSapply(cl=cl,X=mate.list,EQz,ploidy=ploidy,geno=geno,p=p)
        stopCluster(cl)
      } else {
        ans <- sapply(mate.list,EQz,ploidy=ploidy,geno=geno,p=p)
        ans[is.na(ans)] <- 0
      }
      matings$merit <- matings$merit + as.numeric(effects[,2] %*% ans)
    }
    return(list(K=K, ocs.data=ocs.data[,-2], oma.data=matings))
  }
}
