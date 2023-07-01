#' Read data for OCS/OMA
#' 
#' Constructs the G matrix and predicts merit for OCS/OMA
#' 
#' The first column of \code{geno.file} is marker name. The second column contains the additive effects for the breeding value parameterization, and when \code{dominance=TRUE}, the next column contains (digenic) dominance effects. The following column "p.ref" should contain allele frequencies for the reference population, which is needed to control genomic inbreeding. Subsequent columns contain the marker data for the population, coded as allele dosage, from 0 to ploidy. 
#' 
#' The argument \code{matings} can be either (1) a character vector of genotype ids, in which case all possible matings of these individuals are generated (reciprocal matings excluded), or (2) a data.frame with four columns: female, male, min, max. The columns "min" and "max" specify the limits (0-1) for the contribution of each mating. To skip mate allocation entirely, use \code{matings=NULL}.
#'
#' @param geno.file file with marker effects and genotypes
#' @param ploidy even integer
#' @param dominance TRUE/FALSE
#' @param sex optional, data frame with columns id (character) and female (TRUE/FALSE)
#' @param matings see Details 
#' @param n.core multi-core evaluation
#'
#' @return list containing
#' \describe{
#' \item{K}{genomic kinship matrix}
#' \item{parents}{data frame with individual merits and limits}
#' \item{matings}{data frame with mating merits and limits}
#'}
#' @export
#' @importFrom utils read.csv
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster


read_data <- function(geno.file, ploidy, dominance, sex=NULL, matings, n.core=1) {
  
  data <- read.csv(file = geno.file,check.names=F,row.names=1)
  if (dominance) {
    effects <- as.matrix(data[,1:2])
    geno <- as.matrix(data[,-(1:3)])
    p <- as.numeric(data[,3])
  } else {
    effects <- as.matrix(data[,1,drop=FALSE])
    geno <- as.matrix(data[,-(1:2)])
    p <- as.numeric(data[,2])
  }
  rownames(geno) <- rownames(data)
  n <- ncol(geno)
  id <- colnames(geno)
  polymorphic <- which(p >= 1/(ploidy*n))
  m <- length(polymorphic)
  stopifnot(m > 0)
  p <- p[polymorphic]
  geno <- geno[polymorphic,]
  effects <- effects[polymorphic,,drop=FALSE]
  markers <- rownames(geno)
  
  Pmat <- kronecker(matrix(p,nrow=1,ncol=m),matrix(1,nrow=n,ncol=1))
  #p.mat <- matrix(1,nrow=n,ncol=1) %*% matrix(p,nrow=1)
  coeff <- t(geno) - ploidy * Pmat 
  #coeff <- scale(t(geno),scale=F)
  dimnames(coeff) <- list(id,markers)
  coeff[which(is.na(coeff))] <- 0
  scale <- ploidy*sum(p*(1-p))
  K <- tcrossprod(coeff)/scale/ploidy
  w <- 1e-5 
  K <- (1-w)*K + w*mean(diag(K))*diag(n)
  
  #predict merit
  #OCS
  parents <- data.frame(id=id, add=as.numeric(coeff %*% effects[,1,drop=FALSE]))
  
  if (dominance) {
    coeff.D <- -2*choose(ploidy,2)*Pmat^2 + 2*(ploidy-1)*Pmat*t(geno) - t(geno)*(t(geno)-1)
    coeff.D[is.na(coeff.D)] <- 0
    gamma <- (ploidy/2 - 1)/(ploidy - 1)
    parents$dom <- as.numeric(coeff.D %*% effects[,2,drop=FALSE])
    parents$merit <- parents$add + gamma*parents$dom
  } else {
    parents$merit <- parents$add
  }
  
  if (!is.null(sex)) {
    colnames(sex) <- c("id","female")
    parents <- merge(parents,sex)
  }
  
  MPH <- function(parents,ploidy,geno) {
    Xi <- geno[,parents[1]]
    Xj <- geno[,parents[2]]
    ploidy/4/(ploidy-1)*((Xi-Xj)^2 + 2/ploidy*Xi*Xj - (Xi+Xj))
  }
  
  if (is.null(matings)) {
    
    return(list(K=K, parents=parents[,setdiff(colnames(parents),c("add","dom"))]))
    
  } else {
    #OMA
    if (class(matings)=="character") {
      if ("sex" %in% colnames(parents)) {
        females <- intersect(parents$id[parents$female],matings)
        males <- intersect(parents$id[!parents$female],matings)
        stopifnot(length(females)>0 & length(males)>0)
        matings <- data.frame(expand.grid(female=females,male=males,stringsAsFactors = F),
                              min=0,max=1)
      } else {
        id2 <- intersect(matings,id)
        stopifnot(length(id2)>1)
        matings <- data.frame(expand.grid(female=id2,male=id2,min=0,max=1,stringsAsFactors = F))
        matings <- matings[matings$female >= matings$male,] 
      }
    }
    stopifnot(nrow(matings)>0)
    matings$female <- as.character(matings$female)
    matings$male <- as.character(matings$male)
    id2 <- union(matings$female,matings$male)
    stopifnot(id2 %in% parents$id)
    parents <- parents[parents$id %in% id2,]
    matings$merit <- (parents$add[match(matings$female,parents$id)] +
                        parents$add[match(matings$male,parents$id)])/2
    
    if (dominance) {
      matings$merit <- matings$merit + 
                        (parents$dom[match(matings$female,parents$id)] +
                         parents$dom[match(matings$male,parents$id)])/2
      
      mate.list <- split(as.matrix(matings[,1:2]),f=1:nrow(matings))
      if (n.core > 1) {
        cl <- makeCluster(n.core)
        clusterExport(cl=cl,varlist=NULL)
        ans <- parSapply(cl=cl,X=mate.list,MPH,ploidy=ploidy,geno=geno)
        stopCluster(cl)
      } else {
        ans <- sapply(mate.list,MPH,ploidy=ploidy,geno=geno)
        #ans[is.na(ans)] <- 0
      }
      matings$merit <- matings$merit + as.numeric(crossprod(ans,effects[,2]))
    }
    return(list(K=K, parents=parents[,setdiff(colnames(parents),c("add","dom"))], 
                matings=matings[,c("female","male","merit","min","max")]))
  } 
}
