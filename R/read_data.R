#' Read data for OCS/OMA
#' 
#' Constructs the kinship matrix and predicts merit for OCS/OMA
#' 
#' The first column of \code{geno.file} is marker name. The second column contains the additive effects for the breeding value parameterization, and (digenic) dominance effects should be in the third column with the header "dom". When \code{kinship.file=NULL}, the software assumes the following column "p.ref" contains allele frequencies for the reference population to control genomic inbreeding. Subsequent columns contain the marker data for the population, coded as allele dosage, from 0 to ploidy. 
#' 
#' There are several options for argument \code{matings}: (1) "none" = no matings; (2) "all" = all possible matings of the individuals in \code{geno.file} (excluding reciprocals); (3) a character vector of genotype ids to calculate all pairs of matings; (4) a data.frame of desired matings with 2 columns: female, male. 
#'
#' @param geno.file file with marker effects and genotypes
#' @param kinship.file NULL
#' @param ploidy even integer
#' @param sex optional, data frame with columns id and female (T/F)
#' @param matings see Details 
#' @param standardize T/F, standardize merit in parental candidates
#' @param n.core multi-core evaluation
#'
#' @return list containing
#' \describe{
#' \item{K}{genomic kinship matrix}
#' \item{parents}{data frame of individual merits}
#' \item{matings}{data frame of mating merits}
#'}
#' @export
#' @importFrom utils read.csv
#' @importFrom stats sd
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster


read_data <- function(geno.file, kinship.file=NULL, ploidy, sex=NULL, 
                      matings="none", standardize=FALSE, n.core=1) {
  
  data <- read.csv(file = geno.file,check.names=F,row.names=1)
  dominance <- (colnames(data)[2]=="dom")
  geno.start <- 2 + as.integer(dominance) + as.integer(is.null(kinship.file))
  
  effects <- as.matrix(data[,1:(as.integer(dominance)+1),drop=FALSE])
  geno <- as.matrix(data[,geno.start:ncol(data)])
  
  if (is.null(kinship.file)) {
    p <- as.numeric(data[,as.integer(dominance)+2])
  } else {
    p <- apply(geno,1,mean,na.rm=T)/ploidy
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
  coeff <- t(geno) - ploidy * Pmat 
  dimnames(coeff) <- list(id,markers)
  coeff[which(is.na(coeff))] <- 0
  
  if (is.null(kinship.file)) {
    scale <- ploidy*sum(p*(1-p))
    K <- tcrossprod(coeff)/scale/ploidy
    w <- 1e-5 
    K <- (1-w)*K + w*mean(diag(K))*diag(n)
  } else {
    K <- as.matrix(read.csv(kinship.file,row.names=1,check.names=F))
    colnames(K) <- rownames(K)
    stopifnot(id %in% colnames(K))
    K <- K[id,id]
  }
  
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
  mean.merit <- mean(parents$merit)
  sd.merit <- sd(parents$merit)
  if (standardize) 
    parents$merit <- (parents$merit - mean.merit)/sd.merit
  
  if (!is.null(sex)) {
    colnames(sex) <- c("id","female")
    parents <- merge(parents,sex)
  }
  
  MPH <- function(parents,ploidy,geno) {
    Xi <- geno[,parents[1]]
    Xj <- geno[,parents[2]]
    ploidy/4/(ploidy-1)*((Xi-Xj)^2 + 2/ploidy*Xi*Xj - (Xi+Xj))
  }
  
  if (matings=="none") {
    
    return(list(K=K, parents=parents[,setdiff(colnames(parents),c("add","dom"))]))
    
  } else {
    #OMA
    if (inherits(matings,"character")) {
      if (matings=="all")
        matings <- parents$id
      if ("sex" %in% colnames(parents)) {
        females <- intersect(parents$id[parents$female],matings)
        males <- intersect(parents$id[!parents$female],matings)
        stopifnot(length(females)>0 & length(males)>0)
        matings <- data.frame(expand.grid(female=females,male=males,stringsAsFactors = F))
      } else {
        id2 <- intersect(matings,parents$id)
        stopifnot(length(id2)>1)
        matings <- data.frame(expand.grid(female=id2,male=id2,stringsAsFactors = F))
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
    
    if (standardize)
      matings$merit <- (matings$merit-mean.merit)/sd.merit
    
    return(list(K=K, parents=parents[,setdiff(colnames(parents),c("add","dom"))], 
                matings=matings[,c("female","male","merit")]))
  } 
}
