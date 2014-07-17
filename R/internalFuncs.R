##################################################
##################################################
### An internal function to check the inputs.
##################################################

GSReg.Check.input <- function(V,pathways,exprsdata,phenotypes,genenames,prunedpathways)
{
  if( !missing(V))
  {
    if(!is.numeric(V) || !is.matrix(V))
      stop("V must be a numeric matrix. Please transform vectors it to column matrices.")
    if(sum(!is.finite(V))>0)
      stop("V must contain finite numbers. NAs,NaN's, and Inf are unacceptable")
  }
  if(!missing(pathways))
  {
    if(!is.list(pathways) || sum(!sapply(pathways,is.character)))
      stop("pathways must be list of vectors of characters.")
    if(length(names(pathways))==0)
      stop("pathways require names. If you do not know the name of the pathway assign numbers as the names.")
  }
  if(!missing(exprsdata))
    if(!is.numeric(exprsdata) || ! is.matrix(exprsdata) || is.null(rownames(exprsdata)) )
      stop("exprsdata must be a numeric matrix with gene names as the rownames");
  if(!missing(phenotypes))
    if(!is.factor(phenotypes) || length(levels(phenotypes))!=2)
      stop("The phenotype must be factor consisting of only two levels.")
  if(!missing(exprsdata) && !missing(phenotypes))
    if(length(phenotypes)!=dim(exprsdata)[2])
      stop("The length of phenotypes and the column number of exprsdata must match.")
  if(!missing(genenames))
    if(!is.character(genenames) || ! is.vector(genenames))
      stop("genenames must be a vector of character names.")
  if(!missing(prunedpathways))
  {
    if(!is.list(prunedpathways) || sum(!sapply(prunedpathways,is.character)))
      stop("pathways must be list of vectors of characters.")
    if(length(names(prunedpathways)) == 0)
      stop("pathways require names. If you do not know the name of the pathway, please assign numbers as the names.")
    if(sum(sapply(prunedpathways,length)<2))
      stop("The pruned pathways must contain at least two genes. Please remove pathways with less than two genes.")
  }
  if(!missing(genenames) && !missing(prunedpathways))
  {
    if(length(setdiff(unlist(prunedpathways),genenames)))
      stop("The genes in the pathways must be a subset of the genes in the data expression. Use GSReg.Prune. ")
    sum(sapply(prunedpathways,length)<2)
  }  
}

##################################################
##################################################
### internal function for calculating the variance
### of only one pathway. Pathwayexp contains 
### expression data. samplesC1 and sampleC2 are the
### indices for class 1 and class 2.
##################################################
GSReg.Variance <- function(pathwayexp, samplesC1, samplesC2,distFunc =GSReg.kendall.tau.distance ){
  dist1 <- distFunc(pathwayexp[,samplesC1])
  n1 <- length(samplesC1)
  E1 <- sum(dist1)/(n1*(n1-1))
  Dis1P2 <- sum(dist1^2)
  E1P2 <- Dis1P2/(n1*(n1-1))
  E1DCross <- (sum(apply(dist1,1,sum)^2)-Dis1P2)/(n1*(n1-1)*(n1-2))
  VarD1 <- E1P2 - E1^2
  CovE1Cross <- E1DCross - E1^2
  VarEta1 <- 4*((n1*(n1-1))/2*VarD1 + n1*(n1-1)*(n1-2)*CovE1Cross)/((n1*(n1-2))^2)
  
  dist2 <- GSReg.kendall.tau.distance(pathwayexp[,samplesC2])
  n2 <- length(samplesC2)
  E2 <- sum(dist2)/(n2*(n2-1))
  Dis2P2 <- sum(dist2^2)
  E2P2 <- Dis2P2/(n2*(n2-1))
  E2DCross <- (sum(apply(dist2,1,sum)^2)-Dis2P2)/(n2*(n2-1)*(n2-2))
  VarD2 <- E2P2 - E2^2
  CovE2Cross <- E2DCross - E2^2
  VarEta2 <- 4*((n2*(n2-1))/2*VarD2 + n2*(n2-1)*(n2-2)*CovE2Cross)/((n2*(n2-2))^2)
  
  vartotal <- VarEta2 + VarEta1;
  zscore <- (E1-E2)/sqrt(abs(vartotal+0.0000001));
  return(list(E1=E1,E2=E2,zscore=zscore,
              VarEta1=VarEta1,VarEta2=VarEta2,
              sdtotal=sqrt(abs(vartotal)),pvalue=2*(1-pnorm(abs(zscore)))))
  
}
##################################################
##################################################
### Calculate the normalized DIRAC
### columns of V
### Outcome: \frac{\sum_{i,j} max(P(X_i<X_j),P(X_i>X_j))}{\sum_i 1}
##################################################
GSReg.DIRAC.mu<-function(V)
{
  #Checking if V is a numeric matrix 
  GSReg.Check.input(V)
  n <- dim(V)[1]
  m <- dim(V)[2]
  d <-.C("Nij",
        as.double(V),
        as.integer(dim(V)[1]),
        as.integer(dim(V)[2]),
        as.double(matrix(data= 0 ,ncol=n,nrow=n)))
  pij = d[[4]]/m   #pij = p(V_i<V_j)

  dim(pij)<-c(n,n)
  
  return(sum(pmin(pij,t(pij)))/(n*(n-1)))
}
####################################
####################################
##### Implements DIRAC analysis for 
##### one pathway.
####################################
####################################
GSReg.DIRAC.Pathways<-function(geneexpres,pathways,phenotypes){

  samplesC1 <- which(phenotypes==levels(phenotypes)[1])
  samplesC2 <- which(phenotypes==levels(phenotypes)[2])
  
  
  mu1 <- sapply(pathways,function(x) GSReg.DIRAC.mu(geneexpres[x,samplesC1]))
  mu2 <- sapply(pathways,function(x) GSReg.DIRAC.mu(geneexpres[x,samplesC2]))
  
  names(mu1) <- names(pathways)
  names(mu2) <- names(pathways)
  
  return(list(mu1=mu1,mu2=mu2,diffmu=mu1-mu2))
}
##################################################
##################################################
## GSReg.Prune function prunes the pathway list 
## based on the genenames 
####################################################
GSReg.Prune <- function(pathways,genenames, minGeneNum=5)
{
  prunedpathways <- lapply( pathways,function(x) intersect(x,genenames))
  emptypathways <- which(sapply(prunedpathways,length)<minGeneNum)
  if(length(emptypathways)>0)
  {
    message(paste("After pruning the pathways, there exist pathways with zero or one gene!\n Small pathways were deleted. Deleted pathways:","\n",paste(names(emptypathways),collapse='\n')))
    return(prunedpathways[-emptypathways])  
  }else{ return(prunedpathways)}
}

##################################################
##################################################
### Calculate the normalized kendall-tau-distance between
### columns of V
### Outcome: Z(i,j) = #(I(I(V_l^i<V_k^i)!=I(V_l^j<V_k^j)))/(n(n-1)/2)
##################################################

GSReg.kendall.tau.distance <- function(V){
  #Checking if V is a numeric matrix 
  GSReg.Check.input(V)
  
  n <- dim(V)[1]
  m <- dim(V)[2]
  
  d <- .C("kendalltaudist",
        as.double(V),
        as.integer(dim(V)[1]),
        as.integer(dim(V)[2]),
        as.double(matrix(data= 0 ,ncol=m,nrow=m)))
  dist <- d[[4]]
  dim(dist) <- c(m,m)
  return(2*dist/(n*(n-1)))
}