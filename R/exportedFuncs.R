############################################
##########        GSReg Package   ##########
##########        Bahman Afsari   ##########
##########        Elana J. Fertig ##########

###########################################
#### GSVReg Analysis Function
###########################################

GSReg.GeneSets.VReg <- function(geneexpres,pathways,phenotypes,minGeneNum = 5 )
{
  #Checking the input data
  GSReg.Check.input(prunedpathways=pathways,exprsdata=geneexpres,phenotypes=phenotypes)
  #pruning the pathways
  pathways <- GSReg.Prune(pathways,rownames(geneexpres), minGeneNum)
  values <- vector(mode="list",length(pathways))
  names(values) <- names(pathways)
  samplesC1 <- which(phenotypes==levels(phenotypes)[1])
  samplesC2 <- which(phenotypes==levels(phenotypes)[2])
  for(i in seq_along(pathways))
  {
    values[[i]] <- GSReg.Variance(geneexpres[pathways[[i]],], samplesC1, samplesC2,GSReg.kendall.tau.distance)
  }
  return(values)
}

###########################################
###########################################
#### DIRAC Analysis Function
###########################################
###########################################
GSReg.GeneSets.DIRAC <- function(geneexpres,pathways,phenotypes,Nperm=1000,minGeneNum = 5)
{
  GSReg.Check.input(prunedpathways=pathways,exprsdata=geneexpres,phenotypes=phenotypes)
  #prune genes in the pathway
  pathways <- GSReg.Prune(pathways,rownames(geneexpres), minGeneNum)
  #calculate DIRAC variability measure
  mus <- GSReg.DIRAC.Pathways(geneexpres,pathways,phenotypes)
  #Calculating a p-value
  if(Nperm > 0)
  {
    musperm <- vector(mode="list",Nperm)
    for( i in seq_len(Nperm))
    {
      musperm[[i]] <- GSReg.DIRAC.Pathways(geneexpres=geneexpres, pathways=pathways, phenotypes=sample(phenotypes))
    }
    pvaluesperm <- vector(mode="numeric",length=length(pathways))
    names(pvaluesperm) <- names(pathways) 
    for( i in seq_along(pathways))
    {
      z <- sapply(musperm,function(x) x$diffmu[i])
      pvaluesperm[i] <- mean(abs(mus$mu1[i]-mus$mu2[i])<=abs(z)) 
    }
    return(list(mu1=mus$mu1,mu2=mus$mu2,pvalues=pvaluesperm))
  }else{
    return(list(mu1=mus$mu1,mu2=mus$mu2))
  }
}