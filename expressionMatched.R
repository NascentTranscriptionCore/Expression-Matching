expressionMatched <- function(geneset,ctrlgeneset,dat,numQuantGroups,CtrlSetRatio=1,setSeedNum=1){  
  #pull out the relevant columns from the dataframe given
  geneIDs <- as.character(dat[,1])
  expressionValues <- as.numeric(dat[,2])
  
  
  #pull out geneset quantile values
  de_set_quantiles <- quantile(expressionValues[geneIDs %in% as.character(geneset)],seq(0,1,1/numQuantGroups))
  
  #calculate the number of ctrl genes to pull per quantile group
  num_per_quantile <- CtrlSetRatio*round(length(geneset)/numQuantGroups,0)
  
  #create an empty character object to add expression matched genes too
  ctrlGenes <- character()
  
  #loop to pull out genes from each quantile group
  for (i in 1:(length(de_set_quantiles)-1)){
    #use set.seed to make the random sampling reproducible
    set.seed(setSeedNum)
    #Randomly select "num_per_quantile" genes from the ctrlgeneset for each quantile group of geneset expression
    x <- sample(geneIDs[geneIDs %in% as.character(ctrlgeneset) & 
                          expressionValues >= de_set_quantiles[i] & 
                          expressionValues < de_set_quantiles[(i+1)]],num_per_quantile)
    #add the genes to the growing vector of ctrlGenes                              
    ctrlGenes <- c(ctrlGenes,x)
    
  }
  #return the finished list of ctrlGenes
  return(ctrlGenes)
}