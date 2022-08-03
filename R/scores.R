#dt[,bplapply(.SD, fun), by=c(x,y,z), BPPARAM), .SDcols=ids]

.getTemplates <- function(x, pos, m){
  
  cellTable <- data.table(rate=x, pos=pos)
  
  # Mark adjacent CpGs
  isCovered <- !is.na(cellTable[,rate])
  temps <- frollsum(isCovered, m+1, align="right")
  cellTable$isAdj <- fifelse(temps==m+1 & !is.na(temps), T, F)
  
  # Construct (adjacent) templates
  cellTable[,temp_start:=fifelse(shift(isAdj,n=m, type="lead"), pos, NaN, na=NaN)]
  cellTable[,temp_end:=fifelse(!is.na(temp_start), shift(pos, n=m,type="lead"), NaN)] # temp_start+m
  
  tempCols <- as.character(0:m)
  cellTable[,(tempCols):=lapply(0:m, function(i) shift(rate, n=i, type="lead"))]
  
  templates <- cellTable[!is.nan(temp_start),]
  templates[,tempM:=paste(.SD, sep="", collapse=","), by=temp_start, .SDcols=tempCols[1:m]]
  templates[,tempMP:=paste(.SD, sep="", collapse=","), by=temp_start, .SDcols=tempCols[1:(m+1)]]
  
  return(templates)
}

.findPairs <- function(templates, r, measure="maximum")
{
  # Count templates of length m in proximity r
  tempCounts <- as.data.table(table(templates))
  temp <- sapply(tempCounts$templates, 
                  function(x){as.numeric(unlist((tstrsplit(x, ","))))})
  
  # Get distances between templates
  tempComb <- rbind(t(combn(tempCounts$templates, 2)), cbind(tempCounts$templates, 
                                                             tempCounts$templates))
  tempDist <- as.matrix(dist(t(temp), method=measure))
  tempDist <- data.table(tempComb, dist=tempDist[tempComb])
  colnames(tempDist) <- c("from", "to", "dist")
  
  # filter by radius
  tempDist <- subset(tempDist, dist<=r)

  # Merge -> remove duplicates
  tempPairs <- merge(tempDist, tempCounts, by.y=c("templates"), by.x=c("from"))
  setnames(tempPairs, "N", "N_from")
  tempPairs <- merge(tempPairs, tempCounts, by.y=c("templates"), by.x=c("to"))
  setnames(tempPairs, "N", "N_to")
  
  # Count matching pairs
  tempPairs[,matches:=fifelse(from==to, choose(N_from,2), N_from*N_to)]
  matchCount <- sum(tempPairs$matches)
  
  return(matchCount)
}

.keepSampEn <- function(x, pos, m, 
                        r=NULL, measure="maximum"){
  templates <- .getTemplates(x, pos, m)
  
  if(is.null(r)){
    B <- sum(choose(table(templates$tempM), 2))
    A <- sum(choose(table(templates$tempMP), 2))
  }
  else{
    if(r=="sd") r <- 0.2*sd(x, na.rm=TRUE)
    
    # check should give the same as above in the binary case, with r=1 & chebyshev
    B <- .findPairs(templates$tempM, r, measure)
    A <- .findPairs(templates$tempMP, r, measure)
  }
  
  # Sample Entropy
  sampEn <- -log(A/B)
  
  #TODO: implement confidence interval
  #TODO: Normalize
  
  return(sampEn)
}

sampEn <- function(data, 
                   cols=NULL,
                   block=NULL,
                   pos=NULL,
                   m=2,
                   naMode=c("remove", "keep"),
                   r=NULL,
                   measure="maximum",
                   nCores=1){
  data <- as.data.table(data)
  
  # ensure ordering
  if(is.null(pos)) pos <- 1:nrow(data)
  if(is.character(pos)) pos <- c(data[,pos])
  data$pos <- pos
  setorderv(data, c(block, "pos"))
  
  if(is.null(cols)) cols <- setdiff(colnames(data), c("pos", block))

  if(!is.null(block)){
    if(naMode=="keep"){
      scores <- data[, mclapply(.SD, .keepSampEn, pos, m, r, 
                                measure, mc.cores=nCores), 
                       by=block, .SDcols=cols] # Reduce(c,
    }
    else if(naMode=="remove")
    {
      if(is.null(r)) r <- 0
      scores <- data[, mclapply(.SD, function(x){sampleEn(na.omit(x), m, r)}, 
                                mc.cores=nCores), 
                       by=block, .SDcols=cols]
    }
  }
  else{
    if(naMode=="keep"){
      scores <- data[, mclapply(.SD, .keepSampEn, pos, m, r, 
                                measure, mc.cores=nCores),
                       .SDcols=cols]
    }
    else if(naMode=="remove")
    {
      if(is.null(r)) r <- 0
      scores <- data[, mclapply(.SD, function(x){sampleEn(na.omit(x), m, r)},
                                mc.cores=nCores), .SDcols=cols]
    }
  }
  
  setnames(scores, cols, paste("sampEn", cols, sep="_"))
  
  return(scores)
}


shannonEn <- function(data, cols=NULL, block=NULL, normalize=TRUE, 
                      discretize=TRUE, nCores=1)
{
  data <- as.data.table(data)
  if(is.null(cols)) cols <- setdiff(colnames(data), block)
  
  if(!is.null(block)){
    scores <- data[,mclapply(.SD, function(x){shannonEnDiscrete(na.omit(x), 
                                                                normalize, 
                                                                discretize)},
                           mc.cores=nCores), 
                    by=c(block), .SDcols=cols]
    }
  else{
    scores <- data[,mclapply(.SD, function(x){shannonEnDiscrete(na.omit(x), 
                                                                normalize,
                                                                discretize)},
                           mc.cores=nCores), 
                   .SDcols=cols]
  }
  
  setnames(scores, cols, paste("shannonEn", cols, sep="_"))
  
  return(scores)
}

mhl <- function(){
  
}