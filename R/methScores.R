# Wrapper function for qu
#'@author: Emanuel Sonder

#TODO: 
# - ShannonEn height axis scores vs bins
# - data naming & variable naming 
# - cell type /sample type column



methScores <- function(data, # can be GRanges, or bsseq, or data.table / data.frame
                       score=c("sampEn", "shannonEn"), # one of the implemented scores ? but how to ensure correctness ? 
                       axis=c("w", "h"),
                       posCol="pos",
                       seqCol="chr",
                       binMode=c("costum", "tiled", "fixed"),
                       regions=NULL,
                       startCol="start",
                       endCol="end",
                       binSize=1000,
                       ...){
  
    # check input arguments
    if(!(score %in% c("sampEn", "shannonEn"))){
      stop("Please choose one of the following scores: sampEn, shannonEn") # maybe add a constants file
    }
  
    # Retrieve the methylation matrix
    if(is(data, "BSseq")){
      methData <- bsseq::getMeth(data, type="raw", what="perBase")
      methData <- as.data.table(methData)
      methData <- cbind(methData,
                        as.data.table(granges(data))[,c("seqnames", "start")])
      setnames(methData, c("seqnames", "start"), c(seqCol, posCol))
    }
    else if(is(data, "GRanges")){
      methData <- as.data.table(data)[,c("seqnames", "start", 
                                         names(mcols(data))), with=FALSE]
      setnames(methData, c("seqnames", "start"), c(seqCol, posCol))
    }
    else if(is(data, "data.frame") || is(data, "data.table")){
      methData <- as.data.table(data)
      }
    else{
      stop("Data class must be one of the following: 
            BSseq, GRanges, data.frame or data.table (preferred)")
    }
  
    cellIds <- setdiff(colnames(methData), c(seqCol, posCol))
    methData[,c(startCol, endCol):=get(posCol)]
    
    # Binning
    if(binMode=="costum")
    {
      regionsDt <- as.data.table(regions)
      setnames(regionsDt, c("seqnames", "start", "end"), c(seqCol, startCol, endCol), 
               skip_absent=TRUE)
      regionsDt[,bin:=paste0("(", chr, ",", start, ",", end, "]")]
      regionsDt$bin <- factor(regionsDt$bin, levels=unique(regionsDt$bin))
      setkeyv(regionsDt, c(seqCol, startCol, endCol))
      
      methData <- foverlaps(methData, 
                            regionsDt[,c(seqCol, startCol, endCol, "bin"), 
                                    with=FALSE],
                            by.x=c(seqCol, startCol, endCol),
                            by.y=c(seqCol, startCol, endCol),
                            type="within", 
                            nomatch=NULL)
    }
    else if(binMode=="tiled")
    {
      #methData[, bin:=cut(pos, seq(min(get(posCol)), 
      #                             max(get(posCol))+binSize, binSize), 
      #                    include.lowest=TRUE), by=seqCol]
      methData[,maxPos:=max(get(posCol)), by=seqCol]
      methData[, bin:=.binning(max(get(posCol)), binSize, 
                              get(posCol), get(seqCol),
                              mode="tiled"), 
                 by=seqCol]
      methData$maxPos <- NULL
    }
    else if(binMode=="fixed")
    {
      methData[, bin:=.binning(.N, binSize, 
                              get(posCol), get(seqCol),
                              mode="fixed"), 
                 by=seqCol]
    }
    else if(binMode=="moving")
    {
      stop("Moving bin mode not yet implemented")
    }

    methData <- methData[,c(seqCol, posCol, cellIds, "bin"), with=FALSE]
    #TODO: Ensure ordering

    scoreFun <- get(score) # here use keyword args
    if(axis=="h")
    {
      methData <- melt(methData, id.vars=c(seqCol, "bin", posCol))
      setnames(methData, c("variable", "value"), c("cell_id", "met_level"))
    
      # could also do specifically not with posCol
      if(score=="sampEn") block <- c(seqCol, "bin", posCol) else block <- c(seqCol, "bin")
      scores <- scoreFun(methData, cols=c("met_level"), 
                         block=block, ...)
      
      #scores <- shannonEn(methData, cols=c("met_level"), block=c(seqCol, "bin", posCol))
      setnames(scores, c(paste(score, "met_level", sep="_")), 
                       c(paste(score, axis, sep="_")))
    }
    else if(axis=="w")
    {
      scores <- scoreFun(methData, cols=cellIds, 
                         block=c(seqCol, "bin"), ...)
     # scores <- shannonEn(methData, cols=cellIds, block=c(seqCol, "bin"))
      scores <- melt(scores, id.vars=c(seqCol, "bin"), 
                     value.name=paste(score, axis, sep="_"))
      
      scores[,sample_id:=tstrsplit(variable, split="_", keep=2)]
      scores$variable <- NULL
    }
    else
    {
      stop("Axis needs to be either w (width, intra-cell) or h (height, inter-cell)")
    }
    
    # add as a slot if bsseq or granges
    if(is(data, "BSseq"))
    {
      if(axis=="h")
      {
        if(binMode=="costum")
        {
          data <- subsetByOverlaps(data, regions)
        }
        
        metadata(data) <- append(metadata(data), list("entropyScores"=scores))
      }
      else if(axis=="w")
      {
        scoreCol <- setdiff(names(scores), c(seqCol, "bin", "sample_id"))
        scoresWide <- dcast(scores, bin~sample_id, value.var=scoreCol)
        colNames <- paste(score, scoresWide$bin, sep="_")
        scoresWide$bin <- NULL
        
        scoresWide <- DataFrame(t(scoresWide))
        colnames(scoresWide) <- colNames

        colData(data) <- cbind(colData(data), scoresWide)
      }
      
      return(list(data, scores))
    }
    
    return(scores)
}
