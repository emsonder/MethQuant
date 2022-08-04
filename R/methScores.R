# Wrapper function for quantifying methylation patterns
#'@author: Emanuel Sonder

#' @title methScores
#' Wrapper function to calculate different methylation scores.
#'
#' @description  Wrapper function to calculate different entropy scores for different dimensions
#' of methylation data (i.e. inter-cell/height axis (h) or intra-cell/width axis (w)).
#'
#' @name methScores
#' @rdname methScores
#' @param data methylation data, can be either a GRanges object with the methylated levels
#' in the meta columns, a BSseq object, or a data.table/data.frame with
#' the methylation levels of a sample in the columns, a position (posCols) and sequence name column (seqCol).
#' @param score score function to use. Scores currently implemented are sampEn and shannonEn
#' @param axis either "w" or "h". The axis "w" specifies the calculation of the scores across the (CpG) sequence of
#' a single cell/sample. The axis "h" specifies calculation of the scores across cells for a position/region of interest.
#' @param posCol Column specifying the position of a CpG if input is a data.frame/data.table.
#' @param seqCol  Column specifying the chromosome name of a CpG if input is a data.frame/data.table.
#' @param binMode Mode for binning the sequence. Mode "costum" if regions of interest are provided via
#' the regions argument. Mode "tiled" if the sequence should be split in tiled bins based on the positions.
#' Mode "fixed" if the sequence should be split up in bins containing a fixed number of CpGs.
#' Modes fixed and tiled differ in the sense that tiled bins the sequences in equally-sized stretches
#' while fixed ensures that the bins contain the same number of CpGs.
#' @param regions If binMode="costum" GRanges object which contains the regions of interest
#' @param startCol Name of the startCol (will get removed)
#' @param endCols Name of the endCol (will get removed)
#' @param binSize Size of the bins, either the length of the stretches
#' in nucleotides (binMode="tiled") or in CpGs (binMode="fixed").
#' @param ... args specific for the score function used
#' @return If data.table/data.frame or GRanges objects are provided as an input a data.table with the scores
#' for the bins across the axis of interest are provided. If a BSseq object is provided the scores get added to the
#' colData (if axis="w") or the metadata (if axis="h").
#' @export
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
