#source("./R/calcCpGDensity.R")

# Tiling Binning approach (user input bin size in bps)
tileBinning <- function(dt, binSize)
{
  # Bin by fixed window sizes in bps
  dt[,bin:=cut(pos, seq(min(pos), max(pos)+binSize, binSize), 
               include.lowest=TRUE), by=chr]
  
}

# Fixed Binning approach (fixed number of CpGs per bin)
.binning <- function(nTot, nBin, pos, chr, mode)
{
  coords <- data.table(start=pos, end=pos)
  nIntervals <- ceiling(nTot/nBin)
  bins <- rep(1:nIntervals, each=nBin)
  
  #bins <- sapply(1:nIntervals, function(i){rep(i, nBin)})
  
  if(mode=="tiled")
  {
    bins <- data.table(start=1:nTot, end=1:nTot, bin=bins)
    bins[,bin_start:=min(start), by=c("bin")]
    bins[,bin_end:=max(start), by=c("bin")]
    
    setkey(bins, bin_start, bin_end)
    bins <- foverlaps(coords, 
                      bins, 
                      by.x=c("start", "end"),
                      by.y=c("bin_start", "bin_end"),
                      type="within",nomatch=NULL, mult="last")
    bins <- bins[,c("bin_start", "bin_end", "bin"), with=FALSE]
    #setnames(bins, "start", "pos")
  }
  else if(mode=="fixed")
  {
    bins <- data.table(bin=bins[1:nTot], pos=pos)
    bins[,bin_start:=min(pos), by=c("bin")]
    bins[,bin_end:=max(pos), by=c("bin")]
  }
  
  bins[, bin:=paste0("(", chr, ",", bin_start,",", bin_end, "]")]
  bins$bin <- factor(bins$bin, levels=unique(bins$bin))
  
  # tiled binning
  #methData[, bin:=cut(pos, seq(min(get(posCol)), 
  #                             max(get(posCol))+binSize, binSize), 
  #                    include.lowest=TRUE), by=seqCol]

  return(bins$bin)
}

# Distance based Binning 
distanceBinning <- function()
{
  error("Not Implemented yet")
}


# Density based Binning
#'@param byCell determine bins by cell (TRUE) or on bulk data 
densityBinning <- function(metTable, byCell=TRUE, mode="linear", 
                           windowSize=1000, densityTolerance=0.001)
{
  #cellIds <- setdiff(names(metTable), c("chr", "pos"))
  #cpgCell <- calcCpGDensity(cpgCell, mode="linear", windowSize=1000)

  error("Not Implemented yet")
  
  # add decay options
}
