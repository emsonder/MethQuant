#source("./R/calcCpGDensity.R")

# Tiling Binning approach (user input bin size in bps)
tileBinning <- function(dt, binSize)
{
  # Bin by fixed window sizes in bps
  dt[,bin:=cut(pos, seq(min(pos), max(pos)+binSize, binSize), 
               include.lowest=TRUE), by=chr]
  
}

# Fixed Binning approach (fixed number of CpGs per bin)
fixedBinning <- function(nCpGs, nCpGsBin, pos, chr)
{
  nIntervals <- ceiling(nCpGs/nCpGsBin)
  
  bins <- sapply(1:nIntervals, function(i){rep(i, nCpGsBin)})
  bins <- bins[1:nCpGs]
  
  bins <- data.table(bin=bins, pos=pos, chr)

  bins[,bin_start:=min(pos), by=c("bin", "chr")]
  bins[,bin_end:=max(pos), by=c("bin", "chr")]
  bins[, bin:=paste0("(", chr, ",", bin_start,",", bin_end, "]")]
  bins$bin <- factor(bins$bin, levels=unique(bins$bin))

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
