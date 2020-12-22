# Functionalities for calculation Entropies along different axes
#'@author: Emanuel Sonder

library(data.table)
library(Rcpp)

first <- data.table::first
source("./R/binSequence.R")
sourceCpp("./C++/entropies.cpp")

#'@title widthEntropy
#'@description Calculates Sample Entropies along methylation patterns
#'of genomic subsequences of single cells (width axis).
#'Sample Entropy: Richman and Moorman, 2000
#'https://journals.physiology.org/doi/full/10.1152/ajpheart.2000.278.6.h2039
#'Binary Entropy: Croll, 2013
#'https://arxiv.org/ftp/arxiv/papers/1305/1305.0954.pdf
#'@param dataTable Either wide data.table of methylation rates (byCell=F) or
#' or list with two wide data.tables, i.e. a methylation rates table and a table 
#' with cell-specific bins (byCell=T).
#'dimensions: genome coordinates (pos, chr) vs cell. 
#'@param cellIds unique Identifiers of single cells
#'@param binByCell Cell-specific binning
#'@param nCpGsBin Number of CpGs per Bin
#'@param aggregateOn Column of table defining the genomic subsequences to 
#'calculate the entropy scores on. 
#'@return data.table with Entropies per genomic subsequence and cell
widthEntropy <- function(metTable, cellIds, 
                         binByCell=T, nCpGsBin=8,
                         aggregateOn="bin"){
  # Loop proved to be faster than conversion to long
  widthEntropies <- list()
  for(cellId in cellIds)
  {
    if(binByCell){
      columns <- c("pos", "chr", cellId)}
    else{
      columns <- c("pos", "chr", cellId, aggregateOn)}
    
    # Retrieve met values for one cell
    cellTable <- metTable[complete.cases(metTable[,..cellId]),..columns]
    setnames(cellTable, cellId, "rate")
    
    # Make Sure table is ordered
    setorder(cellTable, chr, pos)
    
    nCpGs <- nrow(cellTable)
    if(binByCell) cellTable[, (aggregateOn):=fixedBinning(nCpGs, nCpGsBin)] 
    
    # Calculate Sample Entropy along width axis & aggregate
    widthEntropies[[cellId]] <- cellTable[,.(width_SampleEn=sampleEn(round(rate), 2, 0.2),
                                             width_BiEn=biEn(round(rate), tresBin=F),
                                             nCpGs_width=.N,
                                             methylation_level_cell=mean(rate),
                                             mean_dis_CpGs=(last(pos)-first(pos))/.N,
                                             median_pos=as.integer(median(pos))),
                                          by=aggregateOn]
  }
  
  widthEntropies <- rbindlist(widthEntropies,  idcol="cell_id")
  
  return(widthEntropies)
}

#'@title heightEntropy
#'@description Calculates Shannon Entropies of methylation levels across  
#'different cells (height axis). Aggregated on genomic subsequences
#'Shannon Entropy: Shannon, 1948
#'http://people.math.harvard.edu/~ctm/home/text/others/shannon/entropy/entropy.pdf
#'@param metTable Wide data.table of methylation rates, 
#'dimensions: genome coordinates (pos, chr) vs cell. 
#'@param cellIds unique Identifiers of single cells
#'@param aggregateOn Column of metTable defining the genomic subsequences to 
#'aggregate the entropy scores on. 
#'@return data.table with Shannon Entropies for the genomic subsequences
heightEntropy <- function(metTable, cellIds, aggregateOn="bin"){
  
  # Reshape to long 
  metTableLong <- melt(metTable, id.vars=c("pos", "chr", aggregateOn),
                       variable.name="cell", measure.vars=cellIds, 
                       value.name="rate")
  
  # Calculate Shannon Entropy per position 
  metTableLong[complete.cases(rate), 
               height_entropy:=shannonEnDiscrete(round(rate)), 
               by=list(chr,pos)]
  metTableLong[complete.cases(rate), nCpGs:=.N, by=list(chr,pos)]
  
  # Aggregate
  heightEntropies <- metTableLong[,.(height_entropy=mean(height_entropy, na.rm=T),
                                     mean_nCpGs_height=mean(nCpGs, na.rm=T),
                                     methylation_level_bulk=mean(rate, na.rm=T),
                                     median_pos=as.integer(median(pos))), 
                                     by=aggregateOn]
  
  return(heightEntropies)
}