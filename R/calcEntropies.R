# Functionalities for calculation Entropies along different axes
#'@author: Emanuel Sonder

library(data.table)
library(plyr)
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
#'@param dataTable Either wide data.table of methylation rates.
#'dimensions: genome coordinates (pos, chr) x cells. 
#'@param cellIds unique Identifiers of single cells
#'@param binByCell Cell-specific binning
#'@param nCpGsBin Number of CpGs per Bin
#'@param templateDim Size of methylation templates (embedding dimension for 
#'sample entropy)
#'@param aggregateOn Column of table defining the genomic subsequences to 
#'calculate the entropy scores on. 
#'@return data.table with Entropies per genomic subsequence and cell
widthEntropy <- function(metTable, cellIds, 
                         binByCell=T, nCpGsBin=128,
                         templateDim=2,
                         aggregateOn="bin", scaleFactors=c(2,4,8)){
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
    
    # Make sure table is ordered
    setorder(cellTable, chr, pos)
    
    nCpGs <- nrow(cellTable)
    if(binByCell) cellTable[, (aggregateOn):=fixedBinning(nCpGs, nCpGsBin)] 
    
    # Possible methylation rates: 0, 1, 0.5 (allele-specific methylation)
    cellTable[,rate:=round_any(rate, 0.5)]
    #cellTable[,rate:=fifelse(rate==1 | rate==0, rate, 0.5)]
    
    # Calculate Sample Entropy along width axis & aggregate
    widthEntropies[[cellId]] <- cellTable[,.(width_SampleEn=sampleEn(rate, templateDim, 0),
                                             width_BiEn=biEn(rate, tresBin=F),
                                             width_MsEn_sc1=msEn(rate, scaleFactors[1]), # Make it generic
                                             width_MsEn_sc2=msEn(rate, scaleFactors[2]),
                                             width_MsEn_sc3=msEn(rate, scaleFactors[3]),
                                             nCpGs_width=.N,
                                             bin_start = min(pos),
                                             bin_end = max(pos),
                                             methylation_level_cell=mean(rate),
                                             mean_dis_CpGs=(last(pos)-first(pos))/.N,
                                             median_pos=as.integer(median(pos))),
                                            by=c(aggregateOn, "chr")]
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
  
  # Calculate Shannon Entropy per position & aggregateOn in case there are 
  # several annotations per position
  cols <- c("chr", "pos", aggregateOn)
  metTableLong[complete.cases(rate), 
               height_entropy:=shannonEnDiscrete(round_any(rate, 0.5)), 
               by=cols]
  
  
  metTableLong[complete.cases(rate), nCpGs:=.N, by=cols]
  
  # Aggregate
  if(is.null(aggregateOn)) aggregateOn <- c("chr", "pos")
  
  heightEntropies <- metTableLong[,.(height_entropy=mean(height_entropy, na.rm=T),
                                     mean_nCpGs_height=mean(nCpGs, na.rm=T),
                                     methylation_level_bulk=mean(rate, na.rm=T),
                                     median_pos=as.integer(median(pos))), 
                                     by=aggregateOn]
  
  return(heightEntropies)
}