# Functionalities for calculation Entropies along different axes
#'@author: Emanuel Sonder

library(data.table)
library(Rcpp)

first <- data.table::first
sourceCpp("../C++/entropies.cpp")


#'@title widthEntropy
#'@description Calculates Sample Entropies along methylation patterns
#'of genomic subsequences of single cells (width axis).
#'Sample Entropy: Richman and Moorman, 2000
#'https://journals.physiology.org/doi/full/10.1152/ajpheart.2000.278.6.h2039
#'@param metTable Wide data.table of methylation rates, 
#'dimensions: genome coordinates (pos, chr) vs cell. 
#'@param cellIds unique Identifiers of single cells
#'@param aggregateOn Column of metTable defining the genomic subsequences to 
#'calculate the entropy scores on. 
#'@return data.table with Sample Entropies per genomic subsequence and cell
widthEntropy <- function(metTable, cellIds, aggregateOn="bin"){

  # Loop proved to be faster than conversion to long
  widthEntropies <- list()
  for(cellId in cellIds)
  {
    # Retrieve met values for one cell
    columns <- c("pos", "chr", cellId, aggregateOn)
    cellTable <- metTable[complete.cases(metTable[,..cellId]),..columns]
    setnames(cellTable, cellId, "rate")
    
    # Calculate Sample Entropy along width axis & aggregate
    widthEntropies[[cellId]] <- cellTable[,.(width_SampleEn=sampleEn(round(rate), 2, 0.2),
                                             nCpGs_width=.N,
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
  metTableLong[complete.cases(rate), height_entropy:=shannonEnDiscrete(round(rate)), by=pos]
  metTableLong[complete.cases(rate), nCpGs:=.N, by=pos]
  
  # Aggregate
  heightEntropies <- metTableLong[,.(height_entropy=mean(height_entropy, na.rm=T),
                                     mean_nCpGs_height=mean(nCpGs, na.rm=T),
                                     median_pos=as.integer(median(pos))), 
                                     by=aggregateOn]
  
  return(heightEntropies)
}