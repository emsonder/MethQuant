# Functionality for calculating Cross-Sample Entropies between methylation 
# patterns of cells. 
#'@author: Emanuel Sonder

library(data.table)
library(Rcpp)
library(plyr)

source("./R/binSequence.R")
source("./R/calcEntropies.R")

shift <- data.table::shift

#'@title calcCrossSampEn
#'@description Calculates Cross Sample Entropies (cross-SampEn) between methylation patterns 
#'of different cells. 
#'Sample Entropy: Richman and Moorman, 2000
#'https://journals.physiology.org/doi/full/10.1152/ajpheart.2000.278.6.h2039
#'@param metTable Either wide data.table of methylation rates.
#'dimensions: genome coordinates (pos, chr) x cells. 
#'@param cellIds unique Identifiers of single cells
#'@param templateDim Size of methylation templates (embedding dimension for 
#'sample entropy)
#'@param annotationTable Annotation set containing genomic regions of interest. 
#'Data.table needs to contain chr, start, end and name column. 
#'@param mergeAnnotations Concatenate entries of annotationTables for 
#'calculation of cross-SampEns.  
#'@return Dissimilarity matrix in form of a data.table wwith 
#' pairwise cross-SampEns between cells
calcCrossSampEn <- function(metTable, 
                            cellIds, 
                            templateDim,
                            annotationTable, 
                            minTemps=32,
                            mergeAnnotations=T){
  
  m <- templateDim
  
  # Construct Templates
  templatesSet <- list()
  for(cellId in cellIds)
  {
    # Celltable 
    columns <- c("pos", "chr", cellId)
    
    cellTable <- metTable[,..columns]
    setnames(cellTable, cellId, "rate")
    setorder(cellTable, chr, pos)
    
    # Rounding 
    cellTable[,rate:=round_any(rate, 0.5)]
    
    templatesSet[[cellId]] <- .getTemplates(cellTable, m)
  }
  templatesSet <- rbindlist(templatesSet, idcol="cell_id")

  
  if(mergeAnnotations)
  {
    annotationTable$name <- "merged"
  }
  
  setkey(annotationTable, chr, start, end)
  templatesSet <- foverlaps(templatesSet,
                            annotationTable,
                            by.x=c("chr", "temp_start", "temp_end"),
                            by.y=c("chr", "start", "end"),
                            type="within")

  templatesSubset <- subset(templatesSet, !is.na(name))
  templatesSubset[,n_temp:=.N, by=c("cell_id", "name")]
  
  
  calcPairWiseMatches <- function(tempX, tempY){
    sum(outer(tempX,tempY,function(x,y){return(x==y)}))
  }
  
  # Calculate pairwise matches and aggregate into Dissimilarity table
  disSimTable <- templatesSubset[
    , {
      cellComb <- data.table(t(combn(unique(cell_id), 2)))
      names(cellComb) <- c("cell_x", "cell_y")
      cellComb[
        , {
          if(length(tempMP[cell_id==cell_x])>minTemps & 
             length(tempMP[cell_id==cell_y]>minTemps))
          {
            data.table(B=calcPairWiseMatches(tempM[cell_id==cell_x], 
                                             tempM[cell_id==cell_y]),
                       A=calcPairWiseMatches(tempMP[cell_id==cell_x], 
                                             tempMP[cell_id==cell_y]))
          }
          else
          {
            data.table(B=as.integer(NA), A=as.integer(NA))
          }
        }
        , by = .(cell_x, cell_y)
      ]
    }
    , by = .(name)
  ]

  disSimTable[,cross_SampEn:=-log(A/B)]
  disSimTable[,comp:=paste(cell_x, cell_y, sep="-")]
  
  return(disSimTable)
}



