# Functionalities for calculation Entropies along different axes
#'@author: Emanuel Sonder

library(data.table)
library(plyr)
library(Rcpp)

first <- data.table::first
source("./R/binSequence.R")
sourceCpp("./C++/entropies.cpp")
sourceCpp("./C++/mhl.cpp")

.getTemplates <- function(cellTable, m){
  
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
  templates[,tempM:=paste(.SD, sep="", collapse=""), by=temp_start, .SDcols=tempCols[1:m]]
  templates[,tempMP:=paste(.SD, sep="", collapse=""), by=temp_start, .SDcols=tempCols[1:(m+1)]]
  
  return(templates)
}

#'@title widthKeepSampEn
#'@description Calculates Sample Entropies along methylation patterns
#'of genomic subsequences of single cells (width axis). Slight adaption of Sample
#'Entropy used here to only consider adjacent CpGs for templates, termed KeepSampEn.
#'Sample Entropy: Richman and Moorman, 2000
#'https://journals.physiology.org/doi/full/10.1152/ajpheart.2000.278.6.h2039
#'KeepSampEn: Dong et al., 2019
#'https://www.mdpi.com/1099-4300/21/3/274
#'@param metTable Either wide data.table of methylation rates.
#'dimensions: genome coordinates (pos, chr) x cells. 
#'@param cellIds unique Identifiers of single cells
#'@param nTempsBin Number of templates per Bin. Need for cond. probability 
#'estimation. Should not be too small.
#'@param templateDim Size of methylation templates (embedding dimension for 
#'sample entropy)
#'@param tiledBinning If TRUE tiled (non-overlapping) bins are used across the templates, 
#'if FALSE moving bins (one template at the time) are applied. 
#'@return data.table with Entropies per genomic subsequence and cell
widthKeepSampEn <- function(metTable, 
                            cellIds, 
                            nTempsBin=32,
                            templateDim=2,
                            tiledBinning=T){
  
  m <- templateDim
  
  # Loop proved to be faster than conversion to long
  widthEntropies <- list()
  for(cellId in cellIds)
  {
    columns <- c("pos", "chr", cellId)
    
    cellTable <- metTable[,..columns]
    setnames(cellTable, cellId, "rate")
    setorder(cellTable, chr, pos)
    
    # Rounding 
    cellTable[,rate:=round_any(rate, 0.5)]
    
    # Construct templates
    templates <- .getTemplates(cellTable, m)
    
    if(tiledBinning)
    {
      # Binning
      nTemps <- nrow(templates)
      templates[, bin:=fixedBinning(nTemps, nTempsBin)] 
    
      # count pairwise matches for both options
      templates[,B:=sum(choose(table(tempM), 2)), by=bin]
      templates[,A:=sum(choose(table(tempMP), 2)), by=bin]
      
      # Sample Entropy
      templates[,width_sampEn:=-log(A/B), by=bin]
      
      # CI Calculation
      templates[,sd_matches:=sd(c(rep(1,A[1]),rep(0,B[1]-A[1]))),by=bin]
      templates[,width_SampEn_CI_margin:=sd_matches*qt(0.975,df=B-1)/sqrt(B),
                by=bin]
      templates[,bin_width:=max(temp_end)-min(temp_start), by=bin]
      
      # Calculate methylation rate for Bins
      bins <- templates[,.(bin_start=min(temp_start),
                           bin_end=max(temp_end),
                           chr=unique(chr)), by=bin]
      
      setkey(bins, chr, bin_start, bin_end)
      cellTable[,start_pos:=pos]
      cellTable[,end_pos:=pos]
      metLevelBin <- foverlaps(cellTable, bins,
                               by.x=c("chr", "start_pos", "end_pos"),
                               by.y=c("chr", "bin_start", "bin_end"),
                               type="any")
      metLevelBin <- metLevelBin[!is.na(bin),
                                 .(met_level_width_bin=mean(rate, na.rm=T),
                                   chr=unique(chr)),
                                 by=bin]
      
      # Calculate methylation rate for templates
      setkey(templates, chr, temp_start, temp_end)
      
      # first to prevent that CpGs are counted more than once (as templates move by one position)
      metLevelTemps <- foverlaps(cellTable, templates[,c("chr", "temp_start", 
                                                         "temp_end", "bin"), 
                                                      with=F],
                                 by.x=c("chr", "start_pos", "end_pos"),
                                 by.y=c("chr", "temp_start", "temp_end"),
                                 mult="first",
                                 type="any",
                                 nomatch=NULL)
      metLevelTemps <- metLevelTemps[!is.na(bin),
                                     .(met_level_width_temp=mean(rate, na.rm=T),
                                       nCpGs_temps=.N,
                                       chr=unique(chr)),
                                     by=bin]
      
      # Add methylation rates to templates
      templates <- merge(templates, metLevelBin, 
                         by.x=c("chr", "bin"),
                         by.y=c("chr", "bin"),
                         all.x=T, all.y=F)
      templates <- merge(templates, metLevelTemps, 
                         by.x=c("chr", "bin"),
                         by.y=c("chr", "bin"),
                         all.x=T, all.y=F)
      
      widthEntropies[[cellId]] <- templates[,.(width_SampEn=unique(width_sampEn),
                                               width_SampEn_CI_margin=unique(width_SampEn_CI_margin),
                                               bin_width=unique(bin_width),
                                               bin_start=min(temp_start),
                                               bin_end=max(temp_end),
                                               met_level_width_bin=unique(met_level_width_bin),
                                               met_level_width_temp=unique(met_level_width_temp),
                                               nCpGs_temps=unique(nCpGs_temps),
                                               A=unique(A),
                                               B=unique(B),
                                               chr=unique(chr)), by=bin]
    }
    else 
    {
      # convert templates
      templates[,tempM_fac:=as.numeric(factor(tempM))]
      templates[,tempMP_fac:=as.numeric(factor(tempMP))]
      
      # Calculate Matches in rolling windows 
      templates[,B:=frollapply(tempM_fac, nTempsBin/2, 
                               function(temp){sum(choose(table(temp), 2))},
                               align="center")]
      templates[,A:=frollapply(tempMP_fac, nTempsBin/2, 
                               function(temp){sum(choose(table(temp), 2))},
                               align="center")]
      templates[,bin_start:=shift(temp_start, nTempsBin/2, type="lag")]
      templates[,bin_end:=shift(temp_end, nTempsBin/2, type="lead")]
      templates[,bin:=pos]
      
      # Sample Entropy
      templates[,width_sampEn:=-log(A/B)]
      
      # Calculate methylation rate for Bins
      bins <- templates[!is.na(bin_start) & !is.na(bin_end),
                        .(bin_start=bin_start,
                          bin_end=bin_end,
                          chr=unique(chr)), 
                        by=bin]
      
      setkey(bins, chr, bin_start, bin_end)
      cellTable[,start_pos:=pos]
      cellTable[,end_pos:=pos]
      metLevelBin <- foverlaps(cellTable, bins,
                               by.x=c("chr", "start_pos", "end_pos"),
                               by.y=c("chr", "bin_start", "bin_end"))
      metLevelBin <- metLevelBin[!is.na(bin),
                                 .(met_level_width_bin=mean(rate, na.rm=T),
                                   chr=unique(chr)),
                                 by=bin]
      
      templates <- merge(templates, metLevelBin, 
                         by.x=c("chr", "bin"),
                         by.y=c("chr", "bin"),
                         all.x=T, all.y=F)
      
      # Calculate methylation rate for templates
      
      # biased! & refactor!
      templates[,met_level_width_temp:=(get("0")+get("1")+get("2"))/3, 
                by=seq_len(nrow(templates))]
      templates[,met_level_width_temp:=frollapply(met_level_width_temp, 
                                                  nTempsBin/2, 
                                                  function(met_level_width_temp){
                                                    mean(met_level_width_temp)},
                                                  align="center")]
      
      # CI Calculation
      templates <- subset(templates, !is.na(A) & !is.na(B))
      templates[,sd_matches:=sd(c(rep(1,A),rep(0,B-A))),
                 by=seq_len(nrow(templates))]
      templates[,width_SampEn_CI_margin:=sd_matches*qt(0.975,df=B-1)/sqrt(B),
                 by=seq_len(nrow(templates))]
      templates[,bin_width:=max(temp_end)-min(temp_start)]
      
      widthEntropies[[cellId]] <- templates[,c("chr", "pos", 
                                               "B", "A", 
                                               "width_sampEn", 
                                               "width_SampEn_CI_margin",
                                               "bin_width", 
                                               "met_level_width_bin",
                                               "met_level_width_temp")]
    }
  }
  
  widthEntropies <- rbindlist(widthEntropies,  idcol="cell_id")
  
  return(widthEntropies)
}

#'@title widthMHL
#'@description Calculates MHL score along methylation patterns
#'of genomic subsequences of single cells (width axis). Instead of reads it is
#'calculated along strips of CpGs not intersected by missing values in the data.
#'MHL score: Guo et al., 2017
#'https://pubmed.ncbi.nlm.nih.gov/28263317/
#'@param metTable Either wide data.table of methylation rates.
#'dimensions: genome coordinates (pos, chr) x cells. 
#'@param cellIds unique Identifiers of single cells
#'@return data.table with MHL scores per genomic subsequence and cell.
widthMHL <- function(metTable, cellIds)
{
  
  # Loop proved to be faster than conversion to long
  widthMHLScores <- list()
  for(cellId in cellIds)
  {
    columns <- c("pos", "chr", cellId)
    
    cellTable <- metTable[,..columns]
    setnames(cellTable, cellId, "rate")
    setorder(cellTable, chr, pos)
    
    # Rounding 
    cellTable[,rate:=round_any(rate, 0.5)]
    
    fill <- tidyr::fill
    # Find all adjacent Templates of arbitrary length
    cellTable[,isEnd:=fifelse(is.na(shift(rate,n=1,type="lead")) & !is.na(rate), TRUE, FALSE)]
    cellTable[,tempId:=fifelse(isEnd, pos, as.integer(NA))]
    cellTable[,tempId:=nafill(tempId, type="nocb")]
    
    cellTable <- cellTable[!is.na(rate),]
    cellTable[,N_adj:=.N, by=tempId]
    cellTable <- cellTable[N_adj>=2,]
    
    
    cellTable[,width_mhl:=mhl(rate), by=tempId]
    
    widthMHLScores[[cellId]] <- cellTable[,.(width_mhl=mhl(rate),
                                             nCpGs_width=.N,
                                             bin_start = min(pos),
                                             bin_end = max(pos),
                                             met_level_width=mean(rate),
                                             mean_dis_CpGs=(last(pos)-first(pos))/.N,
                                             median_pos=as.integer(median(pos))),
                                          by=c("tempId", "chr")] # by chr
  }
  
  widthMHLScores <- rbindlist(widthMHLScores,  idcol="cell_id")
  
  return(widthMHLScores)
}


#'@title heightShannonEn
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
heightShannonEn <- function(metTable, cellIds, aggregateOn="bin"){
  
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
  
  # Normalize Height Entropy
  metTableLong[,height_entropy:=(1/log2(nCpGs+1))*height_entropy]
  
  # Aggregate
  if(is.null(aggregateOn)) aggregateOn <- c("chr", "pos")
  
  heightEntropies <- metTableLong[,.(height_entropy=mean(height_entropy, na.rm=T),
                                     mean_nCpGs_height=mean(nCpGs, na.rm=T),
                                     met_level_height=mean(rate, na.rm=T),
                                     median_pos=as.integer(median(pos))), 
                                  by=aggregateOn]
  
  return(heightEntropies)
}

#'@title widthShannonEn
#'@description Calculates Shannon Entropies along methylation patterns
#'of genomic subsequences of single cells (width axis).
#'Shannon Entropy: Shannon, 1948
#'http://people.math.harvard.edu/~ctm/home/text/others/shannon/entropy/entropy.pdf
#'@param metTable Either wide data.table of methylation rates.
#'dimensions: genome coordinates (pos, chr) x cells. 
#'@param cellIds unique Identifiers of single cells
#'@param nCpGsBin Number of CpGs per Bin
#'@param aggregateOn Column of table defining the genomic subsequences to 
#'calculate the entropy scores on. 
#'@return data.table with entropies per genomic subsequence and cell
widthShannonEn <- function(metTable, cellIds, 
                           binByCell=T, nCpGsBin=15,
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
    
    # Make sure table is ordered
    setorder(cellTable, chr, pos)
    
    nCpGs <- nrow(cellTable)
    if(binByCell) cellTable[, (aggregateOn):=fixedBinning(nCpGs, nCpGsBin)] 
    
    # Possible methylation rates: 0, 1, 0.5 (allele-specific methylation)
    cellTable[,rate:=round_any(rate, 0.5)]
    cellTable[,methylation_level_width:=mean(rate), by=c(aggregateOn, "chr")]
    
    # Calculate Sample Entropy along width axis & aggregate
    # TODO: Check ordering!r
    widthEntropies[[cellId]] <- cellTable[,.(met_level_width=mean(rate),
                                             width_ShannonEn=shannonEnDiscrete(rate),
                                             nCpGs_width=.N,
                                             bin_start=min(pos),
                                             bin_end=max(pos),
                                             mean_dis_CpGs=(last(pos)-first(pos))/.N,
                                             median_pos=as.integer(median(pos))),
                                          by=c(aggregateOn, "chr")]
  }
  
  widthEntropies <- rbindlist(widthEntropies,  idcol="cell_id")
  
  return(widthEntropies)
}

