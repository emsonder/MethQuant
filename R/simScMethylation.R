# Method for simulating single-Cell whole genome Bisulfite Sequencing datasets
#'@author: Emanuel Sonder

library(data.table)
library(VGAM)
library(mhsmm)

source("./R/binSequence.R")

.estimateCoveredCpGs <- function(coverageTable, estimateFromData,
                                 alpha, beta){
  if(estimateFromData)
  {
    print("Estimate Coverage from real Data")
    fit <- vglm(cbind(as.integer(n_covered), as.integer(n_uncovered)) ~ 1, 
                betabinomialff(), data = coverageTable, trace = TRUE)
    alpha <- Coef(fit)[1]
    beta <- Coef(fit)[2]
    
    print(paste("Estimated alpha:", round(alpha, 3)))
    print(paste("Estimated beta:", round(beta, 3)))
  }
  else if(is.null(alpha) & is.null(beta))
  {
    print("If parameter are not estimated from data - they have to be provided
           as an argument")
  }

  nCoveredCpGs <- rbetabinom.ab(nCells, nCpGsTotal, 
                                shape1=alpha, 
                                shape2=beta)
  
  return(nCoveredCpGs)
}


.simMetFixed <- function(nCells,
                         metTable,
                         nCoveredCpGs,
                         fullCoverage,
                         pmdProb,
                         meanAlphaNonPMD, # 0.5
                         meanAlphaPMD,
                         binSize){
  
  #TODO: -two params for beta-binomial distribution

  cols <- c("chr", "pos", "bin")
  simCells <- metTable[,..cols,]
  
  # Draw PMD status for each bin
  cellBins <- sapply(1:nCells,function(i) {paste("sim_cell_n_met_bin", i, sep="")})
  simCells[, pmd_status:=rbinom(1,1, pmdProb), by=bin]
  simCells[, alpha:=fifelse(pmd_status==1, meanAlphaPMD, meanAlphaNonPMD), by=bin]
  
  # Draw Number of methylated CpGs per bin (cell specific)
  simCells[, (cellBins):=lapply(1:nCells, function(i){
      rbetabinom.ab(1, binSizeCpGs, shape1=alpha[1], shape2=alpha[1])}), by=bin]
  
  # Draw Methylation values per bin (cell specific)
  cellIds <- sapply(1:nCells,function(i) {paste("sim_cell", i, sep="")})
  simCells[, (cellIds):=lapply(.SD, function(nMet){
    sample(rep(0:1, times=c(binSizeCpGs-nMet[1], nMet[1])), .N)}), 
    by=bin, .SDcols=cellBins]
  
  if(!fullCoverage)
  {
    nUnCoveredCpGs <- abs(nCoveredCpGs-nCpGsTotal)
    
    # randomly replace CpGs with NAs
    simCells[,(cellIds):= mapply(function(met, i){ 
      met[c(sample(1:nCpGsTotal, nUnCoveredCpGs[i]))]=NA
      return(met)}, .SD, 1:nCells, SIMPLIFY=F), 
      .SDcols=cellIds]
  }
  
  return(simCells)
}


.calcAlphas <- function(metTable, binSizeCpGs, nCpGsTotal){
  
  # CpG Binning 
  cols <- c("chr", "pos", "bin")
  cellIds <- setdiff(names(metTable), c("chr", "pos", "bin"))
  
  # Methylation statistics
  sumMet <- metTable[, ..cols]
  sumMet[, met:=rowSums(metTable[,..cellIds], na.rm=T)]
  sumMet[, total:=rowSums(!is.na(metTable[,..cellIds]), na.rm=T)]

  # Alpha Value Binning
  alphas <- seq(0.1, 3, by=0.1)
  probs <- paste("prob_", alphas, sep="")
  
  # MethylSeekr Alpha calculation
  betaLH <- function(met, total, a){(beta(met+a, total-met+a)/beta(a, a))*choose(total, met)}
  
  # Calculate Likelihood of one CpG
  sumMet[,(probs):=lapply(alphas, function(a){betaLH(met, total, a)}),
         by=seq_len(nrow(sumMet))] # probCpGs
  
  # Sliding Bin calculation -------------------------------------------------
  
  # 1. Sliding window calculation for Bins 
  
  #sumMet[,(probs):=frollapply(.SD, 101, prod, na.rm=T, align="center"), 
  #       .SDcols=probs]
  # 2. calc alpha posterior
  #sumMet[,alpha:=sum(alphas*.SD, na.rm=T)/sum(.SD, na.rm=T), 
  #        by=seq_len(nrow(sumMet)), .SDcols=probs] 
  
  #-----------------------------------------------------------------------------

  # Fixed Bin calculation 
  
  # Calculate Bin Likelihoods
  alphaTable <- sumMet[, lapply(.SD, prod, na.rm=T), 
                         by=c("bin", "chr"), .SDcols=probs]
  
  binLoc <- sumMet[,.(bin_start=min(pos),
                      bin_end=max(pos),
                      chr=unique(chr)), by=bin]
  
  alphaTable <- merge(alphaTable, binLoc, on=c("bin", "chr"),  all.x=TRUE)
  
  # Calculate Alpha Posterior Means
  alphaTable[, alpha:=sum(alphas*.SD, na.rm=T)/sum(.SD, na.rm=T), 
              .SDcols=probs, 
               by=seq_len(nrow(alphaTable))] 
  
  cols <- c("chr", "bin_start", "bin_end", "alpha", "bin")
  return(alphaTable[,..cols,])
}


# MethylSeekR Model 
.simMethylationHMM <- function(nCells,
                               metTable,
                               nCpGsTotal,
                               nCoveredCpGs,
                               fullCoverage,
                               meanAlphaNonPMD,
                               meanAlphaPMD,
                               binSizeCpGs,
                               seed)
{
  # calculate Alphas
  print("Estimating alpha parameters for bins")
  alphas <- .calcAlphas(metTable, binSizeCpGs, nCpGsTotal)
  
  # Init. HMM
  initStates <- c(0, 1);
  initTrans <- t(matrix(c(0.5, 0.5, 0.5, 0.5), nrow=2, ncol=2))
  initMu <- list(mu=c(meanAlphaNonPMD, meanAlphaPMD), sigma=c(0.1, 0.1))
  initHmm <- hmmspec(init=initStates, trans=initTrans, parms.emission=initMu, 
                      dens.emission=dnorm.hsmm)
  
  # Train HMM
  print("Training Hmm")
  train <- list(x=alphas$alpha, N=nrow(alphas))
  fittedHmm <- hmmfit(train, initHmm, mstep=mstep.norm)$model

  # Sample Alphas from HMM for each window
  print("Simulating from Hmm")
  nBins <- nrow(alphas)
  sims <- simulate.hmmspec(fittedHmm, 
                           nsim=nBins, 
                           seed=seed,
                           rand.emission=rnorm.hsmm)
  
  alphas$alpha_sim <- sims$x
  alphas$pmd_status <- sims$s
  
  alphas[, alpha_sim:=ifelse(alpha_sim<=0, 0.0001, alpha_sim), 
           by=seq_len(nrow(alphas))]
  
  cols <- c("chr", "pos", "bin")
  simCells <- metTable[,..cols,]
  simCells <- merge(simCells,  alphas, by=c("chr", "bin"), all.x=T)
  
  # Draw Number of methylated CpGs per bin (cell specific)
  cellBins <- sapply(1:nCells,function(i) {paste("sim_cell_n_met_bin", i, sep="")})
  simCells[, (cellBins):=lapply(1:nCells, function(i){
    rbetabinom.ab(1, binSizeCpGs, shape1=alpha_sim, shape2=alpha_sim)}), by=bin]
  
  # Draw Methylation values per bin (cell specific)
  cellIds <- sapply(1:nCells,function(i) {paste("sim_cell", i, sep="")})
  simCells[, (cellIds):=lapply(.SD, function(nMet){
    sample(rep(0:1, times=c(binSizeCpGs-nMet[1], nMet[1])), .N)}), 
             by=bin, .SDcols=cellBins]
  
  simCells[,(cellBins):=NULL]
  
  if(!fullCoverage)
  {
    nUnCoveredCpGs <- abs(nCoveredCpGs-nCpGsTotal)
    
    # randomly replace CpGs with NAs
    simCells[,(cellIds):= mapply(function(met, i){ 
      met[c(sample(1:nCpGsTotal, nUnCoveredCpGs[i]))]=NA
      return(met)}, .SD, 1:nCells, SIMPLIFY=F), 
            .SDcols=cellIds]
  }
  
  return(simCells)
}

#'@title simScMethylation
#'@description Simulation of single cell methylation data. Two simulation 
#'strategies can be used (hmmModel T/F). The HMM model (hmmModel=T) uses an 
#'adapted version of the MethylSeekR model and simulates from it. 
#' https://academic.oup.com/nar/article/41/16/e155/2411535
#'The plain model (hmmModel=F) samples methylation values for each bin 
#'independently from a beta-binomial distribution parameterized by the arguments 
#'meanAlphaNonPMD and meanAlphaPMD.
#'@param nCells Number of cells to simulate methylation data for.
#'@param metTable Real data for parameter estimation. Wide data.table of 
#'methylation rates, dimensions: genome coordinates (pos, chr) vs cell. 
#'@param binSizCpgs Number of CpGs per bin. Each bin gets parametrized.
#'@param meanAlphaNonPMD Shape parameter (alpha) of beta-binomial to sample number
#' of methylated CpGs per Bin for non PMD regions. Only required for plain model.
#'@param meanAlphaPMD Shape parameter (alpha) of beta-binomial to sample number
#' of methylated CpGs per Bin for PMD regions. Only required for plain model.
#'@param hmmModel Use MethylSeekR model (hmmModel=T) or plain model (hmmModel=F)
#'@param seed Seed for random number generator
#'@param pmdProb Probability of a bin belonging to a PMD region. 
#'@param estimateFromData estimate coverage from data
#'@param fullCoverage All CpGs covered for all cells.
#'@param alpha Shape parameter (alpha) for beta-binomial for CpG coverage sampling.
#'Only required if estimateFromData=F.
#'@param beta Shape parameter (beta) for beta-binoial for CpG coverage sampling.
#'Only required if estimateFromData=F.
#'@param subSample Subsample from real input data for speed-up
#'@param subSampleSize Number of cells to sample from real data for parameter
#'estimation.
simScMethylation <- function(nCells, 
                             metTable=NULL,
                             binSizeCpGs=128,
                             meanAlphaNonPMD=0.4, # 0.5
                             meanAlphaPMD=1.3,
                             hmmModel=T,
                             seed=43, 
                             pmdProb=0.7, #  can be estimated on MethylSeekR output
                             estimateFromData=T,
                             fullCoverage=F,
                             # manually inject coverage rate?
                             alpha=NULL,
                             beta=NULL,
                             subSample=T,
                             subSampleSize=100){
  
  set.seed(seed)
  
  # TODO:
  # - simulate number of reads per CpG
  # => beta binomial by cpgsxbin (x reads covering cpgs in bin) methylated
  # - Compare coverage between pmds and non-pmds! (Save plots)
  # - Rename meanAlphaNonPMD/meanAlphaPMD 
  
  # Subsample on the full data
  cellIds <- setdiff(names(metTable), c("chr", "pos")) # or cellIds as an input
  
  if(subSample)
  {
    print("Subsample real Data for Speed-Up")
    cellIds <- cellIds[sample(1:length(cellIds), subSampleSize)]
    cols <- c("chr", "pos", cellIds)
    metTable <- metTable[,..cols]
    
    metTable <- metTable[rowSums(!is.na(metTable[,..cellIds]))>0,]
  }
  
  nCpGsTotal <- nrow(metTable)
  setorder(metTable, chr, pos)

  # Number of CpG-Sites covered per cell
  coverageTable <- metTable[, lapply(.SD, function(x){sum(!is.na(x))}),  
                             .SDcols=cellIds]
  coverageTable <- data.table(cellId=names(coverageTable), 
                              transpose(coverageTable))
  setnames(coverageTable, c("V1"), c("n_covered"))
  
  coverageTable[, n_uncovered:=nCpGsTotal-n_covered]
  
  # Sample Number of covered CpGs by Cell
  # TODO: Implement coverage sampling better for non-data estimated sims
  nCoveredCpGs <- .estimateCoveredCpGs(coverageTable, 
                                       estimateFromData,
                                       alpha,
                                       beta)
  
  # bin sequence 
  metTable$bin <- fixedBinning(nCpGsTotal, binSizeCpGs)
  
  if(hmmModel)
  {
    cat("Simulating single cell methylation data by an HMM model",
        "can take up a few minutes\n")
    simCells <- .simMethylationHMM(nCells,
                                   metTable,
                                   nCpGsTotal,
                                   nCoveredCpGs,
                                   fullCoverage,
                                   meanAlphaNonPMD,
                                   meanAlphaPMD,
                                   binSizeCpGs,
                                   seed)
  }
  else
  {
    simCells <- .simMetFixed(nCells,
                             metTable,
                             nCoveredCpGs,
                             fullCoverage,
                             meanAlphaNonPMD,
                             meanAlphaPMD,
                             binSizeCpGs,
                             seed)
  }
  
  return(simCells)
}
