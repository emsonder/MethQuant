# Create Test objects



#bsData <- readRDS("./tests/testObjects/metScBs.rds")
#bsDataFrame <- readRDS("./tests/testObjects/metDfSmall.rds")
#regions <- GRanges(seqnames="chrSim",
#                   ranges=IRanges(start=c(100, 500, 10000),
#                                  end=c(200, 600, 12000)))

# # use the sim function
#
# # produce different sizes:
# # sample sizes: small (10 samples), large (1000 samples)
# # coverage: low, high
# setwd("/mnt/plger/esonder/R/MethQuant/smp/tests/testObjects")
#
# transMat1 <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow=2, ncol=2)
# metDataSmall <- simMetPattern(nCpGs=1e6,
#                               nCells=10,
#                               mode="markov",
#                               covParams=c(0.2,0.3),
#                               states=c('0', '1.0'),
#                               estimateTransMat=FALSE,
#                               transMat=transMat1,
#                               estimateRate=FALSE,
#                               metData=NULL,
#                               sizeParam=1,
#                               muParam=5.86)
# cellIds <- setdiff(colnames(metDataSmall), c("pos", "chr"))
#
# # data.table -------------------------------------------------------------------
#
# metDtSmall <- as.data.table(metDataSmall)
# saveRDS(metDtSmall, "metDtSmall.rds")
#
# # data.frame -------------------------------------------------------------------
#
# metDfSmall <- as.data.frame(metDataSmall)
# saveRDS(metDfSmall, "metDfSmall.rds")
#
# # GRanges ----------------------------------------------------------------------
#
# metRanges <- GRanges(seqnames=metDataSmall$chr,
#                      ranges=IRanges(start=metDataSmall$pos,
#                                     end=metDataSmall$pos),
#                      met_level=metDataSmall[,cellIds, with=FALSE])
# saveRDS(metRanges, "metRanges.rds")
#
# # bsseq object: sc -------------------------------------------------------------
#
# metCov <- metDataSmall[, lapply(.SD,function(col){fifelse(is.na(col), 0,1)}),
#                        .SDcols=cellIds]
# metLevel <- metDataSmall[, lapply(.SD,function(col){fifelse(is.na(col), 0, col)}),
#                          .SDcols=cellIds]
# metScBs <-   BSseq(chr=metDataSmall$chr,
#                    pos=metDataSmall$pos,
#                    M=as.matrix(metLevel),
#                    Cov=as.matrix(metCov),
#                    sampleNames=cellIds)
# saveRDS(metScBs, "metScBs.rds")
#
# # bsseq object: bulk -----------------------------------------------------------
#
# # full coverage
# # some nonesensical samples
# nMet <- replicate(10, sample(1:20, size=1e6, replace=TRUE, prob=NULL))
# nUnMet <- replicate(10, sample(1:20, size=1e6, replace=TRUE, prob=NULL))
# cov <- nMet+nUnMet
# M <- nMet/cov
#
# metBs <-   BSseq(chr=metDataSmall$chr,
#                  pos=metDataSmall$pos,
#                  M=M,
#                  Cov=cov,
#                  sampleNames=cellIds)
# #saveRDS(metBs, "metBs.rds")
