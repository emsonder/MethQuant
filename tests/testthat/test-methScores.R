# Unit tests for the methScores wrapper function
#'@author Emanuel Sonder
#'
#' Test here rather test for correct output dimensions regarding 
#' different parameter combinationsthan correctness of binning & entropy scores, 
#' for these see test-binning.R & test-scores.R.

# Setup
bsData <- readRDS("/mnt/plger/esonder/R/MethQuant/smp/tests/testObjects/metScBs.rds")
bsDataFrame <- readRDS("/mnt/plger/esonder/R/MethQuant/smp/tests/testObjects/metDfSmall.rds")
regions <- GRanges(seqnames="chrSim",
                   ranges=IRanges(start=c(100, 500, 10000),
                                  end=c(200, 600, 12000)))

# ShannonEn Height -------------------------------------------------------------

test_that("costum binning ShannonEn height", {
  
  scores <- methScores(bsData, 
                       score="shannonEn", 
                       axis="h", 
                       binMode="costum", 
                       regions=regions, 
                       normalize=TRUE)
  
  expect_equal(nrow(scores[[2]]), length(regions))
})

test_that("costum binning ShannonEn height - bsseq metadata", {
  
  scores <- methScores(bsData, 
                       score="shannonEn", 
                       axis="h", 
                       binMode="costum", 
                       regions=regions, 
                       normalize=TRUE)
  
  expect_equal(nrow(metadata(scores[[1]])$entropyScores), length(regions))
})

test_that("tiled binning ShannonEn height", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData,regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="shannonEn", 
                       axis="h", 
                       binMode="tiled", 
                       binSize=binSize,
                       normalize=TRUE)
  
  nBin <- ceiling(width(regions[1])/binSize)
  expect_equal(nrow(scores[[2]]), nBin)
})

test_that("fixed binning ShannonEn height", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData,regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="shannonEn", 
                       axis="h", 
                       binMode="fixed", 
                       binSize=binSize,
                       normalize=TRUE)
  
  nBin <- ceiling(nrow(bsSubData)/binSize)
  expect_equal(nrow(scores[[2]]), nBin)
})

# ShannonEn width --------------------------------------------------------------

test_that("costum binning ShannonEn width", {
  
  scores <- methScores(bsData, 
                       score="shannonEn", 
                       axis="w", 
                       binMode="costum", 
                       regions=regions, 
                       normalize=TRUE)
  
  expect_equal(nrow(scores[[2]]), length(regions)*ncol(bsData))
})

test_that("costum binning ShannonEn width - bsseq metadata", {
  
  scores <- methScores(bsData, 
                       score="shannonEn", 
                       axis="w", 
                       binMode="costum", 
                       regions=regions, 
                       normalize=TRUE)

  expect_equal(nrow(colData(scores[[1]])), ncol(bsData), 
               label="correct number of rows (i.e. samples)")
  expect_equal(ncol(colData(scores[[1]])), length(regions), 
               label="correct number of cols (i.e. bins)")
})

test_that("tiled binning ShannonEn width", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData, regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="shannonEn", 
                       axis="w", 
                       binMode="tiled", 
                       binSize=binSize,
                       normalize=TRUE)
  
  nBins <- ceiling(width(regions[1])/binSize)
  
  expect_equal(nrow(colData(scores[[1]])), ncol(bsSubData), 
               label="correct number of rows (i.e. samples)")
  expect_equal(ncol(colData(scores[[1]])), nBins, 
               label="correct number of cols (i.e. bins)")
})

test_that("fixed binning ShannonEn width", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData, regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="shannonEn", 
                       axis="w", 
                       binMode="fixed", 
                       binSize=binSize,
                       normalize=TRUE)
  
  nBins <- ceiling(nrow(bsSubData)/binSize)
  
  expect_equal(nrow(colData(scores[[1]])), ncol(bsSubData), 
               label="correct number of rows (i.e. samples)")
  expect_equal(ncol(colData(scores[[1]])), nBins, 
               label="correct number of cols (i.e. bins)")
})

# SampEn Height ----------------------------------------------------------------

test_that("costum binning SampEn (keep) height", {
  
  scores <- methScores(bsData, 
                       score="sampEn", 
                       axis="h", 
                       binMode="costum", 
                       regions=regions, 
                       naMode="keep", 
                       m=2)
  
  expect_equal(nrow(scores[[2]]), sum(width(regions)))
})

test_that("costum binning SampEn (remove) height", {
  
  scores <- methScores(bsData, 
                       score="sampEn", 
                       axis="h", 
                       binMode="costum", 
                       regions=regions, 
                       naMode="remove", 
                       m=2)
  
  expect_equal(nrow(scores[[2]]), sum(width(regions)))
})

test_that("costum binning SampEn (keep) height - bsseq metadata", {
  
  scores <- methScores(bsData, 
                       score="sampEn", 
                       axis="h", 
                       binMode="costum", 
                       regions=regions, 
                       naMode="keep", 
                       m=2)
  
  expect_equal(nrow(metadata(scores[[1]])$entropyScores), sum(width(regions)))
})

test_that("costum binning SampEn (remove) height - bsseq metadata", {
  
  scores <- methScores(bsData, 
                       score="sampEn", 
                       axis="h", 
                       binMode="costum", 
                       regions=regions, 
                       naMode="remove", 
                       m=2)
  
  expect_equal(nrow(metadata(scores[[1]])$entropyScores), sum(width(regions)))
})

test_that("tiled binning SampEn (keep) height", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData,regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="sampEn", 
                       axis="h", 
                       binMode="tiled", 
                       binSize=binSize, 
                       naMode="keep", 
                       m=2)
  
  expect_equal(nrow(scores[[2]]), width(regions[1]))
})

test_that("tiled binning SampEn (remove) height", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData,regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="sampEn", 
                       axis="h", 
                       binMode="tiled", 
                       binSize=binSize, 
                       naMode="remove", 
                       m=2)
  
  expect_equal(nrow(scores[[2]]), width(regions[1]))
})

test_that("fixed binning SampEn (keep) height", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData,regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="sampEn", 
                       axis="h", 
                       binMode="fixed", 
                       binSize=binSize, 
                       naMode="keep", 
                       m=2, 
                       nCores=1)
  
  expect_equal(nrow(scores[[2]]), nrow(bsSubData))
})

test_that("fixed binning SampEn (remove) height", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData,regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="sampEn", 
                       axis="h", 
                       binMode="fixed", 
                       binSize=binSize, 
                       naMode="remove", 
                       m=2, 
                       nCores=1)
  
  expect_equal(nrow(scores[[2]]), nrow(bsSubData))
})

# SampEn Width -----------------------------------------------------------------

test_that("costum binning SampEn (keep) width", {
  
  scores <- methScores(bsData, 
                       score="sampEn", 
                       axis="w", 
                       binMode="costum", 
                       regions=regions, 
                       naMode="keep", 
                       m=2)
  
  expect_equal(nrow(scores[[2]]), length(regions)*ncol(bsData))
})

test_that("costum binning SampEn (remove) width", {
  
  scores <- methScores(bsData, 
                       score="sampEn", 
                       axis="w", 
                       binMode="costum", 
                       regions=regions, 
                       naMode="remove", 
                       m=2)
  
  expect_equal(nrow(scores[[2]]), length(regions)*ncol(bsData))
})

test_that("costum binning SampEn (keep) width - bsseq metadata", {

  scores <- methScores(bsData, 
                       score="sampEn", 
                       axis="w", 
                       binMode="costum", 
                       regions=regions, 
                       binSize=1000, 
                       naMode="keep", 
                       m=2, 
                       nCores=1)
  
  expect_equal(nrow(colData(scores[[1]])), ncol(bsData), 
               label="correct number of rows (i.e. samples)")
  expect_equal(ncol(colData(scores[[1]])), length(regions), 
               label="correct number of cols (i.e. bins)")
})

test_that("costum binning SampEn (remove) width - bsseq metadata", {
  
  scores <- methScores(bsData, 
                       score="sampEn", 
                       axis="w", 
                       binMode="costum", 
                       regions=regions, 
                       binSize=1000, 
                       naMode="remove", 
                       m=2, 
                       nCores=1)
  
  expect_equal(nrow(colData(scores[[1]])), ncol(bsData), 
               label="correct number of rows (i.e. samples)")
  expect_equal(ncol(colData(scores[[1]])), length(regions), 
               label="correct number of cols (i.e. bins)")
})

test_that("tiled binning SampEn (keep) width", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData, regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="sampEn", 
                       axis="w", 
                       binMode="tiled", 
                       binSize=binSize, 
                       naMode="keep", 
                       m=2, 
                       nCores=1)
  
  nBins <- ceiling(width(regions[1])/binSize)
  
  expect_equal(nrow(colData(scores[[1]])), ncol(bsSubData), 
               label="correct number of rows (i.e. samples)")
  expect_equal(ncol(colData(scores[[1]])), nBins, 
               label="correct number of cols (i.e. bins)")
})

test_that("tiled binning SampEn (remove) width", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData, regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="sampEn", 
                       axis="w", 
                       binMode="tiled", 
                       binSize=binSize, 
                       naMode="remove", 
                       m=2, 
                       nCores=1)
  
  nBins <- ceiling(width(regions[1])/binSize)
  
  expect_equal(nrow(colData(scores[[1]])), ncol(bsSubData), 
               label="correct number of rows (i.e. samples)")
  expect_equal(ncol(colData(scores[[1]])), nBins, 
               label="correct number of cols (i.e. bins)")
})

test_that("fixed binning SampEn (keep) width", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData, regions[1])
  binSize <- 50

  scores <- methScores(bsSubData, 
                       score="sampEn", 
                       axis="w", 
                       binMode="fixed", 
                       binSize=binSize, 
                       naMode="keep", 
                       m=2, 
                       nCores=1)
  
  nBins <- ceiling(nrow(bsSubData)/binSize)
  
  expect_equal(nrow(colData(scores[[1]])), ncol(bsSubData), 
               label="correct number of rows (i.e. samples)")
  expect_equal(ncol(colData(scores[[1]])), nBins, 
               label="correct number of cols (i.e. bins)")
})

test_that("fixed binning SampEn (remove) width", {
  
  # Subset for fast testing
  bsSubData <- subsetByOverlaps(bsData, regions[1])
  binSize <- 50
  
  scores <- methScores(bsSubData, 
                       score="sampEn", 
                       axis="w", 
                       binMode="fixed", 
                       binSize=binSize, 
                       naMode="remove", 
                       m=2, 
                       nCores=1)
  
  nBins <- ceiling(nrow(bsSubData)/binSize)
  
  expect_equal(nrow(colData(scores[[1]])), ncol(bsSubData), 
               label="correct number of rows (i.e. samples)")
  expect_equal(ncol(colData(scores[[1]])), nBins, 
               label="correct number of cols (i.e. bins)")
})

# data.table/data.frame inputs -------------------------------------------------

test_that("costum binning SampEn (keep) width", {
  
  cellIds <- setdiff(colnames(bsDataFrame), c("pos", "chr"))
  scores <- methScores(bsDataFrame, 
                       score="sampEn", 
                       axis="w", 
                       binMode="costum", 
                       regions=regions, 
                       naMode="keep", 
                       m=2)
  
  expect_equal(nrow(scores), length(regions)*length(cellIds))
})

test_that("costum binning ShannonEn height", {
  
  cellIds <- setdiff(colnames(bsDataFrame), c("pos", "chr"))
  scores <- methScores(bsDataFrame, 
                       score="shannonEn", 
                       axis="h", 
                       binMode="costum", 
                       regions=regions, 
                       normalize=TRUE)
  
  expect_equal(nrow(scores), length(regions))
})
