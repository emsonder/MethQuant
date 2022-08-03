source("/mnt/plger/esonder/R/MethQuant/smp/R/methScores.R")

library(data.table)
library(plyr)
library(Rcpp)
library(BiocParallel)

data <- readRDS("/mnt/plger/esonder/R/MethQuant/smp/tests/testObjects/metScBs.rds")
regions <- GRanges(seqnames="chrSim",
                   ranges=IRanges(start=c(100, 500, 10000),
                                  end=c(200, 600, 12000)))

# Shannon Height
methScores(data, score="shannonEn", axis="h", binMode="costum", regions=regions, normalize=TRUE)
methScores(data, score="shannonEn", axis="h", binMode="tiled", binSize=1000, normalize=TRUE)
methScores(data, score="shannonEn", axis="h", binMode="fixed", binSize=1000, normalize=TRUE)

# Shannon Width
methScores(data, score="shannonEn", axis="w", binMode="costum", regions=regions, normalize=TRUE)
methScores(data, score="shannonEn", axis="w", binMode="tiled", binSize=1000, normalize=TRUE)
methScores(data, score="shannonEn", axis="w", binMode="fixed", binSize=1000, normalize=TRUE)

# SampEn (discrete) height
methScores(data, score="sampEn", axis="h", binMode="costum", regions=regions, naMode="keep", m=2)
methScores(data, score="sampEn", axis="h", binMode="tiled", binSize=1000, naMode="keep", m=2, nCores=8)
methScores(data, score="sampEn", axis="h", binMode="fixed", binSize=1000, naMode="keep", m=2, nCores=1)

methScores(data, score="sampEn", axis="h", binMode="costum", regions=regions, naMode="remove", m=2)
methScores(data, score="sampEn", axis="h", binMode="tiled", binSize=1000, naMode="remove", m=2, nCores=4)
methScores(data, score="sampEn", axis="h", binMode="fixed", binSize=1000, naMode="remove", m=2, nCores=8)

# SampEn (discrete) width
methScores(data, score="sampEn", axis="w", binMode="costum", regions=regions, naMode="keep", m=2)
methScores(data, score="sampEn", axis="w", binMode="tiled", binSize=1000, naMode="keep", m=2, nCores=1)
methScores(data, score="sampEn", axis="w", binMode="fixed", binSize=1000, naMode="keep", m=2, nCores=1)

methScores(data, score="sampEn", axis="w", binMode="costum", regions=regions, naMode="remove", m=2)
methScores(data, score="sampEn", axis="w", binMode="tiled", binSize=1000, naMode="remove", m=2, nCores=1)
methScores(data, score="sampEn", axis="w", binMode="fixed", binSize=1000, naMode="remove", m=2, nCores=1)