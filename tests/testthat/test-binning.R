test_that("binning tiled", {
  
  x <- data.table(pos=c(101:130, 161:200, 1:100),
                  chr=c(rep("chr1", 70), rep("chr2", 100)))
  x[, bin:=.binning(max(pos), 50, pos, chr, mode="tiled"), by=c("chr")]
  
  expect_equal(length(unique(x$bin)), 4, label="correct number of bins")
  expect_equal(as.character(unique(subset(x, pos>161 & pos<170)$bin)), 
               "(chr1,151,200]")
  
  binsExpected <- c("(chr1,101,150]", "(chr1,151,200]",
                    "(chr2,1,50]", "(chr2,51,100]")
  binsExpected <- factor(binsExpected, levels=binsExpected)
  expect_setequal(unique(x$bin), binsExpected)
})

# case where cpgs the same as positions

test_that("binning fixed", {
  
  x <- data.table(pos=c(101:130, 161:200, 1:100),
                  chr=c(rep("chr1", 70), rep("chr2", 100)))
  x[, bin:=.binning(.N, 50, pos, chr, mode="fixed"), by=c("chr")]
  
  expect_equal(length(unique(x$bin)), 4, label="correct number of bins")
  expect_equal(as.character(unique(subset(x, pos>161 & pos<170)$bin)), 
               "(chr1,101,180]")
  
  binsExpected <- c("(chr1,101,180]", "(chr1,181,200]",
                    "(chr2,1,50]", "(chr2,51,100]")
  binsExpected <- factor(binsExpected, levels=binsExpected)
  expect_setequal(unique(x$bin), binsExpected)
})
