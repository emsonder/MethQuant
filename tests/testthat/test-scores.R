test_that("adjacent template construction", {

  x1 <- c(1,NA,0,NA,0,0,1,NA,NA,0,0,0,0,0)
  x2 <- c(1, NA, NA)

  temps2 <- .getTemplates(x1, pos=1:length(x1), m=2)
  temps3 <- .getTemplates(x1, pos=1:length(x1), m=3)
  tempsNone <- .getTemplates(x2, pos=1:length(x2), m=2)

  expect_equal(nrow(temps2), 4)
  expect_equal(nrow(temps3), 2)
  expect_equal(nrow(tempsNone), 0)
})

test_that("find matching pairs Chebyshev",{

  x1 <- c(1,NA,0,NA,0,0,1,NA,NA,0,0,0,0,0)
  temps <- .getTemplates(x1, pos=1:length(x1), m=2)

  pairsM <- .findPairs(temps$tempM, r=1, measure="maximum")
  pairsMP <- .findPairs(temps$tempMP, r=1, measure="maximum")

  expect_equal(pairsM, choose(4,2))
  expect_equal(pairsMP, choose(4,2))
})

test_that("find matching pairs Euclidean",{

  x1 <- c(1,NA,0,NA,0,0,0.5,NA,NA,0,0,0.75,0,0)
  temps <- .getTemplates(x1, pos=1:length(x1), m=2)

  pairsM <- .findPairs(temps$tempM, 1, measure="euclidean")
  pairsMP <- .findPairs(temps$tempMP, 1, measure="euclidean")

  expect_equal(pairsM, 5)
  expect_equal(pairsMP, 3)
})

test_that("Equivalence of NA modes for no NAs (discrete)",{

  x <- data.table(V1=sample(c(0,1), 100, replace=TRUE),
                  pos=1:100)

  yKeep <- sampEn(x, naMode="keep")
  yRemove <- sampEn(x, naMode="remove")

  expect_equal(as.numeric(yKeep), as.numeric(yRemove))
})

test_that("Equivalence of NA modes for no NAs (continuous)",{

  set.seed(42)
  x <- data.table(V1=runif(100, 0, 5),
                  pos=1:100)
  r <- sd(x$V1)*0.2

  yKeep <- sampEn(x, r=r, measure="maximum", naMode="keep")
  yRemove <- sampEn(x, r=r, measure="maximum", naMode="remove")

  expect_equal(as.numeric(yKeep), as.numeric(yRemove))
})

test_that("Sample Entropy - consistency with other implementations",{

  set.seed(42)
  x <- data.table(V1=runif(100, 0, 5),
                  pos=1:100)
  r <- sd(x$V1)*0.2

  # Results from other implementations, commented out to save dependencies
  # library(pracma)
  # yRef <- sample_entropy(x$V1, edim = 2, r = r, tau = 1)

  yRef <- 1.852384
  yKeep <- sampEn(x, m=2, r=r, measure="maximum", naMode="keep")
  yRemove <- sampEn(x, m=2, r=r, measure="maximum", naMode="remove")

  expect_equal(as.numeric(yKeep), yRef, tolerance=1e-03)
  expect_equal(as.numeric(yRemove), yRef, tolerance=1e-03)
})

test_that("Shannon Entropy - consistency with other implementations (discrete)",{

  set.seed(42)
  x <- data.table(V1=sample(c(0,1), 100, replace=TRUE))

  # Results from other implementations, commented out to save dependencies
  # library(entropy)
  # yRef <- entropy(c(table(x$V1)), method="ML", unit="log")
  yRef <- 0.6859298
  y <- shannonEn(x, normalize=FALSE)

  expect_equal(as.numeric(y), yRef, tolerance=1e-03)
})

# test_that("Shannon Entropy - consistency with other implementations (continuous, binning)",{
#
#   set.seed(42)
#   x <- data.table(V1=runif(1000,0,2))
#
#   # Results from other implementations, commented out to save dependencies
#   # library(entropy)
#   bins  <- entropy::discretize(c(x$V1), 40, r=range(c(x$V1)))
#   # yRef <- entropy(bins, method="ML", unit="log")
#   yRef <- 3.664798
#
#   binVars <- unlist(lapply(1:40, function(i){rep(i, a[i])}))
#   y <- shannonEn(binVars, normalize=FALSE, discretize=TRUE)
#
#   expect_equal(as.numeric(y), yRef, tolerance=1e-03)
# })

test_that("Shannon Entropy - consistency with other implementations (continuous)",{

  set.seed(42)
  x <- data.table(V1=runif(1000,0,2))

  # Results from other implementations, commented out to save dependencies
  # library(entropy)
  # yRef <- entropy(table(x$V1), method="ML", unit="log")
  yRef <- 6.907755
  y <- shannonEn(x$V1, normalize=FALSE, discretize=FALSE)

  expect_equal(as.numeric(y), yRef, tolerance=1e-03)
})

