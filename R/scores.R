#' Quantification functions
#' @author Emanuel Sonder

library(wCorr)
library(fuzzyjoin)

.getTemplates <- function(x, pos, m){

  cellTable <- data.table(rate=x, pos=pos)

  # Mark adjacent CpGs
  isCovered <- !is.na(cellTable[,rate])
  temps <- frollsum(isCovered, m+1, align="right")
  cellTable$isAdj <- fifelse(temps==m+1 & !is.na(temps), T, F)

  # Construct (adjacent) templates
  cellTable[,temp_start:=fifelse(data.table::shift(isAdj,n=m, type="lead"), pos, NaN, na=NaN)]
  cellTable[,temp_end:=fifelse(!is.na(temp_start), data.table::shift(pos, n=m,type="lead"), NaN)] # temp_start+m

  tempCols <- as.character(0:m)
  cellTable[,(tempCols):=lapply(0:m, function(i) data.table::shift(rate, n=i, type="lead"))]

  templates <- cellTable[!is.nan(temp_start),]
  templates[,tempM:=paste(.SD, sep="", collapse=","), by=temp_start, .SDcols=tempCols[1:m]]
  templates[,tempMP:=paste(.SD, sep="", collapse=","), by=temp_start, .SDcols=tempCols[1:(m+1)]]

  return(templates)
}

.findPairs <- function(templates, r, measure="maximum", weighDist=FALSE)
{
  # Count templates of length m in proximity r
  tempCounts <- as.data.table(table(templates))
  temp <- sapply(tempCounts$templates,
                  function(x){as.numeric(unlist((tstrsplit(x, ","))))})

  # Get distances between templates
  if(length(tempCounts$templates)==1){
    tempComb <- cbind(tempCounts$templates, tempCounts$templates)}
  else
  {
    tempComb <- rbind(t(combn(tempCounts$templates, 2)), cbind(tempCounts$templates,
                                                               tempCounts$templates))
  }

  tempDist <- as.matrix(dist(t(temp), method=measure))
  tempDist <- data.table(tempComb, dist=tempDist[tempComb])
  colnames(tempDist) <- c("from", "to", "dist")

  # filter by radius
  tempDist <- subset(tempDist, dist<=r)

  # Merge -> remove duplicates
  tempPairs <- merge(tempDist, tempCounts, by.y=c("templates"), by.x=c("from"))
  setnames(tempPairs, "N", "N_from")
  tempPairs <- merge(tempPairs, tempCounts, by.y=c("templates"), by.x=c("to"))
  setnames(tempPairs, "N", "N_to")

  # Count matching pairs
  tempPairs[,matches:=fifelse(from==to, choose(N_from,2), N_from*N_to)]
  matchCount <- sum(tempPairs$matches)

  return(matchCount)
}

# thats much slower
.findWeightedPairs <- function(templates, r=NULL, measure){

  templates[,id:=1:nrow(templates)]

  if(is.null(r))
  {
    matchM <- merge(templates, templates,
                    by.x="tempM", by.y="tempM",
                    allow.cartesian=TRUE) # slow
    matchM <- subset(matchM, id.x!=id.y) # exlcude self-matches
    B <- sum(matchM$w.x*matchM$w.y)
    matchMP <- merge(templates, templates,
                     by.x="tempMP", by.y="tempMP",
                     allow.cartesian=TRUE) # slow
    matchMP <- subset(matchMP, id.x!=id.y) # exclude self-matches
    A <- sum(matchMP$w.x*matchMP$w.y)
  }
  else
  {
    measure <- match.arg(measure, choices=c("euclidean", "manhattan"))
    # get distance between m pairs
    matchM <- fuzzyjoin::distance_full_join(templates,
                                            templates,
                                            by=c("0", "1"),
                                            method=measure)
    matchM <- subset(matchM, id.x!=id.y)
    B <- sum(matchM$w.y*matchM$w.x)

    # get distance between m+1 pairs
    matchMP <- fuzzyjoin::distance_full_join(templates,
                                             templates,
                                             by=c("0", "1", "2"))
    matchMP <- subset(matchMP, id.x!=id.y)
    A <- sum(matchMP$w.y*matchMP$w.x)
  }

  return(list(A=A, B=B))
}

.keepSampEn <- function(x,
                        pos,
                        m,
                        r=NULL,
                        measure="euclidean",
                        weighDist=FALSE, ...){
  templates <- .getTemplates(x, pos, m)

  if(weighDist)
  {
    templates[,dist:=temp_end-temp_start]
    templates[,w:=.weightDecay(dist, isDist=TRUE, ...)]
  }

  if(is.null(r)){
    if(weighDist)
    {
      matches <- .findWeightedPairs(templates, r=r, measure=measure)
      B <- matches$B
      A <- matches$A
    }
    else
    {
      B <- sum(choose(table(templates$tempM), 2))
      A <- sum(choose(table(templates$tempMP), 2))
    }
  }
  else{
    if(r=="sd") r <- 0.2*sd(x, na.rm=TRUE)

    if(weighDist)
    {
      matches <- .findWeightedPairs(templates, r=r, measure=measure)
      B <- matches$B
      A <- matches$A
    }
    else{
      # check should give the same as above in the binary case, with r=1 & chebyshev
      B <- .findPairs(templates$tempM, r, measure)
      A <- .findPairs(templates$tempMP, r, measure)
    }
  }

  # Sample Entropy
  sampEn <- -log(A/B)

  #TODO: implement confidence interval
  #TODO: Normalize

  return(sampEn)
}
#' Sample Entropy function
#'
#' Implementation of Sample Entropy (Richman and Moorman, 2000,
#' https://journals.physiology.org/doi/full/10.1152/ajpheart.2000.278.6.h2039).
#' The complexity measure Sample Entropy is obtained by determining the negative
#' logarithm of the number of matching subseries/templates of length m+1 (A) and length m (B).
#' Matching subseries/templates are within a distance of r each other.
#'
#' @name sampEn
#' @rdname sampEn
#' @param data data for which to calculate Sample Entropies. Format expected is
#' a data.table with the numerical data for which to calculate the Sample Entropies
#' in the columns.
#' @param cols vector with names of the columns for which to calculate the Sample Entropy.
#' @param block If the Sample Entropies should be calculated for subseries of the
#' original data, block is the name of the column(s) which indicate the subseries.
#' @param pos Numerical column which gives the positions of the elements (not necessarily needed, might be removed)
#' @param m Length of the subseries/templates
#' @param naMode How to treat NAs. Mode "remove" all NAs get removed from the series and Sample Entropy
#' is calculated for the series with removed NAs. Mode "keep" subseries/templates get constructed such that no
#' subseries/templates is intersected by a NA.
#' (KeepSampEn, Dong et al. 2019, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7514754/)
#' @param r distance threshold (tolerance) for which subseries/templates are counted as matching pairs.
#' If r=NULL, r is set to 0, thus only exact matches are considered.
#' Makes sense if the data has only few states (e.g. binary series).
#' @param measure distance measure to be used, all distance measures of the dist function can be used:
#' "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski.
#' @param nCores, number of cores to be used, parallelization happens across the cols.
#' @return data.table with entropy scores for the specified columns and subseries
#' @export
sampEn <- function(data,
                   cols=NULL,
                   block=NULL,
                   pos=NULL,
                   weighDist=FALSE,
                   m=2,
                   naMode=c("keep", "remove"),
                   r=NULL,
                   measure="maximum",
                   nCores=1, getBg=FALSE, longFormat=FALSE, ...){

  data <- as.data.table(data)
  naMode <- match.arg(naMode, choices=c("keep", "remove"))

  # ensure ordering
  if(is.null(pos)) pos <- 1:nrow(data)

  if(is.character(pos)){
    pos <- c(data[,pos])
    data$pos <- pos}
  setorderv(data, union(block, "pos"))

  if(is.null(cols)) cols <- setdiff(colnames(data), c("pos", block))

  if(!is.null(block)){
    if(naMode=="keep"){
      scores <- data[, mclapply(.SD, .keepSampEn, pos, m, r,
                                measure, weighDist, ..., mc.cores=nCores),
                       by=block, .SDcols=cols] # Reduce(c,
    }
    else if(naMode=="remove")
    {
      if(is.null(r)) r <- 0
      scores <- data[, mclapply(.SD, function(x){sampleEn(na.omit(x), m, r)},
                                mc.cores=nCores),
                       by=block, .SDcols=cols]
    }
  }
  else{
    if(naMode=="keep"){
      scores <- data[, mclapply(.SD, .keepSampEn, pos, m, r,
                                measure, mc.cores=nCores),
                       .SDcols=cols]
    }
    else if(naMode=="remove")
    {
      if(is.null(r)) r <- 0
      scores <- data[, mclapply(.SD, function(x){sampleEn(na.omit(x), m, r)},
                                mc.cores=nCores), .SDcols=cols]
    }
  }

  if(getBg)
  {
    scores <- .getBackground(data, "sampEn",
                             scores,
                             cols=cols,
                             block=block,
                             nCores=nCores,
                             ...)
  }
  else if(longFormat & !getBg)
  {
    # switch to long
    scores <- melt(scores, id.vars=block,
                   value.name="sampEn")
    setorderv(scores, cols=block)
  }

  return(scores)
}

#' @export
metLevel <- function(data,
                     cols=NULL,
                     block=NULL,
                     pos=NULL,
                     nCores=1,
                     longFormat=FALSE,
                     getBg=FALSE, ...){

  if(!is.null(block))
  {
    scores <- data[, mclapply(.SD, function(x){mean(x, na.rm=TRUE)},
                            mc.cores=nCores), by=block, .SDcols=cols]
  }
  else
  {
    scores <- data[, mclapply(.SD, function(x){mean(x, na.rm=TRUE)},
                              mc.cores=nCores), .SDcols=cols]
  }

  if(getBg)
  {
    longFormat <- TRUE
  }

  if(longFormat)
  {
    # switch to long
    scores <- melt(scores, id.vars=block,
                   value.name="metLevel")
    setorderv(scores, cols=block)
  }

  return(scores)
}

#' Shannon Entropy function
#'
#' Implementation of Shannon Entropy
#' (Shannon, 1948, http://people.math.harvard.edu/~ctm/home/text/others/shannon/entropy/entropy.pdf)
#'
#' @name shannonEn
#' @rdname shannonEn
#' @param data data for which to calculate Shannon Entropies. Format expected is
#' a data.table with the numerical data for which to calculate the Shannon Entropies
#' in the columns.
#' @param cols vector with names of the columns for which to calculate the Sample Entropy.
#' @param block If the Shannon Entropies should be calculated for subseries of the
#' original data, block is the name of the column(s) which indicate the subseries.
#' @param normalize (TRUE/FALSE) If the shannon entropy should be normalized by the number of states.
#' @param discretize Should the values be discretized to the next integer
#' @param nCores, number of cores to be used, parallelization happens across the cols.
#' @return data.table with entropy scores for the specified columns and subseries
#' @export
shannonEn <- function(data, cols=NULL, block=NULL, normalize=TRUE,
                      discretize=TRUE, nCores=1,
                      getBg=FALSE, longFormat=FALSE, ...)
{
  data <- as.data.table(data)
  if(is.null(cols)) cols <- setdiff(colnames(data), block)

  if(!is.null(block)){
    scores <- data[,mclapply(.SD, function(x){shannonEnDiscrete(x[!is.na(x)],
                                                                normalize,
                                                                discretize)},
                           mc.cores=nCores),
                    by=c(block), .SDcols=cols]
    }
  else{
    scores <- data[,mclapply(.SD, function(x){shannonEnDiscrete(x[!is.na(x)],
                                                                normalize,
                                                                discretize)},
                           mc.cores=nCores),
                   .SDcols=cols]
  }

  if(getBg)
  {
    scores <- .getBackground(data, "shannonEn",
                             scores,
                             cols=cols,
                             block=block,
                             nCores=nCores,
                             ...)
  }
  else if(longFormat & !getBg)
  {
    # switch to long
    scores <- melt(scores, id.vars=block,
                   value.name="shannonEn")
    setorderv(scores, cols=block)
  }

  return(scores)
}

.weightDecay <- function(pos,
                         decay=c("exp", "linear"),
                         range=200,
                         maxDist=1000,
                         shift=1,
                         minDist=NULL,
                         isDist=FALSE, ...){

  decay <- match.arg(decay, choices=c("exp", "linear"))

  if(!isDist)
  {
    d <- data.table::shift(pos, type="lead", n=shift)-pos
  }
  else
  {
    d <- pos
  }

  if(decay=="exp")
  {
    w <- exp(-d/range)
  }
  else if(decay=="linear")
  {
    w <- (1-d/maxDist)
  }

  w[w<0] <-0

  if(!is.null(minDist))
  {
    w[d<minDist] <- 0
  }

  return(w)
}

ent <- function(p, norm, ...){
  if(norm)
  {
    etr <- sum(-p*log(p), na.rm=TRUE)/log(length(p))
  }
  else
  {
    etr <- sum(-p*log(p), na.rm=TRUE)
  }
  return(etr)
}

.trProb <- function(x, pos,
                    weighDist=TRUE,
                    type=c("ctcm", "w"),
                    w=NULL,
                    aggFun=sum,
                    norm=TRUE){

  type <- match.arg(type, choices=c("ctcm", "w"))
  x <- round(as.numeric(x))

  xs <- data.table::shift(x, type="lead", n=1)
  if(weighDist)
  {
    if(type=="w")
    {
      dt <- data.table(x=x, xs=xs, w=w)
      dt <- dt[complete.cases(dt),]
      tr <- dt[,.(w_t=sum(w)),by=.(x,xs)]
      tr <- dt[,.(p_t=w_t/sum(w_t)),by=.(x)]
    }
    else
    {
      d <- data.table::shift(pos, type="lead", n=1)-pos
      dt <- data.table(x=x, xs=xs, dist=d)
      dt <- dt[complete.cases(dt),]

      tr <- dt[,.(r_t=.N/sum(dist)),by=.(x,xs)]
      tr[,r_tot:=sum(r_t), by=x]
      tr[,p_t:=r_t/r_tot]
    }
  }
  else
  {
    dt <- data.table(x=x, xs=xs)
    dt <- dt[complete.cases(dt),]

    tr <- dt[,.(n_t=.N), by=.(x,xs)]
    tr[,p_t:=n_t/sum(n_t), by=x]
  }

  tr <- matrix(c(subset(tr, x==1 & xs==1)$p_t[1],
                 subset(tr, x==0 & xs==1)$p_t[1],
                 subset(tr, x==1 & xs==0)$p_t[1],
                 subset(tr, x==0 & xs==0)$p_t[1]), nrow=2, ncol=2)

  score <- aggFun(apply(tr, 1, ent, norm), na.rm=TRUE)
  score
}

#' @export
trProb <- function(data,
                   cols=NULL,
                   block=NULL,
                   weighDist=TRUE,
                   type=c("ctcm", "w"),
                   nCores=1,
                   getBg=FALSE,
                   longFormat=TRUE,
                   aggFun=sum,
                   norm=TRUE, ...){

  if(is.null(cols)) cols <- setdiff(colnames(data), block)

  if(weighDist){
    pos <- data[["pos"]]
    w <- .weightDecay(pos, ...)}

  if(!is.null(block)){
    scores <- data[,mclapply(.SD, function(x){.trProb(x,
                                                      pos,
                                                      weighDist=weighDist,
                                                      type=type,
                                                      w=w,
                                                      aggFun=aggFun)},
                             mc.cores=nCores),
                   by=c(block), .SDcols=cols]
  }
  else{
    scores <- data[,mclapply(.SD, function(x){.trProb(x,
                                                      pos,
                                                      weighDist=weighDist,
                                                      type=type,
                                                      w=w,
                                                      aggFun=aggFun)},
                             mc.cores=nCores),
                   .SDcols=cols]
  }

  scores[,(cols):=lapply(.SD, as.numeric), .SDcols=cols]

  if(getBg)
  {
    scores <- .getBackground(data, "trProb",
                             scores,
                             cols=cols,
                             block=block,
                             nCores=nCores,
                             ...)
  }
  else if(longFormat & !getBg)
  {
    # switch to long
    scores <- melt(scores, id.vars=block,
                   value.name="trProb")
    setorderv(scores, cols=block)
  }

  return(scores)
}

# weighted correlation for bulk data
#' @export
corr <-  function(data,
                  cols=NULL,
                  block=NULL,
                  weighDist=TRUE,
                  nCores=1,
                  shift=1,
                  getBg=FALSE,
                  longFormat=FALSE, ...){

  if(weighDist){
    pos <- data[["pos"]]
    w <- .weightDecay(pos, ...)}

  if(!is.null(block)){
    scores <- data[,mclapply(.SD,function(x,w){
      xs <- data.table::shift(x, type="lead", n=shift)
      if(weighDist)
      {
        dt <- data.table(x=x, xs=xs, w=w)
        dt <- dt[complete.cases(dt),]
        s <- weightedCorr(dt$x, dt$xs, weights=dt$w, method="Pearson")
      }
      else
      {
        dt <- data.table(x=x, xs=xs)
        dt <- dt[complete.cases(dt),]
        s <- weightedCorr(dt$x, dt$xs, method="Pearson")
      }
      s}, w, mc.cores=nCores), by=c(block), .SDcols=cols]}
  else
  {
    scores <- data[,mclapply(.SD, function(x,w){
      xs <- data.table::shift(x, type="lead", n=shift)
      if(weighDist)
      {
        dt <- data.table(x=x, xs=xs)
        dt <- dt[complete.cases(dt),]
        s <- weightedCorr(dt$x, dt$xs, weights=w)
      }
      else
      {
        dt <- data.table(x=x, xs=xs)
        dt <- dt[complete.cases(dt),]
        s <- weightedCorr(dt$x, dt$xs)
      }
      s}, w, mc.cores=nCores), .SDcols=cols]
  }

  if(getBg)
  {
    scores <- .getBackground(data, "corr",
                             scores,
                             cols=cols,
                             block=block,
                             nCores=nCores,
                             ...)
  }
  else if(longFormat & !getBg)
  {
    # switch to long
    scores <- melt(scores, id.vars=block,
                   value.name="corr")
    setorderv(scores, cols=block)
  }

  return(scores)
}


.smc <- function(x, xs, w=NULL){
  if(is.null(w))
  {
    w <- rep(1, length(x))
  }
  if(sum(!is.na(x) & !is.na(xs))>10 & (sum(w>0, na.rm=TRUE)>10)){
    s <- as.numeric(sum(w*as.numeric(x==xs), na.rm=TRUE)/sum(w, na.rm=TRUE))
  }
  else
  {
    s <- NA
  }
}

#' @export
smc <- function(data,
                cols=NULL,
                block=NULL,
                weighDist=TRUE,
                shift=1,
                nCores=1,
                getBg=FALSE,
                longFormat=FALSE, ...){

      if(weighDist){
        pos <- data[["pos"]]
        w <- .weightDecay(pos, ...)}

      if(!is.null(block)){
      scores <- data[,mclapply(.SD,function(x,w){
          xs <- data.table::shift(x, type="lead", n=shift)
          if(weighDist)
          {
            wo <- w
            wo[is.na(x) | is.na(xs)] <- NA
            s <- .smc(x, xs, wo)
          }
          else
          {
            s <- .smc(x, xs)
          }
          s}, w, mc.cores=nCores), by=c(block), .SDcols=cols]}
      else
      {
        scores <- data[,mclapply(.SD, function(x,w){
          xs <- data.table::shift(x, type="lead", n=shift)
          if(weighDist)
          {
            wo <- w
            wo[is.na(x) | is.na(xs)] <- NA
            s <- .smc(x, xs, wo)
          }
          else
          {
            s <- .smc(x, xs)
          }
          s}, w, mc.cores=nCores), .SDcols=cols]
      }

      if(getBg)
      {
        scores <- .getBackground(data, "smc",
                                 scores,
                                 cols=cols,
                                 block=block,
                                 nCores=nCores,
                                 ...)
      }
      else if(longFormat & !getBg)
      {
        # switch to long
        scores <- melt(scores, id.vars=block,
                       value.name="smc")
        setorderv(scores, cols=block)
      }

      return(scores)
    }

shuffleByPos <- function(x, pos){
  covPos <- which(!is.na(x))
  covX <- x[!is.na(x)]
  covX <- sample(covX, length(covX), replace=FALSE)
  sX <- rep(NA, length(x))
  sX[covPos] <- covX

  return(sX)
}

.getBackground <- function(data,
                          score,
                          scores,
                          cols=NULL,
                          block=NULL,
                          nIterations=1000,
                          posConst=TRUE,
                          byBlock=TRUE,
                          nBins=50,
                          seed=42,
                          nCores=1, ...){
  set.seed(seed)
  scoreFun <- get(score)

  if(!byBlock)
  {
    subBlocks <- sample(unique(data[[block]]), nBins)
    data <- subset(data, block %in% subBlocks)
  }

  bg <- mclapply(1:nIterations, function(i, ...){
    #data <- data[sample(nrow(data)),]

    if(posConst)
    {
      # only shuffle on original non-NA positions
      data[, c(cols) := lapply(.SD, shuffleByPos, pos), .SDcols=cols, by=block]
    }
    else{
      data[, c(cols) := .SD[sample(.I, .N)], .SDcols=cols, by=block]}

    bg <- scoreFun(data, cols=cols, block=block, nCores=1,
                   getBg=FALSE, longFormat=FALSE, ...)
    gc()
    bg
  }, ..., mc.cores=nCores)

  bg <- rbindlist(bg)

  # TODO: Add min & max
  # get percentiles of score vs background
  bg$type <- "bg"
  scores$type <- "obs"

  bgm <- rbind(scores, bg, use.names=TRUE)
  qts <- bgm[,(cols):=lapply(.SD, frank, na.last="keep"), .SDcols=cols, by=block]
  qts <- subset(qts, type=="obs")
  qts$type <- NULL
  scores$type <- NULL
  qts <- melt(qts, id.vars=block, value.name="rank")
  qts[,(paste("percentile", score, sep="_")):=rank/(nIterations+1)]

  # switch also observed scores to long
  scores <- melt(scores, id.vars=block,
                 value.name=score)

  # if rows are ordered that could also just be a cbind
  scores <- merge(scores, qts,
                  by.x=c(block, "variable"),
                  by.y=c(block, "variable"),
                  all.x=TRUE)
  return(scores)
}

