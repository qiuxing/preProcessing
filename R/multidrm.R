######################################################################
## Useful functions for SC analysis
######################################################################

## the 5-parameter log-logistic function
LL5 <- function(x, params) {
  b <- params[1]; c1 <- params[2]; d <- params[3]; e <- params[4]; f <- params[5]
  return(pmax(0, c1 +(d-c1)/(1 +exp(b*(log(x) -log(e))))^f))
}

## its inverse
LL5inv <- function(y, params) {
  b <- params[1]; c1 <- params[2]; d <- params[3]; e <- params[4]; f <- params[5]
  xhat <- e*( ((d-c1)/(y-c1))^(1/f) -1)^(1/b)
  ## now replace the NAs resulted by y <= c1 by zero
  xhat[y<=c1] <- 0
  return(xhat)
}

## to fit a family of SCs by batches, with trimming based on a 3-STD rule
multidrm <- function(formula, batchid, dataset, nMAD=3, ngrid=101, logdist=FALSE, ...) {
  batch <- dataset[,batchid]; BB <- unique(batch)
  ## Individual SCs
  SCs <- lapply(BB, function(bb) drm(formula, fct=LL.5(), data=dataset[batch==bb,]))
  names(SCs) <- BB
  ## outlier detection
  yx <- all.vars(terms(formula))
  ## xgrid is always in log-scale
  L <- log10(max(1,min(dataset[,yx[2]]))); U <- log10(min(2e5, max(dataset[,yx[2]])))
  xgrid <- 10^(seq(L, U, length.out=ngrid)); dx <- (U-L)/(ngrid-1)
  FImat <- sapply(SCs, function(sc) LL5(xgrid, params=coef(sc)))
  if (logdist) FImat <- log10(FImat)
  colnames(FImat) <- BB
  ## matplot(xgrid, FImat, type="l", col=1, lty=1, log="x")
  DistVec <- sqrt(colSums(sweep(FImat, 1, rowMeans(FImat))^2)*dx)
  ## out.ids <- which(DistVec-mean(DistVec) > nSTD*sd(DistVec))
  ## Using the Hampel filter to identify outliers 
  out.ids <- which(DistVec-median(DistVec) > nMAD*mad(DistVec))
  ## hist(DistVec, breaks=21); rug(DistVec); abline(v=median(DistVec), lwd=2); abline(v=median(DistVec)+nMAD*mad(DistVec), lwd=2, lty=2)
  dataset2 <- dataset[!(batch %in% BB[out.ids]),]
  ## the mean curve
  MeanSC <- drm(formula, fct=LL.5(), data=dataset2)
  return(list(SCs=SCs, MeanSC=MeanSC, out.ids=out.ids, DistVec=DistVec, xmax=10^U, ymax=max(FImat)))
}

