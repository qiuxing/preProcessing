## trimmed STD; adapted from the same function in chemometrics,
## ver. 1.4.2. Note that this function compute the trimmed STD for all
## columns in X.
## Notable changes made by Xing: (a) the trim proportion can be set to
## any value in [0, 0.8); (b) gives warning if too few observations
## are available after trimming; (c) NA safe.
sd_trim <- function(x,trim=0.2, const=TRUE){
  if (trim < 0 | trim > .8) stop("trim must be a value in (0, .8)")
  x <- as.matrix(x)
  if (const){
    ## cc is the attenuation factors of trimmed STD associated with
    ## trim values seq(.1, .8, by=.1), obtained by Monte Carlo method
    ## (see trimSD1.r).
    trim.grid <- seq(0, .8, by=0.1)
    cc <- c(1, 0.7892, 0.6615, 0.5564, 0.4633, 0.3777, 0.2972, 0.2202, 0.1457)
    const <- approx(trim.grid, cc, xout=trim, rule=2)$y
  } else {
    const <- 1
  }
  m <- apply(x,2, function(x) mean(x, trim=trim, na.rm=TRUE))
  res <- x-rep(1,nrow(x))%*%t(m)
  qu <- apply(abs(res),2, function(x) quantile(x,1-trim, na.rm=TRUE))
  ## compute the trimmed SS and effective Ns
  ns <- ss <- rep(0, ncol(x))
  for (j in 1:ncol(x)){
    res.j <- res[,j]
    remain.j <- which(abs(res.j)<=qu[j])
    ns[j] <- length(remain.j)
    ss[j] <- sum(res.j[remain.j]^2)
  }
  if (min(ns) <3 ) warning("Less than 3 observations are available after trimming. Consider using a smaller trim value. ")
  sdtrim <- sqrt(ss/(ns-1))/const
  return(sdtrim)
}

## Lot specific mean/STD normalization. to="overall" means to
## normalize all groups to the overall mean/median; otherwise, to must
## be set to a specific level in grp.
normalize <- function(x, grp, method=c("meanSTD", "medianIQR"), to="overall", trim=0, max.b1=3) {
  method <- match.arg(method)
  if (method=="meanSTD"){
    locfun <- function(x) mean(x, trim=trim, na.rm=TRUE)
    scalefun <- function(x) sd_trim(x, trim=trim)
  } else if (method=="medianIQR") {
    locfun <- function(x) median(x, na.rm=TRUE)
    scalefun <- function(x) IQR(x, na.rm=TRUE)
  } else {
    stop("We only implemented meanSTD and medianIQR normalization at this moment.")
  }
  ## drop nonexistinig levels in grp
  if (is.factor(grp)) grp <- droplevels(grp)
  xloc <- tapply(x, grp, locfun); xscale <- tapply(x, grp, scalefun)
  if (to=="overall"){
    toloc <- locfun(x); toscale <- scalefun(x)
  } else {
    if (to %in% names(xloc)) {
      toloc <- xloc[to]; toscale <- xscale[to]
    } else {
      stop(paste('Please ensure that "to" is one of the valid levels in grp:', paste(names(xloc), collapse=", ")))
    }
  }
  b1 <- pmin(toscale/xscale, max.b1)
  b0 <- toloc -b1*xloc
  ## 06/23/22 fix
  beta0 <- b0[grp]; beta1 <- b1[grp]
  betas <- cbind(beta0=beta0, beta1=beta1)
  return(list(bs=cbind(b0=b0, b1=b1), betas=betas, xnormed=beta0 +x*beta1))
}

## Xvar: the outcome variable (e.g., "Conc.mean.log2"). analytevar:
## name of the analyte in dataset (e.g., "Analyte"). batchvar: name of
## the batch variable in dataset (e.g., "Lot.Panel.Plate").
## 12/22/22: New feature: it can be used without "analytevar"
varpropfun <- function(dataset, Xvar, batchvar, analytevar=NULL) {
  if (is.null(analytevar)) {
    ## by default (analytevar==NULL), all observations are considered
    ## as one analyte.
    X <- dataset[[Xvar]]; Batch <- dataset[[batchvar]]
    ss <- summary(aov(X~factor(Batch)))[[1]]
    varprops <- ss[1,"Sum Sq"]/(var(X, na.rm=TRUE)*(sum(!is.na(X))-1))
  } else {
    analytes <- dataset[[analytevar]]
    VARS <- unique(analytes)
    varprops <- sapply(VARS, function(v) {
      ids <- analytes==v
      X <- dataset[[Xvar]][ids]; Batch <- dataset[[batchvar]][ids]
      ss <- summary(aov(X~factor(Batch)))[[1]]
      return(ss[1,"Sum Sq"]/(var(X, na.rm=TRUE)*(sum(!is.na(X))-1)))
    }); names(varprops) <- VARS
  }
  return(varprops)
}

## identify outlying batches and only normalize those batches
BatchNormalize <- function(x, grp, method=c("meanSTD", "medianIQR"), nMAD=1.5, trim=0, max.b1=3, ...){
  method <- match.arg(method)
  if (method=="meanSTD"){
    locfun <- function(x) mean(x, trim=trim, na.rm=TRUE)
    scalefun <- function(x) sd_trim(x, trim=trim)
  } else if (method=="medianIQR") {
    locfun <- function(x) median(x, na.rm=TRUE)
    scalefun <- function(x) IQR(x, na.rm=TRUE)
  } else {
    stop("We only implemented meanSTD and medianIQR normalization at this moment.")
  }
  ## drop nonexistinig levels in grp
  if (is.factor(grp)) grp <- droplevels(grp)
  xloc <- tapply(x, grp, locfun); xscale <- tapply(x, grp, scalefun)
  toloc <- locfun(x); toscale <- scalefun(x)
  out.grps <- union(names(which(is.na(GrpHampel(xloc, nMAD=nMAD, ...)))),
                    names(which(is.na(GrpHampel(xscale, nMAD=nMAD, ...)))))
  for (g in out.grps) {
    ids <- which(grp==g)
    b1 <- pmin(toscale/xscale[g], max.b1)
    b0 <- toloc -b1*xloc[g]
    x[ids] <- b0 +x[ids]*b1
  }
  return(list(xnormed=x, out.grps=out.grps))
}
