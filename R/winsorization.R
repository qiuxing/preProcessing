## apply observation-level Hampel filter to subgroups of data. When
## Grp is not given, it takes the entire X as one group. When
## winsorize=TRUE, it replaces the outliers by upper/lower bounds;
## otherwise, with NA. If with.out.ids=TRUE, it returns a list of two
## objects, X is the vector of results after outlier
## removal/winsorization, out.ids are the indices of outliers.
GrpHampel <- function(X, Grp=NULL, nMAD=3, type=c("two.sided", "less", "greater"), winsorize=FALSE, with.out.ids=FALSE) {
  type <- match.arg(type)
  if (is.null(Grp)) Grp <- rep(0, length(X))
  out.ids <- c()
  for (gp in unique(Grp)){
    ids <- which(Grp==gp); X.gp <- X[ids]
    Xmed <- median(X.gp, na.rm=TRUE); Xmad <- mad(X.gp, na.rm=TRUE)
    U <- Xmed+nMAD*Xmad; L <- Xmed-nMAD*Xmad
    out.lower.ids <- which(X.gp<L); out.upper.ids <- which(X.gp>U)
    if (type=="two.sided"){
      if (winsorize) {
        X.gp[out.lower.ids] <- L; X.gp[out.upper.ids] <- U
      } else {
        X.gp[out.lower.ids] <- NA; X.gp[out.upper.ids] <- NA
      }
      out.ids <- c(out.ids, ids[out.lower.ids], ids[out.upper.ids])
    } else if (type=="less") {
      if (winsorize) {
        X.gp[out.lower.ids] <- L
      } else {
        X.gp[out.lower.ids] <- NA
      }
      out.ids <- c(out.ids, ids[out.lower.ids])
    } else if (type=="greater") {
      if (winsorize) {
        X.gp[out.upper.ids] <- U
      } else {
        X.gp[out.upper.ids] <- NA
      }
      out.ids <- c(out.ids, ids[out.upper.ids])
    } else {
      stop("Valid choice of type are: two.sided, less, and greater.")
    }
    X[ids] <- X.gp
  }
  if (with.out.ids) {
    return(list(X=X, out.ids=sort(out.ids)))
  } else {
    return(X)
  }
}
