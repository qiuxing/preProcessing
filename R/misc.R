## misc. useful functions

## this is a function that plots boxplots together with the actual
## data points
boxplot2 <- function(y, grp, pch=1, ...) {
  boxplot(y~grp, boxwex=.5, outpch=NA, ylim=range(y), ...)
  stripchart(y~grp, vertical=TRUE, method="jitter", pch = pch, add=TRUE)
}
