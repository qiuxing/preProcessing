\name{GrpHampel}
\alias{GrpHampel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generating block design randomizations used for clinical trials.}
\description{
  Generate block design randomizations used for clinical trials.
}
\usage{
GrpHampel(X, Grp=NULL, nMAD=3, type=c("two.sided", "less", "greater"),
          winsorize=FALSE, with.out.ids=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ A vector of input. }
  \item{Grp}{ A vector of group membership. By default (Grp=NULL), the entire X is considered as one group.  }
  \item{nMAD}{ Number of median absolute deviation (MAD) used in the
    Hampel filter.}
  \item{type}{ two-sided, less, or greater. }
  \item{winsorize}{ If winsorize=TRUE, it replace outliers with the
    threshold(s) instead of NA. }
  \item{with.out.ids}{ If with.out.ids=TRUE, return a list of both the
    outlier removed/winsorized data and the indices of outliers. }
}
\details{
  TBD.
}
\value{
  \item{X}{ Outlier removed/winsorized data. }
  \item{out.ids}{ Indices of outliers (useful when winsorize=TRUE).}
}
\references{Some text book on Hampel filter.}
\author{Xing Qiu}
\examples{

## TBD. Write some examples to show the usage of this funtion.

}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{manip}% __ONLY ONE__ keyword per line
