\name{geterrors}
\alias{geterrors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
get errors of a fit obtained from a purity and ploidy combination
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
get errors of a fit obtained from a purity and ploidy combination
}
\usage{
geterrors(rho, phi, meansSeg, weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{
%%     ~~Describe \code{rho} here~~
double value: purity
}
  \item{phi}{
%%     ~~Describe \code{phi} here~~
double value: ploidy
}
  \item{meansSeg}{
%%     ~~Describe \code{meansSeg} here~~
segments observed values for logr
}
\item{weights}{
%%     ~~Describe \code{weights} here~~
segments weights
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rho, phi, meansSeg, weights)
{
    signal <- mytransform(meansSeg, rho, phi)
    mean(((round(signal) - signal)^2 * weights/1000), na.rm = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
