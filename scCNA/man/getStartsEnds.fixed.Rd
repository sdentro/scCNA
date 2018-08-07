\name{getStartsEnds.fixed}
\alias{getStartsEnds.fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
get starts and ends of fixed non-overlapping neighbouring bins of a
chromosome given its length
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
get starts and ends of fixed non-overlapping neighbouring bins of a
chromosome given its length
}
\usage{
getStartsEnds.fixed(window, lengthChr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{window}{
    %%     ~~Describe \code{window} here~~
    size of the bins in base pairs
}
  \item{lengthChr}{
    %%     ~~Describe \code{lengthChr} here~~
    length of the chromosome in base pairs
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
function (window, lengthChr)
{
    divideChr <- seq(0, lengthChr, window)
    starts <- divideChr[-c(length(divideChr))] + 1
    ends <- divideChr[-c(1)]
    list(starts = starts, ends = ends)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
