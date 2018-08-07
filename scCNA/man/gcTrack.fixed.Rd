\name{gcTrack.fixed}
\alias{gcTrack.fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
derive gc content of the bins (fixed size)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gcTrack.fixed(chr, pas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{pas}{
%%     ~~Describe \code{pas} here~~
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
function (chr, pas)
{
    gc <- rowSums(letterFrequencyInSlidingView(dna[[chr]], pas,
        c("G", "C")))/pas
    lGC <- length(gc)
    gc <- gc[seq(pas, lGC, pas) - round(pas/2)]
    names(gc) <- paste("bin", seq(pas, lGC, pas) - pas + 1, sep = "-")
    gc
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
