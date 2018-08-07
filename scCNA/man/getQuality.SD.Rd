\name{getQuality.SD}
\alias{getQuality.SD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
getQuality of sample from its profile
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getQuality.SD(alltracks, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alltracks}{
%%     ~~Describe \code{alltracks} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (alltracks, plot = F)
{
    quality <- colSums(sapply(alltracks, function(x) sapply(x,
        function(y) mytry(sum(y$num.mark * y$CNsd)/sum(y$num.mark)))),
        na.rm = T)
    quality <- quality/max(quality)
    if (plot) {
        par(mar = c(8, 5, 1, 1))
        barplot(quality, las = 2)
    }
    quality
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
