\name{indexBams}
\alias{indexBams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
wrapper to index all bams in a directory using samtools (must be installed)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
indexBams(bams, mc.cores = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bams}{
%%     ~~Describe \code{bams} here~~
}
  \item{mc.cores}{
%%     ~~Describe \code{mc.cores} here~~
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
function (bams, mc.cores = 10)
{
    require(parallel)
    mclapply(bams, function(x) {
        cmd <- paste0("samtools index ", x)
        system(cmd, wait = T)
    }, mc.cores = mc.cores)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
