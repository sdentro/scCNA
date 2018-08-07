\name{plotSeg}
\alias{plotSeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
plot track of a single sample
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotSeg(tracksSingle, colSeg = rgb(0.5, 0.2, 0.5, 0.7), lwdSeg = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tracksSingle}{
%%     ~~Describe \code{tracksSingle} here~~
}
  \item{colSeg}{
%%     ~~Describe \code{colSeg} here~~
}
  \item{lwdSeg}{
%%     ~~Describe \code{lwdSeg} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (tracksSingle, colSeg = rgb(0.5, 0.2, 0.5, 0.7), lwdSeg = 2,
    ...)
{
    medianS <- median(unlist(lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed
    })))
    tracksSingle$lCTS <- lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed - medianS
        return(x)
    })
    breaks <- c(0, cumsum(sapply(tracksSingle$lSegs, function(x) max(x$output$loc.end))))/1e+06
    plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
        xlim = c(0, max(breaks)), ylim = c(-2, 4), xlab = "Genomic Position",
        ylab = "relative copy number", frame = F, ...)
    axis(side = 1)
    axis(side = 2)
    meansSeg <- lapply(1:length(tracksMerged$lSegs), function(i) {
        out <- tracksMerged$lSegs[[i]]$output
        means <- lapply(1:nrow(out), function(x) {
            isIn <- tracksSingle$lCTS[[i]]$start > out$loc.start[x] &
                tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            if (sum(isIn) < 2)
                return(list(roundmu = NA, mu = NA, sd = NA))
            mu <- median(tracksSingle$lCTS[[i]]$smoothed[isIn],
                na.rm = T)
            sd <- mad(tracksSingle$lCTS[[i]]$smoothed[isIn],
                na.rm = T)
            list(roundmu = round(2^(mu + log2(2)) - 2), mu = mu,
                sd = sd)
        })
    })
    for (i in 1:length(tracksSingle$lSegs)) {
        segments(tracksSingle$lCTS[[i]]$start/1e+06 + breaks[i],
            tracksSingle$lCTS[[i]]$smoothed, tracksSingle$lCTS[[i]]$end/1e+06 +
                breaks[i], tracksSingle$lCTS[[i]]$smoothed, col = rgb(0.7,
                0.7, 0.7, 0.6))
        segments(tracksMerged$lSegs[[i]]$output$loc.start/1e+06 +
            breaks[i], sapply(meansSeg[[i]], function(x) x$mu),
            tracksMerged$lSegs[[i]]$output$loc.end/1e+06 + breaks[i],
            sapply(meansSeg[[i]], function(x) x$mu), lwd = lwdSeg,
            col = rgb(0.4, 0.4, 0.4, 0.4))
        segments(tracksMerged$lSegs[[i]]$output$loc.start/1e+06 +
            breaks[i], round(sapply(meansSeg[[i]], function(x) x$roundmu)),
            tracksMerged$lSegs[[i]]$output$loc.end/1e+06 + breaks[i],
            round(sapply(meansSeg[[i]], function(x) x$roundmu)),
            lwd = 2.5, col = rgb(1, 0.5, 0.5))
    }
    abline(h = 0, v = breaks, lwd = 1, lty = 2, col = rgb(0.6,
        0.6, 0.6, 0.4))
    text(x = breaks[2:length(breaks)] - 25, y = 4, names(breaks)[2:length(breaks)],
        cex = 0.4)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
