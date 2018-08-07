\name{plotComputeWithSeg}
\alias{plotComputeWithSeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
plot and compute total copy number states given a segmented track.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
plot and compute total copy number states given a segmented track.
}
\usage{
plotComputeWithSeg(tracksSingle, colSeg = rgb(0.5, 0.2, 0.5, 0.7),
lwdSeg = 2, REFs = c(1:23), rhopsi=c(NA, NA),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tracksSingle}{
    %%     ~~Describe \code{tracksSingle} here~~
    segmented track
}
  \item{colSeg}{
    %%     ~~Describe \code{colSeg} here~~
    colour
}
  \item{lwdSeg}{
%%     ~~Describe \code{lwdSeg} here~~
width
}
  \item{REFs}{
%%     ~~Describe \code{REFs} here~~
chromosomes taken as reference to derive purity and ploidy (default = all)
}
\item{rhopsi}{
%%     ~~Describe \code{REFs} here~~
vector of two double values: purity and ploidy
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
passed to plotting function
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
    REFs = c(1:23), ...)
{
    medianS <- mean(unlist(lapply(REFs, function(x) {
        oo <- tracksSingle$lSegs[[x]]$output
        log2(10) * inverse.rle(list(values = oo$seg.mean, lengths = oo$num.mark))
    })), na.rm = T)
    tracksSingle$lCTS <- lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed - medianS
        return(x)
    })
    tracksSingle$lSeg <- lapply(tracksSingle$lSeg, function(x) {
        x$output$seg.mean <- log2(10) * x$output$seg.mean - medianS
        return(x)
    })
    tracksMerged <- tracksSingle
    meansSeg <- unlist(lapply(1:length(tracksMerged$lSegs), function(i) {
        out <- tracksMerged$lSegs[[i]]$output
        means <- unlist(lapply(1:nrow(out), function(x) {
            isIn <- tracksSingle$lCTS[[i]]$start > out$loc.start[x] &
                tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            if (sum(isIn) < 2)
                return(NA)
            mu <- mean(tracksSingle$lCTS[[i]]$smoothed[isIn],
                na.rm = T)
            mu
        }))
    }))
    weights <- unlist(lapply(1:length(tracksMerged$lSegs), function(i) {
        out <- tracksMerged$lSegs[[i]]$output$num.mark
    }))
    purs <- seq(0.05, 1, 0.01)
    ploidies <- seq(1.7, 6, 0.02)
    errs <- matrix(NA, length(purs), length(ploidies))
    rownames(errs) <- purs
    colnames(errs) <- ploidies
    for (pp in 1:length(purs)) {
        for (pl in 1:length(ploidies)) {
            errs[pp, pl] <- geterrors(rho = purs[pp], phi = ploidies[pl],
                meansSeg, weights)
        }
    }
    mins <- arrayInd(which.min(errs), dim(errs))
    purity <- purs[mins[1]]
    ploidy <- ploidies[mins[2]]
    if (F) {
        plot(0, 0, col = rgb(0, 0, 0, 0), xlab = "", ylab = "",
            xaxt = "n", yaxt = "n", frame = F, xlim = c(0, 1),
            ylim = c(0, 1))
        rasterImage(errs/max(errs), 0, 0, 1, 1)
        points(length(ploidies)/mins[2], 1 - length(purs)/mins[1],
            col = "dark blue", pch = 19)
    }
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
            list(roundmu = mytransform(mu, purity, ploidy), mu = mu,
                sd = sd)
        })
    })
    breaks <- c(0, cumsum(sapply(tracksSingle$lSegs, function(x) max(x$output$loc.end))))/1e+06
    plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
        xlim = c(0, max(breaks)), ylim = c(0, 8), xlab = "Genomic Position",
        ylab = "relative copy number", frame = F, ...)
    axis(side = 1)
    axis(side = 2)
    for (i in 1:length(tracksMerged$lSegs)) {
        segments(tracksSingle$lCTS[[i]]$start/1e+06 + breaks[i],
            mytransform(tracksSingle$lCTS[[i]]$smoothed, purity,
                ploidy), tracksSingle$lCTS[[i]]$end/1e+06 + breaks[i],
            mytransform(tracksSingle$lCTS[[i]]$smoothed, purity,
                ploidy), col = rgb(0.7, 0.7, 0.7, 0.6))
        segments(tracksMerged$lSegs[[i]]$output$loc.start/1e+06 +
            breaks[i], mytransform(sapply(meansSeg[[i]], function(x) x$mu),
            purity, ploidy), tracksMerged$lSegs[[i]]$output$loc.end/1e+06 +
            breaks[i], mytransform(sapply(meansSeg[[i]], function(x) x$mu),
            purity, ploidy), lwd = lwdSeg, col = rgb(0.4, 0.4,
            0.4, 0.4))
        segments(tracksMerged$lSegs[[i]]$output$loc.start/1e+06 +
            breaks[i], round(sapply(meansSeg[[i]], function(x) x$roundmu)),
            tracksMerged$lSegs[[i]]$output$loc.end/1e+06 + breaks[i],
            round(sapply(meansSeg[[i]], function(x) x$roundmu)),
            lwd = 2.5, col = rgb(1, 0.5, 0.5))
    }
    abline(h = 0, v = breaks, lwd = 1, lty = 2, col = rgb(0.6,
        0.6, 0.6, 0.4))
    text(x = breaks[2:length(breaks)] - 25, y = 5.5, names(breaks)[2:length(breaks)],
        cex = 0.4)
    mtext(side = 3, paste0("guess purity=", purity, "; guess ploidy=",
        ploidy))
    list(purity, ploidy)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
