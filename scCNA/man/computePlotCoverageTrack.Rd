\name{computePlotCoverageTrack}
\alias{computePlotCoverageTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
derive copy number states and plot results for one sample
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computePlotCoverageTrack(cT, starts, ends, chr, lwdSeg = 1, colSeg = rgb(0, 0, 0, 0.4), counts = NA, seg = T, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cT}{
%%     ~~Describe \code{cT} here~~
}
  \item{starts}{
%%     ~~Describe \code{starts} here~~
}
  \item{ends}{
%%     ~~Describe \code{ends} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{lwdSeg}{
%%     ~~Describe \code{lwdSeg} here~~
}
  \item{colSeg}{
%%     ~~Describe \code{colSeg} here~~
}
  \item{counts}{
%%     ~~Describe \code{counts} here~~
}
  \item{seg}{
%%     ~~Describe \code{seg} here~~
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
function (cT, starts, ends, chr, lwdSeg = 1, colSeg = rgb(0,
    0, 0, 0.4), counts = NA, seg = T, plot = F)
{
    if (is.na(counts[1]))
        counts <- log10(cT$records + 1)
    if (plot) {
        plot(0, 0, col = rgb(0, 0, 0, 0), frame = F, xlim = c(0,
            max(cT$end))/1000, ylim = c(min(counts), max(counts)),
            main = cT$space[1], xaxt = "n", yaxt = "n", xlab = "genomic position (Mb)",
            ylab = "log(counts)")
        axis(side = 1, at = seq(0, max(cT$end), max(cT$end)/10)/1000,
            signif(seq(0, max(cT$end), max(cT$end)/10)/1e+06,
                2))
        axis(side = 2, ylab = "log(counts)")
        segments(cT$start/1000, counts, cT$end/1000, counts,
            lwd = lwdSeg, col = colSeg)
    }
    if (seg) {
        require(DNAcopy)
        segments <- segmentTrack(counts, chr, starts, ends)
        segs <- segments$output
        if (plot) {
            for (i in 1:length(segments$output$ID)) {
                segments(segs$loc.start/1000, segs$seg.mean,
                  segs$loc.end/1000, segs$seg.mean, lwd = 2,
                  col = rgb(0.6, 0, 0, 0.5))
            }
        }
        return(segs)
    }
    return(NULL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
