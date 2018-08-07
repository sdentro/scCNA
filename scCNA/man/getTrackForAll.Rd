\name{getTrackForAll}
\alias{getTrackForAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
get track for all samples
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getTrackForAll(bamfile, window, lSe = NULL, lGCT = NULL, allchr = 1:22, sdNormalise = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamfile}{
%%     ~~Describe \code{bamfile} here~~
}
  \item{window}{
%%     ~~Describe \code{window} here~~
}
  \item{lSe}{
%%     ~~Describe \code{lSe} here~~
}
  \item{lGCT}{
%%     ~~Describe \code{lGCT} here~~
}
  \item{allchr}{
%%     ~~Describe \code{allchr} here~~
}
  \item{sdNormalise}{
%%     ~~Describe \code{sdNormalise} here~~
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
function (bamfile, window, lSe = NULL, lGCT = NULL, allchr = 1:22,
    sdNormalise = 0)
{
    print("get Start-End of segments")
    if (is.null(lSe))
        lSe <- lapply(allchr, function(chr) getStartsEnds(window,
            paste0(chr)))
    print("get Coverage Track")
    lCT <- lapply(allchr, function(chr) getCoverageTrack(bamPath = bamfile,
        chr = paste0(chr), lSe[[chr]]$starts, lSe[[chr]]$ends))
    print("get GC content")
    if (is.null(lGCT))
        lGCT <- lapply(allchr, function(chr) gcTrack(chr, lSe[[chr]]$starts,
            lSe[[chr]]$ends))
    print("correct for GC content")
    lCTS <- smoothCoverageTrackAll(lCT, lSe, lGCT)
    gc()
    print("segment Tracks")
    lSegs <- lapply(1:length(lCTS), function(x) {
        require(DNAcopy)
        segments <- segmentTrack(lCTS[[x]]$smoothed, chr = paste0(x),
            sd = sdNormalise, lSe[[x]]$starts, lSe[[x]]$ends)
    })
    names(lSegs) <- paste0(1:length(lCT))
    tracks <- list(lCTS = lCTS, lSegs = lSegs)
    return(tracks)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
