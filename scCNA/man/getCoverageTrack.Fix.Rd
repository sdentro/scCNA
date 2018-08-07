\name{getCoverageTrack.Fix}
\alias{getCoverageTrack.Fix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
get Coverage Track (fixed bins)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
get coverage track of uniquely-mapped-read counts in genomic
bins of a given length on a given chromosome and its total length from a specified BAM file
}
\usage{
getCoverageTrack.Fix(bamPath, chr, lengthChr, pas, CHRSTRING="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamPath}{
%%     ~~Describe \code{bamPath} here~~
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{lengthChr}{
%%     ~~Describe \code{lengthChr} here~~
}
  \item{pas}{
%%     ~~Describe \code{pas} here~~
}
  \item{CHRSTRING}{
    %%     ~~Describe \code{pas} here~~
    character string: the character string defining chromosomes in the BAM.
    Typically either empty (i.e. 1,2,3,...,X,Y) or "chr" (i.e. chr1,
    chr2,...,chrX,chrY)
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
   This function relies on package Rsamtools to count the uniquely mapped
   reads in the derived bins, and will only count the
   non-duplicated reads. It also requires package GenomicRanges.
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
  maxime.tarabichi@crick.ac.uk
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
}
