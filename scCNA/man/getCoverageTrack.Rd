\name{getCoverageTrack}
\alias{getCoverageTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get Coverage Track
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
get coverage track of uniquely-mapped-read counts in predefined genomic
bins on a given chromosome from a specified BAM file given
starts and ends of the bins
}
\usage{
getCoverageTrack(bamPath, chr, starts, ends, CHRSTRING = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamPath}{
    %%     ~~Describe \code{bamPath} here~~
    character string: full path to the BAM file
}
  \item{chr}{
%%     ~~Describe \code{chr} here~~
character string: chromosome to count reads from
}
  \item{starts}{
    %%     ~~Describe \code{starts} here~~
    numeric vector: starts of the bins
}
  \item{ends}{
%%     ~~Describe \code{ends} here~~
numeric vector: ends of the bins (follows order in object starts)
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
  reads in the given bins, and will only count the
  non-duplicated reads. It also requires package GenomicRanges.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns numeric vector of the records, which are the raw read counts
(integer values).
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

}
