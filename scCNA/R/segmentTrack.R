segmentTrack <-
function(covtrack,
                         chr,
                         starts,
                         ends=NA,
                         sd=0,
                         min.width=5,
                         alpha = 0.01)
{
    covtrack <- covtrack*rnorm(length(covtrack),mean=1,sd=sd)
    cna <- CNA(covtrack,chr=rep(chr,length(covtrack)),
               maploc=starts,data.type="logratio")
    cna.smoothed <- smooth.CNA(cna)
    segment(cna.smoothed,min.width=min.width,alpha=alpha)
}
