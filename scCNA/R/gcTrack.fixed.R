gcTrack.fixed <-
function(chr,
                          pas)
{
    gc <- rowSums(letterFrequencyInSlidingView(dna[[chr]],
                                               pas,
                                               c("G","C")))/pas
    lGC <- length(gc)
    gc <- gc[seq(pas,lGC,pas)-round(pas/2)]
    names(gc) <- paste("bin",seq(pas,lGC,pas)-pas+1,sep="-")
    gc
}
