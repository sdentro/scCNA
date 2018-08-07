computePlotCoverageTrack <-
function(cT,
                                     starts,
                                     ends,
                                     chr,
                                     lwdSeg=1,
                                     colSeg=rgb(0,0,0,.4),
                                     counts=NA,
                                     seg=T,
                                     plot=F)
{
    if(is.na(counts[1]))
        counts <- log10(cT$records+1)
    if(plot)
    {
        plot(0,0,col=rgb(0,0,0,0),frame=F,
             xlim=c(0,max(cT$end))/1000,
             ylim=c(min(counts),max(counts)),main=cT$space[1],
             xaxt="n",yaxt="n",xlab="genomic position (Mb)",ylab="log(counts)")
        axis(side=1,
             at=seq(0,max(cT$end),max(cT$end)/10)/1000,
             signif(seq(0,max(cT$end),max(cT$end)/10)/1000000,2))
        axis(side=2,ylab="log(counts)")
        ##    for(i in 1:nrow(cT))
        segments(cT$start/1000,counts,cT$end/1000,
                 counts,
                 lwd=lwdSeg,
                 col=colSeg)
    }
    if(seg)
    {
        require(DNAcopy)
        segments<- segmentTrack(counts,
                                chr,
                                starts,
                                ends)
        segs <- segments$output
        if(plot)
        {
            for(i in 1:length(segments$output$ID))
            {
                segments(segs$loc.start/1000,
                         segs$seg.mean,
                         segs$loc.end/1000,
                         segs$seg.mean,lwd=2,col=rgb(.6,0,0,.5))
            }
        }
        return(segs)
    }
    return(NULL)
}
