getTrackForAll <-
function(bamfile,
                           window,
                           lSe=NULL,
                           lGCT=NULL,
                           allchr=1:22,
                           sdNormalise=0)
{
    print("get Start-End of segments")
    if(is.null(lSe)) lSe <- lapply(allchr,function(chr) getStartsEnds(window,paste0(chr)))
    ## ##################################################
    print("get Coverage Track")
    lCT <- lapply(allchr, function(chr) getCoverageTrack(bamPath=bamfile,
                                                         chr=paste0(chr),
                                                         lSe[[chr]]$starts,
                                                         lSe[[chr]]$ends))
    print("get GC content")
    if(is.null(lGCT)) lGCT <- lapply(allchr,function(chr) gcTrack(chr,lSe[[chr]]$starts,lSe[[chr]]$ends))
    ## ##################################################
    print("correct for GC content")
    lCTS <- smoothCoverageTrackAll(lCT,lSe,lGCT)
    gc()
    ## ##################################################
    print("segment Tracks")
    lSegs <- lapply(1:length(lCTS),function(x)
    {
        require(DNAcopy)
        segments<- segmentTrack(lCTS[[x]]$smoothed,
                                chr=paste0(x),
                                sd=sdNormalise,
                                lSe[[x]]$starts,
                                lSe[[x]]$ends)
    })
    names(lSegs) <- paste0(1:length(lCT))
    tracks <- list(lCTS=lCTS,lSegs=lSegs)
    return(tracks)
}
