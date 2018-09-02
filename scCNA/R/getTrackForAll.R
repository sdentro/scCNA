getTrackForAll <-
function(bamfile,
                           window,
                           lSe=NULL,
                           lGCT=NULL,
                           mappa=NULL,
                           repli=NULL,
                           allchr=1:22,
                           sdNormalise=0,
                           segmentation_alpha=0.01,
			   isDuplicate=F, 
			   isSecondaryAlignment=F,
			   isNotPassingQualityControls=T,
			   isUnmappedQuery=T, 
			   mapqFilter=0)
{
    
    if(is.null(lSe)) { print("get Start-End of segments"); lSe <- lapply(allchr,function(chr) getStartsEnds(window,paste0(chr))) } 
    ## ##################################################
    print("get Coverage Track")
    lCT <- lapply(allchr, function(chr) getCoverageTrack(bamPath=bamfile,
                                                         chr=paste0(chr),
                                                         lSe[[chr]]$starts,
                                                         lSe[[chr]]$ends,
							 isDuplicate=isDuplicate, 
							 isSecondaryAlignment=isSecondaryAlignment,
							 isNotPassingQualityControls=isNotPassingQualityControls,
							 isUnmappedQuery=isUnmappedQuery,
							 mapqFilter=mapqFilter))
    
    if(is.null(lGCT)) { print("get GC content"); lGCT <- lapply(allchr,function(chr) gcTrack(chr,lSe[[chr]]$starts,lSe[[chr]]$ends)) }
    ## ##################################################
    if (!is.null(mappa)) {
      print("correct for mappability")
      lCTS <- smoothCoverageTrackAll(lCT,lSe,mappa)
    }
    
    print("correct for GC content")
    lCTS <- smoothCoverageTrackAll(lCT,lSe,lGCT)
    
    if (!is.null(repli)) {
      print("correct for replication timing")
      lCTS <- smoothCoverageTrackAll(lCT,lSe,repli)
    }
    
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
                                lSe[[x]]$ends,
                                alpha=segmentation_alpha)
    })
    names(lSegs) <- paste0(1:length(lCT))
    tracks <- list(lCTS=lCTS,lSegs=lSegs)
    return(tracks)
}
