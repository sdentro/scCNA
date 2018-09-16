getTrackForAll <-
function(bamfile,
                           window,
                           lSe=NULL,
                           lGCT=NULL,
                           normalBamfile=NULL,
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
    
    if (!is.null(normalBamfile)) {
      lCTn <- lapply(allchr, function(chr) getCoverageTrack(bamPath=normalBamfile,
                                                            chr=paste0(chr),
                                                            lSe[[chr]]$starts,
                                                            lSe[[chr]]$ends,
                                                            isDuplicate=isDuplicate, 
                                                            isSecondaryAlignment=isSecondaryAlignment,
                                                            isNotPassingQualityControls=isNotPassingQualityControls,
                                                            isUnmappedQuery=isUnmappedQuery,
                                                            mapqFilter=mapqFilter))
    }
    
    if(is.null(lGCT)) { print("get GC content"); lGCT <- lapply(allchr,function(chr) gcTrack(chr,lSe[[chr]]$starts,lSe[[chr]]$ends)) }
    ## ##################################################
    print("correct for GC content")
    lCTS <- smoothCoverageTrackAll(lCT,lSe,lGCT)
    
    if (!is.null(normalBamfile)) {
      lCTSn <- smoothCoverageTrackAll(lCTn,lSe,lGCT)
    }
    
    gc()
    
    ## ##################################################
    if (!is.null(normalBamfile)) {
      print("correct for matched normal")
      for (i in 1:length(lCTS)) {
        tum = lCTS[[i]]
        norm = lCTSn[[i]]
        
        tum$smoothed_nolog = 10^tum$smoothed
        norm$smoothed_nolog = 10^norm$smoothed
        tum$smoothed_corrected_nolog = tum$smoothed_nolog / norm$smoothed_nolog
        tum$smoothed = log10(tum$smoothed_corrected_nolog)
        lCTS[[i]] = tum
      }
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
