getCoverageTrack <- function (bamPath, chr, starts, ends, CHRSTRING = "", isDuplicate=F, isSecondaryAlignment=F, isNotPassingQualityControls=T, isUnmappedQuery=T, mapqFilter=0)
{
    require(Rsamtools)
    require(GenomicRanges)
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate=isDuplicate, isSecondaryAlignment=isSecondaryAlignment, isNotPassingQualityControls=isNotPassingQualityControls, isUnmappedQuery=isUnmappedQuery),
    which = GRanges(paste0(CHRSTRING, chr), IRanges(starts,ends)), mapqFilter=mapqFilter)
    coverageTrack <- countBam(bamPath, param = sbp)
    coverageTrack$records <- coverageTrack$records
    return(coverageTrack)
}
