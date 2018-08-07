getCoverageTrack <- function (bamPath, chr, starts, ends, CHRSTRING = "")
{
    require(Rsamtools)
    require(GenomicRanges)
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE,isSecondaryAlignment=FALSE),
    which = GRanges(paste0(CHRSTRING, chr), IRanges(starts,
    ends)))
    coverageTrack <- countBam(bamPath, param = sbp)
    coverageTrack$records <- coverageTrack$records
    return(coverageTrack)
}