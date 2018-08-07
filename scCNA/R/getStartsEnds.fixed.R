getStartsEnds.fixed <-
function(window,
        lengthChr)
{
    divideChr <- seq(0, lengthChr, window)
    starts <- divideChr[-c(length(divideChr))] + 1
    ends <- divideChr[-c(1)]
    list(starts=starts,ends=ends)
}
