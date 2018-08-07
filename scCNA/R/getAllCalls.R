getAllCalls <-
function(lTracks,
                        MC.CORES=10)
{
    alltracks <- lapply(lTracks,function(tracksSingle)
    {
        res <- mytry(getCalls(tracksSingle,
                              window=window,
                              mc.cores=MC.CORES))
        res
    })
}
