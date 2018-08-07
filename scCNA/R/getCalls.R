getCalls <- function (tracksSingle, mc.cores = 10)
{
    require(parallel)
    medianS <- median(unlist(lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed
    })))
    tracksSingle$lCTS <- lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed - medianS
        return(x)
    })
    meansSeg <- lapply(1:length(tracksSingle$lSegs), function(i) {
        out <- tracksSingle$lSegs[[i]]$output
        means <- mclapply(1:nrow(out), function(x) {
            cond1 <- tracksSingle$lCTS[[i]]$start >= out$loc.start[x]
            cond2 <- tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            isIn <- cond1 & cond2
            if (sum(isIn) < 2)
            return(list(roundmu = NA, mu = NA, sd = NA))
            mu <- median(tracksSingle$lCTS[[i]]$smoothed[isIn],
            na.rm = T)
            sd <- mad(tracksSingle$lCTS[[i]]$smoothed[isIn],
            na.rm = T)
            list(roundmu = round(2^(mu + log2(2)) - 2), mu = mu,
            sd = sd)
        }, mc.cores = mc.cores)
        t <- cbind(out, sapply(means, function(y) y$roundmu),
        sapply(means, function(y) y$mu), sapply(means, function(y) y$sd))
        colnames(t) <- c(colnames(out), "CNtot", "CNseg", "CNsd")
        t
    })
    return(meansSeg)
}
