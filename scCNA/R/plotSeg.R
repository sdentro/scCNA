plotSeg <- function (tracksSingle, colSeg = rgb(0.5, 0.2, 0.5, 0.7), lwdSeg = 2,
...)
{
    tracksMerged <- tracksSingle
    medianS <- median(unlist(lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed
    })))
    tracksSingle$lCTS <- lapply(tracksSingle$lCTS, function(x) {
        x$smoothed <- log2(10) * x$smoothed - medianS
        return(x)
    })
    breaks <- c(0, cumsum(sapply(tracksSingle$lSegs, function(x) max(x$output$loc.end))))/1e+06
    plot(0, 0, col = rgb(0, 0, 0, 0), xaxt = "n", yaxt = "n",
    xlim = c(0, max(breaks)), ylim = c(-2, 4), xlab = "Genomic Position",
    ylab = "relative copy number", frame = F, ...)
    axis(side = 1)
    axis(side = 2)
    meansSeg <- lapply(1:length(tracksSingle$lSegs), function(i) {
        out <- tracksMerged$lSegs[[i]]$output
        means <- lapply(1:nrow(out), function(x) {
            isIn <- tracksSingle$lCTS[[i]]$start > out$loc.start[x] &
            tracksSingle$lCTS[[i]]$start <= out$loc.end[x]
            if (sum(isIn) < 2)
            return(list(roundmu = NA, mu = NA, sd = NA))
            mu <- median(tracksSingle$lCTS[[i]]$smoothed[isIn],
            na.rm = T)
            sd <- mad(tracksSingle$lCTS[[i]]$smoothed[isIn],
            na.rm = T)
            list(roundmu = round(2^(mu + log2(2)) - 2), mu = mu,
            sd = sd)
        })
    })
    for (i in 1:length(tracksSingle$lSegs)) {
        segments(tracksSingle$lCTS[[i]]$start/1e+06 + breaks[i],
        tracksSingle$lCTS[[i]]$smoothed, tracksSingle$lCTS[[i]]$end/1e+06 +
        breaks[i], tracksSingle$lCTS[[i]]$smoothed, col = rgb(0.7,
        0.7, 0.7, 0.6))
        segments(tracksMerged$lSegs[[i]]$output$loc.start/1e+06 +
        breaks[i], sapply(meansSeg[[i]], function(x) x$mu),
        tracksMerged$lSegs[[i]]$output$loc.end/1e+06 + breaks[i],
        sapply(meansSeg[[i]], function(x) x$mu), lwd = lwdSeg,
        col = rgb(0.4, 0.4, 0.4, 0.4))
        segments(tracksMerged$lSegs[[i]]$output$loc.start/1e+06 +
        breaks[i], round(sapply(meansSeg[[i]], function(x) x$roundmu)),
        tracksMerged$lSegs[[i]]$output$loc.end/1e+06 + breaks[i],
        round(sapply(meansSeg[[i]], function(x) x$roundmu)),
        lwd = 2.5, col = rgb(1, 0.5, 0.5))
    }
    abline(h = 0, v = breaks, lwd = 1, lty = 2, col = rgb(0.6,
    0.6, 0.6, 0.4))
    text(x = breaks[2:length(breaks)] - 25, y = 4, names(breaks)[2:length(breaks)],
    cex = 0.4)
}
