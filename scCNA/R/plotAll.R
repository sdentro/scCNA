plotAll <- function (alltracks, window, QualThresh=1, Nchromosomes = 22, scaleY=0.9)
{
    nms <- names(alltracks)
    keep <- which(sapply(alltracks, length) == Nchromosomes)
    alltracks <- lapply(keep, function(x) alltracks[[x]])
    names(alltracks) <- nms[keep]
    quality <- getQuality.SD(alltracks,plot=F)
    mat <- lapply(alltracks, function(x) unlist(lapply(x, function(y) {
        v <- rep(y$CNtot, y$num.mark)
        v
    })))
    nms <- unlist(lapply(alltracks[[1]], function(y) {
        starts <- unlist(lapply(1:nrow(y), function(z) {
            pas <- (-y$loc.start[z] + y$loc.end[z])/y$num.mark[z]
            s <- seq(y$loc.start[z], y$loc.end[z] - pas, length.out = y$num.mark[z])
        }))
        ends <- unlist(lapply(1:nrow(y), function(z) {
            pas <- (-y$loc.start[z] + y$loc.end[z])/y$num.mark[z]
            s <- seq(y$loc.start[z] + pas, y$loc.end[z], length.out = y$num.mark[z])
        }))
        paste0(y$chrom[1], ":", starts, "-", ends)
    }))
    isNull <- sapply(mat, is.null)
    mat <- sapply(which(!isNull), function(x) mat[[x]])
    rownames(mat) <- nms
    chrom <- as.numeric(gsub("(.*):(.*)-(.*)", "\\1", rownames(mat)))
    starts <- as.numeric(gsub("(.*):(.*)-(.*)", "\\2", rownames(mat)))
    ends <- as.numeric(gsub("(.*):(.*)-(.*)", "\\3", rownames(mat)))
    myheat(mat[, names(quality[quality <= QualThresh])], scaleY = scaleY)
    NULL
}
