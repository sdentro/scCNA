getCols.nonrounded <-
function(mmm,channel="R")
{
    poss <- as.character(c(-2,-1,0,1,2))
    Rs <- c(.6,.3,0,0,0)
    Gs <- c(0,0,0,.3,.6)
    Bs <- c(.2,.2,0,.5,.5)
    names(Rs) <- names(Gs) <- names(Bs) <- poss
    rmmm <- round(mmm)
    val <- get(paste0(channel,"s"))[as.character(rmmm)]
    pc <- (mmm-rmmm)/.5
    pc[is.nan(pc)] <- 0
    pc[is.infinite(pc)] <- 0
    val. <- val+as.vector(pc)*val
    val.[val.>1] <- 1
    val.[val.<0] <- 0
    val
}
