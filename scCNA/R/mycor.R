mycor <-
function(x,...)
{
    x[,colSums(x)==0] <- rnorm(nrow(x)*sum(colSums(x)==0),sd=.01)
    cc <- cor(x,...)
    cc[is.na(cc)] <- 0
    cc
}
