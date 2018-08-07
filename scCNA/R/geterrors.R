geterrors <- function(rho,phi,meansSeg,weights, sds)
{
    signal <- mytransform(meansSeg,rho,phi)
    mean(((round(signal)-signal)/sds)^2*weights/1000,na.rm=T)
}
