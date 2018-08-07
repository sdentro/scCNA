getQuality.SD <-
function(alltracks,plot=F)
{
    quality <- colSums(sapply(alltracks,function(x) sapply(x,function(y)
        mytry(sum(y$num.mark*y$CNsd)/sum(y$num.mark)))),na.rm=T)
    quality <- quality/max(quality)
    if(plot)
    {
        par(mar=c(8,5,1,1))
        barplot(quality,las=2)
    }
    quality
}
