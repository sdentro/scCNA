myheat <-
function(mat,
                   keep1=rep(T,ncol(mat)),
                   scaleY=.3,
                   fundist="manhattan",
                   funclust="ward.D2")
{
    require(gplots)
    layout(mat=cbind(c(3,4),c(2,1),c(2,1)),
           widths=c(.4,1,1),
           heights=c(5,1.4))
    par(mar=c(4,4,1,1))
    xx <- plotAllGenome(mat,scaleY=scaleY)
    mmm <- t(mat)
    hcc <- hclust(dist(mmm,method=fundist),met=funclust)
    ord <- hcc$order
    mmm[mmm>1] <- 2
    mmm[mmm< -1] <- -2
    mmm <- mmm[ord[length(ord):1],]
    im <- array(0.8,dim=c(nrow(mmm),ncol(mmm),3))
    im[,,1] <- getCols(mmm,"R")
    im[,,2] <- getCols(mmm,"G")
    im[,,3] <- getCols(mmm,"B")
    par(mar=c(0,4,0,1))
    plot(0,0,col=rgb(0,0,0,0),xlim=c(0,ncol(mmm)),ylim=c(0,1),frame=F,xaxt="n",yaxt="n",xlab="",ylab="")
    rasterImage(im,0,0,ncol(mmm),1)
    abline(v=which(diff(as.numeric(gsub("(.*):(.*)","\\1",
                                        colnames(mmm))))!=0)+1,
           col=rgb(1,1,1,.8))
    axis(side=2,at=(1:nrow(mmm)-.5)/nrow(mmm),rownames(mmm),
         cex.axis=.6,col=rgb(.5,.5,.5,.5),las=2,tick=F,hadj=.4)
    par(mar=c(1.5,0,1.5,0))
    gplots:::plot.dendrogram(as.dendrogram(hcc),
                             horiz = TRUE,
                             axes = FALSE,
                             yaxs = "i",
                             leaflab = "none")
    plot(0,0,col=rgb(0,0,0,0),xlim=c(0,1),ylim=c(0,1),
         frame=F,
         xaxt="n",
         yaxt="n",
         xlab="",
         ylab="")
    legend("center",col=c(rgb(getCols(-2,"R"),getCols(-2,"G"),
                              getCols(-2,"B")),
                          rgb(getCols(-1,"R"),getCols(-1,"G"),
                              getCols(-1,"B")),
                          rgb(getCols(0,"R"),getCols(0,"G"),getCols(0,"B")),
                          rgb(getCols(1,"R"),getCols(1,"G"),getCols(1,"B")),
                          rgb(getCols(2,"R"),getCols(2,"G"),getCols(2,"B"))),
           legend=c("HD","LOH","N","Gain","Amp"),
           pch=19,cex=1.5,box.col=rgb(0,0,0,0))
}
