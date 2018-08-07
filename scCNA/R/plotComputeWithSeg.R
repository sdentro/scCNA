plotComputeWithSeg <- function(tracksSingle,
                               colSeg=rgb(.5,.2,.5,.7),
                               lwdSeg=2,
                               REFs=c(1:23),
                               rhopsi=c(NA,NA),
                               ...)
{
    medianS <- mean(unlist(lapply(REFs,function(x)
    {
        oo <- tracksSingle$lSegs[[x]]$output
        log2(10)*inverse.rle(list(values=oo$seg.mean,lengths=oo$num.mark))
    })),na.rm=T)
    tracksSingle$lCTS <- lapply(tracksSingle$lCTS,function(x)
    {
        x$smoothed <- log2(10)*x$smoothed-medianS
        return(x)
    })
    tracksSingle$lSeg <- lapply(tracksSingle$lSeg,function(x)
    {
        x$output$seg.mean <- log2(10)*x$output$seg.mean-medianS
        return(x)
    })
    purity <- rhopsi[1]
    ploidy <- rhopsi[2]
    errs <- NULL
    if(is.na(rhopsi[1]))
    {
        meansSeg <- unlist(lapply(1:length(tracksSingle$lSegs),function(i)
        {
            out <- tracksSingle$lSegs[[i]]$output
            means <- unlist(lapply(1:nrow(out),function(x)
            {
                isIn <- tracksSingle$lCTS[[i]]$start>out$loc.start[x] & tracksSingle$lCTS[[i]]$start<=out$loc.end[x]
                if(sum(isIn)<2) return(NA)
                mu <- mean(tracksSingle$lCTS[[i]]$smoothed[isIn],na.rm=T)
                mu
            }))
        }))
        weights <- unlist(lapply(1:length(tracksSingle$lSegs),function(i)
        {
            out <- tracksSingle$lSegs[[i]]$output$num.mark
        }))
        sdsSeg <- unlist(lapply(1:length(tracksSingle$lSegs),function(i)
        {
            out <- tracksSingle$lSegs[[i]]$output
            sds <- unlist(lapply(1:nrow(out),function(x)
            {
                isIn <- tracksSingle$lCTS[[i]]$start>out$loc.start[x] & tracksSingle$lCTS[[i]]$start<=out$loc.end[x]
                if(sum(isIn)<2) return(NA)
                sds <- sd(tracksSingle$lCTS[[i]]$smoothed[isIn],na.rm=T)
                sds
            }))
        }))
        purs <- seq(0.05,1,.01)
        ploidies <- seq(1.7,6,.02)
        errs <- matrix(NA, length(purs),length(ploidies))
        rownames(errs) <- purs
        colnames(errs) <- ploidies
        for(pp in 1:length(purs))
        {
            for(pl in 1:length(ploidies))
            {
                errs[pp,pl] <- geterrors(rho=purs[pp],phi=ploidies[pl],meansSeg,weights,(pl/pp)*0+1)
            }
        }
        mins <- arrayInd(which.min(errs),dim(errs))
        purity <- purs[mins[1]]
        ploidy <- ploidies[mins[2]]
        plot(0,0,col=rgb(0,0,0,0),xlab="",ylab="",xaxt="n",yaxt="n",frame=F,xlim=c(0,1),ylim=c(0,1))
        rasterImage(errs/max(errs),0,0,1,1)
        points(length(ploidies)/mins[2],1-length(purs)/mins[1],col="dark blue",pch=19)
    }
    meansSeg <- lapply(1:length(tracksSingle$lSegs),function(i)
    {
        out <- tracksSingle$lSegs[[i]]$output
        means <- lapply(1:nrow(out),function(x)
        {
            isIn <- tracksSingle$lCTS[[i]]$start>out$loc.start[x] & tracksSingle$lCTS[[i]]$start<=out$loc.end[x]
            if(sum(isIn)<2) return(list(roundmu=NA,mu=NA,sd=NA,start=out$loc.start[x],end=out$loc.end[x]))
            mu <- median(tracksSingle$lCTS[[i]]$smoothed[isIn],na.rm=T)
            sd <- mad(tracksSingle$lCTS[[i]]$smoothed[isIn],na.rm=T)
            list(roundmu=mytransform(mu,purity,ploidy),mu=mu,sd=sd,start=out$loc.start[x],end=out$loc.end[x])
        })
    })
    if(is.na(rhopsi[1]))
    {
        breaks <- c(0,cumsum(sapply(tracksSingle$lSegs,function(x) max(x$output$loc.end))))/1000000
        plot(0,0,col=rgb(0,0,0,0),
             xaxt="n",
             yaxt="n",
             xlim=c(0,max(breaks)),
             ylim=c(0,8),
             xlab="Genomic Position",
             ylab="relative copy number",
             frame=F,...)
        axis(side=1)
        axis(side=2)
        for(i in 1:length(tracksSingle$lSegs))
        {
            segments(tracksSingle$lCTS[[i]]$start/1000000+breaks[i],
                     mytransform(tracksSingle$lCTS[[i]]$smoothed,purity,ploidy),
                     tracksSingle$lCTS[[i]]$end/1000000+breaks[i],
                     mytransform(tracksSingle$lCTS[[i]]$smoothed,purity,ploidy),
                     col=rgb(.7,.7,.7,.6))
            segments(tracksSingle$lSegs[[i]]$output$loc.start/1000000+breaks[i],
                     mytransform(sapply(meansSeg[[i]],function(x) x$mu),purity,ploidy),
                     tracksSingle$lSegs[[i]]$output$loc.end/1000000+breaks[i],
                     mytransform(sapply(meansSeg[[i]],function(x) x$mu),purity,ploidy),
                     lwd=lwdSeg,
                     col=rgb(.4,.4,.4,.4))
            segments(tracksSingle$lSegs[[i]]$output$loc.start/1000000+breaks[i],
                     round(sapply(meansSeg[[i]],function(x) x$roundmu)),
                     tracksSingle$lSegs[[i]]$output$loc.end/1000000+breaks[i],
                     round(sapply(meansSeg[[i]],function(x) x$roundmu)),
                     lwd=2.5,
                     col=rgb(1,.5,.5))
        }
        abline(h=0,v=breaks,lwd=1,lty=2,col=rgb(.6,.6,.6,.4))
        text(x=breaks[2:length(breaks)]-25,y=5.5,names(breaks)[2:length(breaks)],cex=.4)
        mtext(side=3,paste0("guess purity=",purity,"; guess ploidy=",ploidy))
    }
    list(purity=purity,ploidy=ploidy,profile=meansSeg, errs=errs)
}
