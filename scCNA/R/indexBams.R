indexBams <-
function(bams,mc.cores=10)
{
    require(parallel)
    mclapply(bams,function(x) {
        cmd <- paste0("samtools index ", x)
        system(cmd,wait=T)
    },mc.cores=mc.cores)
}
