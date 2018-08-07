getRefGenome <- function(fasta=FASTA, CHRS=paste("",c(1:22,"X","Y","MT")))
{
    dna <- readDNAStringSet(fasta,format="fasta")
    dna <- lapply(1:length(CHRS),function(x) dna[[x]])
    names(dna) <- CHRS ##paste(c(1:22,"X","Y","MT"),sep="")
    return(dna)
}
