getChrLengths <-
function(dna)
{
    sapply(names(dna),function(x) length(dna[[x]]))
}
