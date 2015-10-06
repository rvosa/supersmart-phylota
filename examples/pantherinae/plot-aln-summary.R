## This script is called using 'Rscript' with argument being a file
## listing alignment locations. It produces a plot showing
## the following properties of the set of alignments:
## 1) number of sequences in the alignments
## 2) sequence lengths
## 3) per cent of gaps per sequence
## 4) indel sizes
## Note that alignments containing < 2 sequences are filtered out.
## Output is a png file, named as the input file with '.png' appended.

library('ape')

aln.list <- commandArgs(TRUE)[1]

files <- scan(aln.list, 'character')
alignments <- lapply(files, read.FASTA)

## filter out alignments with <3 seqs
alignments <- Filter(Negate(is.null), lapply(alignments, function(x)if (length(x)>2) x))

png(paste0(aln.list, '.png'))
par(mfrow=c(2,2))

## 1) 
## plot size distribution of alignments (number of seqs in alignment)
sizedist <- sapply(alignments, length)
numseqs <- sum(sizedist)
highest <- max(sizedist)
hist(sizedist, main=paste0("Alignment size ditribution, \n(" , numseqs, "seqs in ", length(alignments), " alignments)\nmean al. size: ",round(mean(sizedist),2)), xlab="sequences in alignment", xaxt='n', breaks=seq(1, highest+1)+0.5)
axis(side=1, at=seq(1, highest), labels=seq(1, highest))

## 2)
## plot distribution of alignment lengths
lengthdist <- sapply(alignments, function(x)length(x[[1]]))
hist(lengthdist, breaks=10, main=paste0("Sequence lengths\nmean: ", round(mean(lengthdist),2) ), xlab="length of sequence in alignment")


## 3)
## plot distribution of percent gaps per sequence
gap.percent.dist <- unlist(sapply(alignments, function(a){sapply(as.character(a), function(seq) sum(seq=='-')*100/length(seq))}))
hist(gap.percent.dist, breaks=100, main=paste0("% of gaps per sequence\nmean: ", round(mean(gap.percent.dist),2), " %"), xlab="% gaps")

## 4)
## plot distribution of indels
## count indels
indel.sizes <- unlist(lapply( alignments, function(a)  {    
    unname(unlist(sapply(as.character(a), function(seq) {
        gap <- seq=="-"        
        ## Do not count trailing and leading gaps, replace by 'X'
        trail.lead.gap <- which(cumsum(!gap) == 0 | rev(cumsum(rev(!gap))) == 0)
        seq[trail.lead.gap] <- 'X'
        rl <- rle(seq=='-')
        rl$lengths[rl$values]
    })))
}))
hist(indel.sizes, breaks=100, main=paste0("Sizes of indels (", length(indel.sizes), ")\nmean: ", round(mean(indel.sizes)),2), xlab="indel size")

dev.off()
