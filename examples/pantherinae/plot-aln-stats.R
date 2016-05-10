stats.files <- commandArgs(TRUE)

data <- lapply(stats.files, read.table, header=T)

names(data) <- c('data', paste0('sim-', 1:(length(data)-1)))

## scale ietms with percentage values
for (i in seq(1, length(data))) {
    # turn gap frequency in percentage
    data[[i]][,'gap_freq'] <- data[[i]][,'gap_freq'] * 100
    data[[i]][,'prop_invar'] <- data[[i]][,'prop_invar'] * 100
}

vars <- colnames(data[[1]])[- ( which(colnames(data[[1]]) %in% c('file', 'species')))]

pdf('alnstats.pdf', width=12, height=10)
par(mfrow=c(4,3))

# make name mapping for plot titles and axis labels
mapping <- data.frame(row.names=colnames(data[[1]]))
mapping$title <- 0;
mapping$ylab <- 0;
mapping['del_avg_size', 'title'] <- 'Average size of deletions in alignment'
mapping['del_avg_size', 'ylab'] <- 'size (nucleotides)'
mapping['del_count', 'title'] <- 'Number of deletions in alignment'
mapping['del_count', 'ylab'] <- 'nucleotides'
mapping['gap_freq', 'title'] <- '% gaps per sequence'
mapping['gap_freq', 'ylab'] <- '%'
mapping['gaps_per_seq', 'title'] <- 'Number of gaps per sequence'
mapping['gaps_per_seq', 'ylab'] <- 'gap count'
mapping['ins_avg_size', 'title'] <- 'Average size of insertions in alignment'
mapping['ins_avg_size', 'ylab'] <- 'size (nucleotides)'
mapping['ins_count', 'title'] <- 'Number of insertions in alignment'
mapping['ins_count', 'ylab'] <- 'nucleotides'
mapping['nchar', 'title'] <- 'Number of nucleotides per sequence'
mapping['nchar', 'ylab'] <- 'nucleotides'
mapping['ntax', 'title'] <- 'Number of taxa in alignment'
mapping['ntax', 'ylab'] <- 'taxa'
mapping['prop_invar', 'title'] <- '% of invariant sites in alignment'
mapping['prop_invar', 'ylab'] <- '%'



## make boxplots for all valiables
for ( v in vars ) {    
    a <- boxplot( lapply(data, '[[', v), main=mapping[v, 'title'], ylab=mapping[v, 'ylab'], names=names(data), outline=F, col=c('grey', rep('white', length(data)-1)), cex.main=1.5, cex.axis=1.2, cex.lab=1.5)

    ## code for plotting histograms
    ##    colors <- sapply(palette(), scales::alpha, .5)
    ##    for (i in seq(1, length(data))) { 
    ##        hist( data[[i]][,v], prob=T, add=!i==1, col=colors[i], breaks=20)
    ##    }
    ##    colors <- palette()
    ##    for (i in seq(1, length(data))) {        
    ##        lines( density(data[[i]][,v], adjust=2), col=colors[i] )
    ##    }   
}

## calculate and plot number of alignments per species and story in list
spec.stats <- list()
for (i in seq(1, length(data))) {
    specs <- as.numeric(unlist(sapply(as.character(data[[i]]$species), strsplit, split=',')))
    unique.specs <- unique(specs)
    df <- data.frame('species'=unique.specs, 'counts'=sapply(unique.specs, function(x)sum(specs==x)))
    spec.stats[[i]] <- df
}
boxplot(lapply(spec.stats, '[[', 'counts'), main="Number of alignments per species", ylab="alignments", names=names(data), cex.axis=1.5)


dev.off()
