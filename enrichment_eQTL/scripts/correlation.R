library(ggplot2)
library(reshape2)
library(scales)
library(gtools)


args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)

# eqtl    bin     job     cell    annotation      feature overlap expected_overlap        pval

d$bin <- gsub("lclESIbin","",d$bin)

d$enrichment <- d$overlap/d$expected_overlap

d$enrichment[is.na(d$enrichment)] <- 0

write.table(read.table(text = unlist(lapply(unique(d$cell), function(i){
    d1 <- subset(d, cell == i)
    lapply(unique(d1$annotation), function(j){
        d2 <- subset(d1, annotation == j)
        print(d2)
        print(paste(i,j,cor(as.numeric(d2$bin), d2$enrichment), sep='\t'))
        ## return(i,j,cor(as.numeric(d2$bin), d2$log2_enrichment, method="spearman"))
    })
}), use.names=FALSE)), args[2], sep="\t", quote=FALSE, row.names=FALSE, col.names=c("cell","annotation","correlation"))



