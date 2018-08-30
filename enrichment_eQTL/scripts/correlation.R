library(ggplot2)
library(reshape2)
library(scales)


args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)

# eqtl    bin     job     cell    annotation      feature overlap expected_overlap        pval

d$bin <- gsub("0_lclESIbin","",d$bin)

d$enrichment <- d$overlap/d$expected_overlap

d$enrichment[is.na(d$enrichment)] <- 0

write.table(read.table(text = unlist(lapply(unique(d$cell), function(i){
    d1 <- subset(d, cell == i)
    lapply(unique(d1$annotation), function(j){
        d2 <- subset(d1, annotation == j)
        print(d2)
        t = cor.test(as.numeric(d2$bin), d2$enrichment)
        print(paste(i, j, t$estimate, t$p.value, sep='\t'))
    })
}), use.names=FALSE)), args[2], sep="\t", quote=FALSE, row.names=FALSE, col.names=c("cell", "annotation", "correlation", "corr_pval"))



