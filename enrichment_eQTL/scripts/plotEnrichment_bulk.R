library(ggplot2)
library(stringr)

args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T, sep='\t')


# chrom   proxyStart      proxyEnd        indexPos        snp     GENE_ID TSSpos  rvalue  P       log10pvalue     Alleles cell    annotation


renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}

makePlot <- function(d){
    p <- ggplot(d, aes(x=annotation, y=enrichment)) +
        geom_bar(aes(fill=annotation), width=0.8, stat="identity", colour="black") +
        labs(y="LCL eQTL fold enrichment", x="") +
        scale_fill_brewer(name="", palette="Set1", guide=FALSE) +
        facet_wrap(~cell, nrow=1) +
        theme(axis.text.x=element_text(size=9, angle=30, hjust=1), strip.text.x = element_text(size = 8), panel.background = element_rect(fill = 'white', colour='black'), panel.grid=element_blank())
    return(p)
}


d <- d[d$annotation %in% c("hotRegions","broadDomains","stretchEnhancer","superEnhancers","typicalEnhancers"),]
d$enrichment <- d$overlap/d$expected_overlap
d <- renameAnnotations(d)
pdf(args[2], height=3.5, width=7)
makePlot(d)
dev.off()

