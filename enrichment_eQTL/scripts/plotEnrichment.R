library(ggplot2)
library(stringr)

args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T, sep='\t')

                                        # eqtl bin job cell annotation feature overlap expected_overlap pval

renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}

makePlot <- function(d){
    p <- ggplot(d, aes(x=as.factor(bin), y=Enrich)) +
        geom_point(aes(fill=annotation), size=2, shape=21, colour="black") +
        geom_line(aes(group=annotation, colour=annotation)) +
        labs(y="eQTL Fold enrichment", x="LCL ESI bin") +
        scale_fill_brewer(name="", palette="Set1") +
        scale_colour_brewer(name="", palette="Set1") +
        facet_wrap(~cell, nrow=1) +
        theme(axis.text.x=element_text(size=9), strip.text.x = element_text(size = 9), panel.background = element_rect(fill = 'white', colour='black'), panel.grid=element_blank(), legend.position="bottom")
    return(p)
}

d <- subset(d, ! annotation %in% c("ActiveTss","typicalEnhancerChromatin") )
d$Enrich <- (d$overlap/d$expected_overlap)
d$bin <- gsub("0_lclESIbin", "", d$bin)

pdf(args[2], height=3, width=7)
makePlot(renameAnnotations(d))
dev.off()

