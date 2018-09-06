library(ggplot2)
library(reshape2)
library(scales)



args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)

#cell    annotation      correlation

renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}

makePlot <- function(d){
    p <- ggplot(d, aes(x=annotation, y=correlation)) +
        geom_bar(aes(fill=annotation), stat="identity", position="dodge", width=0.8, colour='black') +
        geom_text(data=d[d$corr_pval <= 0.05,], label="*") +
        facet_wrap(~cell, nrow=1) +
        theme(strip.text.x=element_text(size=9),
              axis.text.x = element_blank(),
              axis.text.y=element_text(size=7),
              axis.title.y=element_text(size=7),
              panel.background = element_rect(fill="white", colour="black"),
              panel.grid=element_blank(),
              legend.position="bottom",
              axis.ticks.x=element_blank(),
              strip.background=element_rect(fill="white", colour="black")
              ) +
        scale_fill_brewer(palette="Set1", name="") +
        labs(y="Correlation of eQTL fold enrichment\nwith lclESI bin number", x="")
    return(p)
}


d <- subset(d, ! annotation %in% c("ActiveTss","typicalEnhancerChromatin") )

pdf(args[2], height=2.5, width=7)
makePlot(renameAnnotations(d))
dev.off()

