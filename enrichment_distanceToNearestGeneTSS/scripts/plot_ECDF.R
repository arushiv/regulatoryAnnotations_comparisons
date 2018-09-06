library(ggplot2)
library(plyr)

args <- commandArgs(TRUE)
d <- read.table(snakemake@input[[1]], header=T)

                                        # ecdf_y  cell    annotation        xvals   ecdf_y_shuffleMean      ecdf_y_shuffleStd

## d <- subset(d, cell=="GM12878")

renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}


makeEcdfPlot <- function(d){
    p <- ggplot(d, aes(x=xvals, y=ecdf_y)) +
        geom_line(aes(colour = annotation), size=0.5) +
        facet_wrap(~cell, nrow=1) +
        labs(y="Fraction", x="log10(Distance(bp) + 1)") +
        theme(strip.text.x = element_text(size = 8),
              panel.background = element_rect(fill = 'white', colour='black'),
              axis.text.x=element_text(size=8),
              axis.text.y=element_text(size=8),
              axis.title=element_text(size=7),
              text=element_text(size=8),
              panel.grid.major.x=element_line(colour="grey", linetype="dashed", size=0.1),
              panel.grid.major.y=element_line(colour="grey", linetype="dashed", size=0.1),
              legend.position="bottom",
              legend.key=element_rect(fill="white"),
              legend.key.size=unit(7,"mm"),
              strip.background=element_rect(fill="white", colour="black")) +
        scale_colour_brewer(palette="Set1", name="") 
    return(p)

}

d <- renameAnnotations(d)
pdf(snakemake@output[[1]], height=2.5, width=7)
makeEcdfPlot(d)
dev.off()

