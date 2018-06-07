library(ggplot2)
library(plyr)

args <- commandArgs(TRUE)
d <- read.table(snakemake@input[[1]], header=T)

                                        # ecdf_y  cell    annotation        xvals   ecdf_y_shuffleMean      ecdf_y_shuffleStd

sqn <- (100)^0.5
d <- subset(d, cell=="GM12878")
d$enrichment <- d$ecdf_y/d$ecdf_y_shuffleMean

d$ymax <- d$ecdf_y/(d$ecdf_y_shuffleMean - 1.96 * d$ecdf_y_shuffleSem)
d$ymin <- d$ecdf_y/(d$ecdf_y_shuffleMean + 1.96 * d$ecdf_y_shuffleSem)

d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)

dmain <- d[ , which(names(d) %in% c("cell","annotation","ecdf_y","xvals","bin"))]
dshuffle <- d[ , which(names(d) %in% c("cell","annotation","xvals","ecdf_y_shuffleMean","ecdf_y_shuffleSem"))]
dshuffle <- rename(dshuffle, c("ecdf_y_shuffleMean"="ecdf_y"))

makeEnrichmentPlot <- function(d){
    p <- ggplot(d, aes(x=xvals, y=enrichment, colour=as.factor(bin))) +
        geom_errorbar(aes(ymax=ymax, ymin=ymin), alpha=0.4) +
        geom_line(size=0.4) +
        facet_wrap(~annotation, nrow=1) +
        labs(y="TSS proximity enrichment", x="log10(Distance to gene TSS (bp) + 1)") +
        theme(strip.text.x = element_text(size = 8), panel.background = element_rect(fill = 'white', colour='black'), axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), axis.title=element_text(size=7), text=element_text(size=8), panel.grid=element_blank(), legend.position="bottom", legend.key.size=unit(7,"mm")) +
        geom_hline(yintercept=1, size=0.3, colour="black") +
        scale_colour_manual(values=c("purple","blue","lightblue","orange","red"), name="lclESI bin for genes")##  +
        ## scale_colour_brewer(palette="Oranges", name="lclESI bin for genes")
    return(p)
}

makeEcdfPlot <- function(dmain, dshuffle){
    dmain <- subset(dmain, cell %in% c("GM12878","H1"))
    dshuffle <- subset(dshuffle, cell %in% c("GM12878", "H1"))
    p <- ggplot(dshuffle, aes(x=xvals, y=ecdf_y)) +
        geom_line(data=dmain, aes(colour = as.factor(bin)), size=0.5) +
        geom_line(data=dshuffle, colour="black", size=0.4) +
        geom_errorbar(aes(ymax=(ecdf_y + ecdf_y_shuffleSem), ymin=(ecdf_y - ecdf_y_shuffleSem)), colour="grey", alpha=0.6) +
        facet_grid(cell~annotation) +
        labs(y="Fraction", x="log10(Distance to gene TSS (bp) + 1)") +
        theme(strip.text.x = element_text(size = 8), panel.background = element_rect(fill = 'white', colour='black'), axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), axis.title=element_text(size=7), text=element_text(size=8), panel.grid=element_blank(), legend.position="bottom", legend.key.size=unit(7,"mm")) +
        scale_colour_manual(values=c("purple","blue","lightblue","orange","red"), name="lclESI bin for genes")##  +
        ## scale_colour_brewer(palette="Spectral", name="lclESI bin for genes") 
    return(p)

}

pdf(snakemake@output[[1]], height=2.5, width=7)
makeEnrichmentPlot(d)
makeEcdfPlot(dmain, dshuffle)
dev.off()

