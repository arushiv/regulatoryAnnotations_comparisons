library(ggplot2)
library(scales)

args <- commandArgs(TRUE)
d <- read.table(snakemake@input[[1]], header=T)
# cell1   feature1        lenght1 cell2   feature2        length2 overlap

renameAnnotations <- function(d, colname){
    d[, colname] <- gsub("hotRegions", "HOT\nRegions", d[, colname])
    d[, colname] <- gsub("broadDomains", "Broad\nDomains", d[, colname])
    d[, colname] <- gsub("stretchEnhancer", "Stretch\nEnhancers", d[, colname])
    d[, colname] <- gsub("superEnhancers", "Super\nEnhancers", d[, colname])
    d[, colname] <- gsub("typicalEnhancers", "Typical\nEnhancers", d[, colname])
    return(d)
}

###### Overlap fraction heatmap
makePlot <- function(d){
    d$feature2 <- factor(d$feature2, levels=c("Typical\nEnhancers","Super\nEnhancers","Stretch\nEnhancers","HOT\nRegions","Broad\nDomains"))
    colvector <- c(0, 0.33*max(d$fraction1), 0.66*max(d$fraction1), max(d$fraction1))
    colvector

    p <- ggplot(d[order(d$feature2),], aes(x = cell1, y = cell2)) +
        geom_tile(aes(fill=fraction1)) +
        facet_grid(feature2~feature1) +
        labs(x="Feature 1", y="Feature 2") +
        scale_fill_gradientn(colours=c("white","yellow", "red", "darkred"), values=rescale(colvector), breaks=c(colvector), labels=c(round(colvector,2)), name="Intersection/\nlength(Feature 1))") +
        theme(strip.text=element_text(size=7), text = element_text(size=10), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=8), axis.text.y = element_text(size=8), panel.background = element_rect(fill="grey", colour="black"), panel.grid=element_blank(), legend.position="bottom") 
    return(p)
}

d$fraction1 <- d$overlap/d$lenght1
d$fraction2 <- d$overlap/d$length2
d$jaccard <- d$overlap/(d$lenght1 + d$length2 - d$overlap)
d <- renameAnnotations(d, "feature1")
d <- renameAnnotations(d, "feature2")

pdf(snakemake@output[[1]], height=4.5, width=4)
makePlot(d)
dev.off()

