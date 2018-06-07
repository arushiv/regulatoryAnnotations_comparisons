library(ggplot2)
library(reshape2)
library(scales)
library(gplots)
library(ggpubr)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)

findmin <- function(d, colname){
    return(min(d[,c(colname)], na.rm=T))
}
findmax <- function(d, colname){
    return(max(d[,c(colname)], na.rm=T))
}

makeHeatmap <- function(d, colname){
    minE <- findmin(d, colname)
    maxE <- findmax(d, colname)
    colscale <- c(minE, 0, maxE*0.33, 0.66*maxE, maxE)
    print(colscale)
    d$annotation2 <- factor(d$annotation2, levels=c("Typical\nEnhancers","Super\nEnhancers","Stretch\nEnhancers","HOT\nRegions","Broad\nDomains"))
    p <- ggplot(d[order(d$annotation2),], aes(x = cell1, y = cell2)) +
        geom_tile(aes(fill=l2fold)) +
        facet_grid(annotation2 ~ annotation1) +
        scale_fill_gradientn(colours=c("blue", "white", "orange", "red", "darkred"), values=rescale(colscale), breaks=colscale, labels=colscale, name="log2(Fold enrichment)") +
        theme(text = element_text(size=7), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), panel.background = element_rect(fill="grey"), panel.grid = element_blank()) ## +
        ## guides(fill=guide_legend(title="log2(Fold enrichment)"))
    return(p)
}

renameAnnotations <- function(d, colname){
    d[, colname] <- gsub("hotRegions", "HOT\nRegions", d[, colname])
    d[, colname] <- gsub("broadDomains", "Broad\nDomains", d[, colname])
    d[, colname] <- gsub("stretchEnhancer", "Stretch\nEnhancers", d[, colname])
    d[, colname] <- gsub("superEnhancers", "Super\nEnhancers", d[, colname])
    d[, colname] <- gsub("typicalEnhancers", "Typical\nEnhancers", d[, colname])
    return(d)
}

d <- renameAnnotations(d, "annotation1")
d <- renameAnnotations(d, "annotation2")

d <- d[! ((d$cell1 == d$cell2) & (d$annotation1 == d$annotation2)), ]

d$minuslog10p <- -log10(d$pval)
bonferroni_threshold <- 0.05/400
d <- d[d$pval <= bonferroni_threshold,]

pdf(args[2], height=4.5, width=5.5)
makeHeatmap(d, "l2fold")
dev.off()


