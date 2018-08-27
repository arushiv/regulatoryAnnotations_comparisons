library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(viridis)
require(reshape)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)
# cell1   feature1        lenght1 cell2   feature2        length2 overlap

renameAnnotations <- function(d, colname){
    d[, colname] <- gsub("hotRegions", "HOT\nRegions", d[, colname])
    d[, colname] <- gsub("broadDomains", "Broad\nDomains", d[, colname])
    d[, colname] <- gsub("stretchEnhancer", "Stretch\nEnhancers", d[, colname])
    d[, colname] <- gsub("superEnhancers", "Super\nEnhancers", d[, colname])
    d[, colname] <- gsub("typicalEnhancers", "Typical\nEnhancers", d[, colname])
    return(d)
}

renameAnnotations_space <- function(d, colname){
    d[, colname] <- gsub("hotRegions", "HOT Regions", d[, colname])
    d[, colname] <- gsub("broadDomains", "Broad Domains", d[, colname])
    d[, colname] <- gsub("stretchEnhancer", "Stretch Enhancers", d[, colname])
    d[, colname] <- gsub("superEnhancers", "Super Enhancers", d[, colname])
    d[, colname] <- gsub("typicalEnhancers", "Typical Enhancers", d[, colname])
    return(d)
}

makePlot_across <- function(d){
    maxE <- max(d$value, na.rm=TRUE )
    colvector <- c(0, 0.33*maxE, 0.66*maxE, maxE)
    colvector

    p <- ggplot(d, aes(x = cell1, y = variable)) +
        geom_tile(aes(fill=value)) +
        facet_wrap(~feature, nrow=1) +
        labs(x="Feature 1", y="Feature 2") +
        scale_fill_viridis(name="Jaccard") +
        theme(strip.text=element_text(size=7),
              text = element_text(size=10),
              axis.text.x = element_text(angle=60, hjust=1, size=8),
              axis.text.y = element_text(size=8),
              panel.background = element_rect(fill="grey", colour="black"),
              panel.grid=element_blank(),
              legend.position="bottom") 
    return(p)
}


makePlot_same <- function(d){
    maxE <- max(d$value, na.rm=TRUE )
    colvector <- c(0, 0.33*maxE, 0.66*maxE, maxE)
    colvector

    p <- ggplot(d, aes(x = feature1, y = variable)) +
        geom_tile(aes(fill=value)) +
        facet_wrap(~cell, nrow=1) +
        labs(x="Feature 1", y="Feature 2") +
        scale_fill_viridis(name="Jaccard") +
        theme(strip.text=element_text(size=7),
              text = element_text(size=10),
              axis.text.x = element_text(angle=60, hjust=1, size=8),
              axis.text.y = element_text(size=8),
              panel.background = element_rect(fill="grey", colour="black"),
              panel.grid=element_blank(),
              legend.position="bottom") 
    return(p)
}


d$fraction1 <- d$overlap/d$lenght1
d$fraction2 <- d$overlap/d$length2
d$jaccard <- d$overlap/(d$lenght1 + d$length2 - d$overlap)


# Fix matrix with jaccard for overlaps across annotations within cell types
samecell <- d[d$cell1 == d$cell2,]

get_lower_tri <- function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

get_same_cell_lower_tri <- function(d){
    lower_tri <- get_lower_tri(dcast(d, feature1 ~ feature2, value.var="jaccard"))
    return(lower_tri)
}


ndf = do.call("rbind", by(samecell, samecell$cell1, get_same_cell_lower_tri))
ndf <- data.frame(add_rownames(ndf, "cell"))
ndf$cell <- gsub(".", "@", ndf$cell, fixed=TRUE)
ndf$cell <- gsub("@.*", "", ndf$cell)
plotdf <- melt(ndf)
plotdf <- renameAnnotations_space(plotdf, "feature1")
plotdf <- renameAnnotations_space(plotdf, "variable")

get_across_cell_lower_tri <- function(d){
    lower_tri <- get_lower_tri(dcast(d, cell1 ~ cell2, value.var="jaccard"))
    return(lower_tri)
}

# Fix matrix with jaccard for overlaps across cell types within annotations
diffcell <- d[d$feature1 == d$feature2,]
diffcell <- diffcell[,c('feature1', 'cell1', 'cell2', 'jaccard')]
ddf = do.call("rbind", by(diffcell, diffcell$feature1, get_across_cell_lower_tri))
ddf <- data.frame(add_rownames(ddf, "feature"))

ddf$feature <- gsub(".", "@", ddf$feature, fixed=TRUE)
ddf$feature <- gsub("@.*", "", ddf$feature)
plotdfd <- melt(ddf)
plotdfd <- renameAnnotations(plotdfd, "feature")



# plot
pdf(args[2], height=3, width=6)
makePlot_same(plotdf)
makePlot_across(plotdfd)
dev.off()

