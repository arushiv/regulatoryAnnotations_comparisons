library(ggplot2)
library(stringr)
library(ggpubr)
library(gridExtra)
library(powerEQTL)
library(reshape2)

args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T, sep='\t')

# annotation      count   mean    std     min     p_25    p_50    p_75    max

renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}

makePlot <- function(d, colname, n){
    p <- ggplot(d, aes_string(x="variable", y=colname, group="annotation")) +
        geom_point(aes(fill=annotation), stat="identity", shape=21, size=0.4) +
        geom_line(aes(colour=annotation)) +
        labs(y="Power to detect eQTL", x="Observed eQTL effect size in annotations", title=paste("N=", n, ", MAF=", maf, ", typeI=", type1, ", nTests=", nTests, ", error_stddev=1", sep="")) +
        scale_fill_brewer(name="", palette="Set1") +
        scale_colour_brewer(name="", palette="Set1") +
        scale_x_continuous(limits=c(10, 90), breaks=seq(10, 90, by=10)) +
        theme(title=element_text(size=6), axis.title=element_text(size=8), axis.text=element_text(size=7), strip.text.x = element_text(size = 8), panel.background = element_rect(fill = 'white', colour='black'), panel.grid=element_blank(), legend.position="bottom", legend.key.size=unit(3, "mm"), legend.text=element_text(size=6)) +
        guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow=2))
    return(p)
}

d <- d[, !(names(d) %in% c('min','max','std','mean','X50.0.1'))]
d1 <- melt(d, id=c('annotation','count'))
maf = 0.2
type1 = 0.05
nTests = 100000

d1$power_100 = powerEQTL.SLR(maf, typeI = type1, nTests = nTests, slope=d1$value, myntotal = 100, verbose=FALSE, mystddev=1)
## d1$power_100 = powerEQTL.SLR(0.5, typeI = 0.05, nTests = 200000, slope=d1$value, myntotal = 100, verbose=FALSE, mystddev=1)
d1$power_250 = powerEQTL.SLR(maf, typeI = type1, nTests = nTests, slope=d1$value, myntotal = 250, verbose=FALSE, mystddev=1)
d1$power_500 = powerEQTL.SLR(maf, typeI = type1, nTests = nTests, slope=d1$value, myntotal = 500, verbose=FALSE, mystddev=1)

d2 <- renameAnnotations(d1)
d2$variable <- as.numeric(gsub("X", "", d2$variable))

pdf(args[2], height=2.5, width=3)
makePlot(d2, "power_100", 100)
makePlot(d2, "power_250", 250)
makePlot(d2, "power_500", 500)
dev.off()

