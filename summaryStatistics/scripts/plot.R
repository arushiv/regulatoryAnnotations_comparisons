library(ggplot2)
library(gridExtra)

args <- commandArgs(TRUE)
d1 <- read.table(args[1], header=T) # Number of segments dataframe
d2 <- read.table(args[2], header=T) # Length distribution dataframe
d3 <- read.table(args[3], header=T) # % genome coverage
                                        #d1: # cell    annotation      number
                                        #d2: # cell    annotation      length
                                        #d3: # length_genome   cell    annotation      length_annotation       

renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}

makePlot <- function(d1, d2, d3){
                                        # Plot number of segments
    p1 <- ggplot(d1, aes(x=annotation, y=number, fill=annotation)) +
        geom_bar(position="dodge", stat="identity", width=0.6, colour="black", size=0.3) +
        facet_wrap(~cell, nrow=1) +
        labs(x="Cell type", y="# of segments") +
        theme(strip.text.x = element_text(size = 6),
              panel.background = element_rect(fill = 'white', colour='black'),
              axis.text.x=element_blank(),
              axis.text=element_text(size=6),
              panel.grid=element_blank(),
              legend.key.size=unit(4, "mm"),
              axis.title=element_text(size=7),
              axis.title.x=element_blank(),
              axis.ticks.x=element_blank()) +
        scale_fill_brewer(name="", palette="Set1", guide=FALSE) 
    
    p2 <- ggplot(d2, aes(annotation,length, fill=annotation)) +
        geom_violin(size=0.5) +
        geom_boxplot(width=0.1, size=0.1, outlier.size=0.01, outlier.shape=NA) +
        facet_wrap(~cell, nrow=1) +
        labs(x="Annotation", y="Length (bp)") +
        theme(strip.text.x = element_blank(),
              panel.background = element_rect(fill = 'white', colour='black'),
              axis.text.x=element_blank(),
              axis.text=element_text(size=6),
              panel.grid=element_blank(),
              legend.key.size=unit(4, "mm"),
              axis.title=element_text(size=7),
              axis.title.x=element_blank(),
              axis.ticks.x=element_blank()) +
        scale_fill_brewer(name="", palette="Set1", guide=FALSE) +
        scale_y_log10(breaks=c(100,1000,10000,100000))
    
    p3 <- ggplot(d3, aes(x=annotation, y=percent_coverage, fill=annotation)) +
        geom_bar(position="dodge", stat="identity", width=0.6, colour="black", size=0.3) +
        facet_wrap(~cell, nrow=1) +
        labs(x="Annotation", y="%genome\ncoverage") +
        theme(strip.text.x = element_blank(),
              panel.background = element_rect(fill = 'white', colour='black'),
              axis.text.x=element_text(size=5.5, angle=60, hjust=1),
              axis.text.y=element_text(size=4.5),
              panel.grid=element_blank(),
              legend.key.size=unit(4, "mm"),
              axis.title=element_text(size=6)) +
        scale_fill_brewer(name="", palette="Set1", guide=FALSE) +
        ylim(0,max(d3$percent_coverage)+0.05)
    
    
    p <- grid.arrange(p1,p2,p3, nrow=3, ncol=1, heights=c(.3,.25,.45))
    return(p)
}

anovaTest <- function(d){
    stat <- lm(formula = length ~ annotation, data=d)
    return(anova(stat))

}

wilcoxTest <- function(d){
    stat <- pairwise.wilcox.test(d$length, d$annotation, paired=FALSE, p.adjust="bonferroni", na.rm=TRUE)
    return(stat)
}

by(d2, d2$cell, anovaTest)
by(d2, d2$cell, wilcoxTest)

d3$percent_coverage <- (d3$length_annotation/d3$length_genome)*100
d1 <- renameAnnotations(d1)
d2 <- renameAnnotations(d2)
d3 <- renameAnnotations(d3)

pdf(args[4], height=3, width=7)
makePlot(d1, d2, d3)
dev.off()


