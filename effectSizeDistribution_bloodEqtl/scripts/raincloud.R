library(ggplot2)
library(ggpubr)
library(gridExtra)

library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

args <- commandArgs(TRUE)

d <- read.table(gzfile(args[1]), header=T, sep='\t')


anovaTest <- function(d){
    stat <- lm(formula = absoluteEffect ~ annotation, data=d)
    return(anova(stat))

}

wilcoxTest <- function(d){
    stat <- pairwise.wilcox.test(d$absoluteEffect, d$annotation, paired=FALSE, p.adjust="bonf", na.rm=TRUE)
    return(stat)
}

renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}

makeRaincloud <- function(d){
    raincloud_theme = theme(
        text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position = "bottom",
        plot.title = element_text(lineheight=.8, face="bold", size = 16),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

    d$annotation <- paste(d$annotation, "\n(N = ", d$count, ")", sep="")   
    g <- ggplot(data = d, aes(y = absoluteEffect, x = annotation, fill = annotation)) +
        geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
        geom_point(aes(y = absoluteEffect, color = annotation), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
        geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
        expand_limits(x = 5.25) +
        guides(fill = FALSE) +
        guides(color = FALSE) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
                                        # coord_flip() +
        labs(y="Effect size (eQTL)", x="") +   
        theme_bw() +
        raincloud_theme
    return(g)

}

d2 <- d[!duplicated(d[c("chrom","SNP","annotation")]),]   # count SNP (index or proxy) only once for each annotation

d2$absoluteEffect <- abs(d2$slope)
d2 <- data.frame(d2 %>% group_by(annotation) %>% mutate(count = n()))

anovaTest(d2)
wilcoxTest(d2)
d2 <- renameAnnotations(d2)


pdf(args[2], height=2.5, width=3)
makeRaincloud(d2)
dev.off()



