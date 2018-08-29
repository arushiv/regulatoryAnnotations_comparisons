library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(stats)


args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)

renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}



makePlot <- function(d){
    p <- ggplot(d, aes(y = minuslog_pval, x = enrichment)) +
        geom_point(aes(colour=annotation, shape=annotation), stat="identity", size=3) +
        scale_shape_manual(values=c(15,16,17,18,25)) +
        facet_grid(cell~trait, labeller = label_wrap_gen(), scales="free_y") + #, scales="free") +
        theme(axis.text=element_text(size=8),
              axis.title.x=element_text(size=9),
              axis.title.y=element_text(size=9),
              strip.text.x=element_text(size=7),
              panel.background = element_rect(fill="white", colour="black"),
              legend.position="bottom",
              legend.justification="center",
              legend.text=element_text(size=8),
              legend.key.size=unit(4,"mm"),
              panel.grid.major.x=element_line(colour="grey", linetype="dashed", size=0.1),
              panel.grid.major.y=element_line(colour="grey", linetype="dashed", size=0.1),
              strip.background=element_rect(fill="white", colour="black")
              ) +
        labs(y="-log(p value)", x="Fold enrichment") +
        scale_color_brewer(palette="Set1", name="") +
        geom_point(data=subset(d, minuslog_pval < bonferroni_threshold),
                   aes(shape=annotation),
                   colour="grey", size=3) +
        geom_hline(yintercept=bonferroni_threshold, colour="red") +
        geom_vline(xintercept=1, colour="black")
    return(p)
}


d <- d[d$overlap >= as.numeric(args[3]),]
d$enrichment <- d$overlap/d$expected_overlap
d$log2_enrichment <- log2(d$overlap/d$expected_overlap)
d$minuslog_pval <- -(log10(d$pval))

d$trait <- gsub("Crohn_s", "Crohn's", d$trait)
d$trait <- gsub("_", " ", d$trait)

mytraits <- c("Celiac disease","Crohn's disease","Type 1 diabetes","Rheumatoid arthritis","Systemic lupus erythematosus","BMI","Type 2 diabetes")
mytraits1 <- c("Type 2 diabetes","Rheumatoid arthritis")

bonferroni_threshold <- c(-log10(0.05/30))


d <- renameAnnotations(d)

## pdf(args[2], height=4, width=100)
## makePlot(d)
## dev.off()

d1 <- d[((d$trait %in% mytraits) & (d$cell %in% c("GM12878","HepG2"))),]
d1$trait <- factor(d1$trait, levels=mytraits)
pdf(args[2], height=4, width=8)
makePlot(d1)
dev.off()
