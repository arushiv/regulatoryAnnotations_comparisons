library(ggplot2)
library(stringr)
library(ggpubr)
library(gridExtra)
library(dplyr)

args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T, sep='\t')
# gene_id gene_name       gene_chr        gene_start      gene_end        strand  num_var beta_shape1     beta_shape2     true_df pval_true_df    variant_id      tss_distance    chr     pos     ref     alt     num_alt_per_site        rs_id_dbSNP147_GRCh37p13        minor_allele_samples    minor_allele_count      maf     ref_factor      pval_nominal    slope   slope_se        pval_perm       pval_beta       qval    pval_nominal_threshold


anovaTest <- function(d){
    stat <- lm(formula = absoluteEffect ~ annotation, data=d)
    return(anova(stat))

}

wilcoxTest <- function(d){
    stat <- pairwise.wilcox.test(d$absoluteEffect, d$annotation, paired=FALSE, p.adjust="none", na.rm=TRUE)
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

makePlot <- function(d){
    my_comparisons <- list(c("HOT Regions", "Stretch Enhancers"), c("Broad Domains", "Stretch Enhancers"))
    ptop <- ggplot(d, aes(annotation, absoluteEffect)) +
        theme(axis.text=element_blank(), panel.background = element_rect(fill = 'white'), panel.grid=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
        stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)

    pbottom <- ggplot(d, aes(annotation, absoluteEffect)) +
        geom_boxplot(aes(fill=annotation), width=0.5, outlier.shape=NA, alpha=0.7) +
        scale_y_continuous(limits = quantile(d$absoluteEffect, c(0, .99))) +
        ## geom_jitter(size=0.8, alpha=0.3, shape=16, colour="blue") +
        labs(y="Effect size (eQTL)", x="") +
        scale_fill_brewer(name="", palette="Set1", guide=FALSE) +
        theme(axis.text.x=element_text(size=9, angle=30, hjust=1), strip.text.x = element_text(size = 8), panel.background = element_rect(fill = 'white', colour='black'), panel.grid=element_blank()) 

    p <- grid.arrange(ptop, pbottom, nrow=2, ncol=1, heights=c(2,8))
    ## gA <- ggplotGrob(ptop)
    ## gB <- ggplotGrob(pbottom)
    ## p <- grid::grid.draw(rbind(gA, gB))

    ## resize_heights <- function(g, heights = rep(1, length(idpanels))){
    ##     idpanels <- unique(g$layout[grepl("panel",g$layout$name), "t"])
    ##     g$heights <- grid:::unit.list(g$heights)
    ##     g$heights[idpanels] <- lapply(heights, unit, "null")
    ##     g
    ## }

    ## ## grid.newpage()
    ## p <- grid::grid.draw(resize_heights(rbind(gA, gB), c(3,1)))
    return(p)
}

makePlotWithPvalInside <- function(d){
    my_comparisons <- list(c("HOT Regions", "Stretch Enhancers"), c("Broad Domains", "Stretch Enhancers"))
    ptop <- ggplot(d, aes(annotation, absoluteEffect)) +
        theme(axis.text=element_blank(), panel.background = element_rect(fill = 'white'), panel.grid=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
        stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)

    pbottom <- ggplot(d, aes(annotation, absoluteEffect)) +
        geom_boxplot(aes(fill=annotation), width=0.5, outlier.shape=NA, alpha=0.7) +
        scale_y_continuous(limits = quantile(d$absoluteEffect, c(0, .99))) +
        ## geom_jitter(size=0.8, alpha=0.3, shape=16, colour="blue") +
        labs(y="Effect size (eQTL)", x="") +
        scale_fill_brewer(name="", palette="Set1", guide=FALSE) +
        theme(axis.text.x=element_text(size=9, angle=30, hjust=1), strip.text.x = element_text(size = 8), panel.background = element_rect(fill = 'white', colour='black'), panel.grid=element_blank()) 

    ## p <- grid.arrange(ptop, pbottom, nrow=2, ncol=1, heights=c(2,8))
    return(p)
}

makePlotSimple <- function(d){
    my_comparisons <- list(c("HOT Regions", "Stretch Enhancers"), c("Broad Domains", "Stretch Enhancers"))
    d$annotation <- paste(d$annotation, "\n(N = ", d$count, ")", sep="")    
    p <- ggplot(d, aes(annotation, absoluteEffect)) +
        geom_boxplot(aes(fill=annotation), width=0.5, outlier.shape=NA, alpha=0.7) +
        scale_y_continuous(limits = quantile(d$absoluteEffect, c(0, .99))) +
        ## geom_jitter(size=0.8, alpha=0.3, shape=16, colour="blue") +
        labs(y="Effect size (eQTL)", x="") +
        scale_fill_brewer(name="", palette="Set1", guide=FALSE) +
        stat_compare_means(method = "wilcox.test", comparisons = my_comparisons) +
        theme(axis.text.x=element_text(size=7, angle=30, hjust=1), strip.text.x = element_text(size = 8), panel.background = element_rect(fill = 'white', colour='black'), panel.grid=element_blank()) 

    return(p)
}

distTssAnnot <- function(d){
    d$absDist <- abs(d$tss_distance)
    colvector <- c(rgb(255,0,0,maxColorValue=255), rgb(255,195,77,maxColorValue=255), rgb(194,255,5,maxColorValue=255))
    dtest <- data.frame(compare_means(absDist ~ annotation,  data = d))
    my_comparisons <- list(c("HOT Regions", "Stretch Enhancers"), c("Broad Domains", "Stretch Enhancers"))
    print(dtest)
    p <- ggplot(d, aes(x = annotation, y = absDist)) +
        geom_jitter(aes(colour=annotation), size=0.8) +
        geom_boxplot(aes(fill=annotation), colour="black", width=0.2, outlier.shape=NA, alpha=0.6) +
        labs(y="absolute eQTL dist to TSS", x="ATAC-seq peaks in chromatin states") +
        ## scale_colour_manual(values=colvector, guide=FALSE) +
        scale_y_log10() +
        scale_fill_brewer(name="", palette="Set1", guide=FALSE) +
        scale_color_brewer(name="", palette="Set1", guide=FALSE) +
        theme(axis.text.x=element_text(size=9, angle=30, hjust=1), panel.background = element_rect(fill = 'white', colour='black'), panel.grid=element_blank(), axis.text.y=element_text(size=9)) +
        stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)
    return(p)
}
## d2 <- d[!duplicated(d[c("chrom","indexPos","snp","annotation")]),]   # count SNP (index or proxy) only once for each annotation
d2 <- d[!duplicated(d[c("chrom","SNP","annotation")]),]   # count SNP (index or proxy) only once for each annotation

d2$absoluteEffect <- abs(d2$slope)
d2 <- data.frame(d2 %>% group_by(annotation) %>% mutate(count = n()))

anovaTest(d2)
wilcoxTest(d2)
d2 <- renameAnnotations(d2)


pdf(args[2], height=2.5, width=3)
## makePlot(renameAnnotations(d2))
makePlotSimple(d2)
distTssAnnot(d2)
dev.off()

