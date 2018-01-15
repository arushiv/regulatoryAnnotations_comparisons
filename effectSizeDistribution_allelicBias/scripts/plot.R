library(ggplot2)
library(stringr)
library(ggrepel)
library(scales)
require(plyr)
library(dplyr)
library(ggpubr)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)

# chromosome      start   end     allele_1        allele_2        ref     total_coverage  ref_coverage    fraction_ref    p.value significant     neg_log_10_p.value      SNP_pair        cell    annotation

subsetColumn <- function(d, colname, factorlist){
    d <- d[d[,colname] %in% factorlist,]
    return(d)
}

filterMaf <- function(d){
    dtemp <- data.frame(str_split_fixed(d$altFreq, ",", 2))
    colnames(dtemp) <- c('f1','f2')
    dtemp$f2 <- as.numeric(as.character(dtemp$f2))
    dtemp$maxval <- apply(dtemp, 1, max, na.rm=TRUE)
    d <- cbind(d, dtemp)

    d <- subset(d, maxval >= 0.2 & maxval <= 0.8)  # index SNP MAF >= 0.2
    return(d)
}

anovaTest <- function(d){
    stat <- lm(formula = deviation ~ annotation, data=d)
    return(anova(stat))
}

wilcoxTest <- function(d){
    stat <- pairwise.wilcox.test(d$deviation, d$annotation, paired=FALSE, p.adjust="none", na.rm=TRUE)
    return(stat)
}

binomTest <- function(x, n, p){binom.test(x, n, p=p, alternative="two.sided", conf.level=0.95)}

renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}

makePlot <- function(d, filter){
    d <- data.frame(d %>% group_by(annotation) %>% mutate(count = n()))
    d$annotation <- paste(d$annotation, "\n(N = ", d$count, ")", sep="") 
    p <- ggplot(d, aes(x = annotation, y = deviation)) +
        geom_boxplot(aes(fill=annotation), width=0.5, outlier.shape=NA) +
        labs(y="Effect size (Allelic Bias)", x="", title=paste("coverage filter", filter, sep=" ")) +
        scale_fill_brewer(name="", palette="Set1", guide=FALSE) +
        theme(axis.text.x=element_text(size=7, angle=30, hjust=1), strip.text.x = element_text(size = 8), panel.background = element_rect(fill = 'white', colour='black'), panel.grid=element_blank())
    return(p)
}

nominalPval_coverageFilters <- function(){
    for ( filter in c(10, 15,20)){
        d <- subsetColumn(d, "significant", c("yes"))   # Nominal pval
        d1 <- subset(d, (total_coverage >= filter)) ## & (ref_coverage >= 2 ) & ((total_coverage - ref_coverage) >= 2))
        print(filter)
        print(wilcoxTest(d1))
        print(makePlot(renameAnnotations(d1), filter))
    }
}

fdrBY_coverageFileters <- function(){
    for ( filter in c('qval_10')){
        d1 <- d[ d[[filter]] <= 0.05 ,]
        d1 <- d1[complete.cases(d1[,filter]),]
        print(filter)
        print(wilcoxTest(d1))
        print(makePlot(renameAnnotations(d1), filter))
    }
}


d <- subsetColumn(d, "cell", c("GM12878"))
d <- subsetColumn(d, "annotation", c('hotRegions','broadDomains','stretchEnhancer','superEnhancers','typicalEnhancers'))

d$deviation <- abs(d$expectedFracRef - d$fraction_ref)
d <- filterMaf(d)
d <- d[!duplicated(d[,c("chrom","end","annotation")]),]

pdf(args[2], height=2.5, width=3)
fdrBY_coverageFileters()
dev.off()



testBinomPval <- function(){
    ndf <- t(data.frame(mapply(binomTest, d$ref_coverage, d$total_coverage, d$expectedFracRef)))[,c('statistic','parameter','p.value')]
    colnames(ndf) <- c('statistic','parameter','newPval')
    d <- cbind(d, ndf)
    print(d[,c('p.value','newPval')])
}
