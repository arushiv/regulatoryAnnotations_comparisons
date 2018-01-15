library(ggplot2)
library(stringr)
library(ggrepel)
library(scales)
require(plyr)
library(dplyr)

args <- commandArgs(TRUE)
## d <- read.table(args[1], header=T)
d <- read.table(snakemake@input[[1]], header=T)

# chromosome      start   end     allele_1        allele_2        ref     total_coverage  ref_coverage    fraction_ref    p.value significant     neg_log_10_p.value      SNP_pair   expectedFracRef

anovaTest <- function(d){
    stat <- lm(formula = deviation ~ annotation, data=d)
    return(anova(stat))
}

wilcoxTest <- function(d){
    stat <- pairwise.wilcox.test(d$deviation, d$annotation, paired=FALSE, p.adjust="BY")
    return(stat)
}

binomTest <- function(x, n, p){binom.test(x, n, p=p, alternative="two.sided", conf.level=0.95)}


makeQvalColumns <- function(d, filterList){
    for (filter in filterList){
        passed <- which(d$total_coverage >= filter)
        name <- paste("qval",filter,sep="_")
        d[[name]] = NA
        d[[name]][passed] = p.adjust(d$p.value[passed], method="BY") 
    }
    return(d)
}


filterList <- c(10, 15, 20)
d1 <- makeQvalColumns(d, filterList)

## write.table(d1, file = args[2], sep='\t', quote=FALSE)
write.table(d1, file = snakemake@output[[1]], sep='\t', quote=FALSE, row.names=FALSE)    
