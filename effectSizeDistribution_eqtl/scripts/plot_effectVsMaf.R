library(ggplot2)
library(reshape2)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)

                                        # gene_id gene_name       gene_chr        gene_start      gene_end        strand  num_var beta_shape1     beta_shape2     true_df pval_true_df    variant_id      tss_distance    chr     pos     ref     alt     num_alt_per_site        rs_id_dbSNP147_GRCh37p13        minor_allele_samples    minor_allele_count      maf     ref_factor      pval_nominal    slope   slope_se        pval_perm       pval_beta       qval    pval_nominal_threshold

d$absEffect <- abs(d$slope)
## d$minuslogpval <- -log10(d$pval)

densityEffect <- function(d){
    p <- ggplot(d, aes(absEffect), cex.lab=1) +
        geom_density() +
        theme(panel.background=element_rect(fill="white", color="black"), panel.grid=element_blank())
    return(p)
}

makeplot <- function(){
    p <- ggplot(d, aes(y=absEffect, x=maf), cex.lab=1) +
        geom_point(size=1, alpha=0.4) +
        theme(panel.background=element_rect(fill="white", color="black"), panel.grid=element_blank())
    return(p)
}



pdf(args[2], height=6, width=6)
densityEffect(d)
makeplot()
dev.off()
