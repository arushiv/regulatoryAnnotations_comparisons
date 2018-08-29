library(ggplot2)
library(reshape2)
args <- commandArgs(TRUE)
d <- read.table(snakemake@input[[1]], header=T)

# cell1 annotation cell2 chromatinState annotation_length stateOverlapLength"

renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}

makePlot <- function(d){
    p <- ggplot(d, aes(x = cell2, y=fraction, fill=chromatinState)) +
        geom_bar(colour = "black", stat = "identity") +
        facet_grid(cell1 ~ annotation, labeller = label_wrap_gen()) +
        theme(axis.text = element_text(size=7),
          text = element_text(size=7),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(),
          legend.position = "bottom",
          strip.text = element_text(size=8),
          axis.text.x = element_text(size=7, angle=60, hjust=1),
          axis.title=element_text(size=9),
          strip.background=element_rect(fill="white", colour="black")) +
        scale_fill_manual(values = colvector, guide = guide_legend(reverse = TRUE), name="Chromatin States") +
        labs(x = "Chromatin States defined in cell type",
             y = "Fraction of annotations overlapped by chromatin states")
    return(p)
}


d$fraction <- d$stateOverlapLength/d$annotation_length

d$chromatinState <- factor(d$chromatinState, levels = rev(c("1_Active_TSS",
                                                            "2_Weak_TSS",
                                                            "3_Flanking_TSS",
                                                            "5_Strong_transcription",
                                                            "6_Weak_transcription",
                                                            "8_Genic_enhancer",
                                                            "9_Active_enhancer_1",
                                                            "10_Active_enhancer_2",
                                                            "11_Weak_enhancer",
                                                            "14_Bivalent_poised_TSS",
                                                            "16_Repressed_polycomb",
                                                            "17_Weak_repressed_polycomb",
                                                            "18_Quiescent_low_signal")))

colvector <- rev(c(rgb(255,0,0,maxColorValue=255),
                   rgb(255,69,0,maxColorValue=255),
                   rgb(255,69,0,maxColorValue=255),
                   rgb(0,128,0,maxColorValue=255),
                   rgb(0,100,0,maxColorValue=255),
                   rgb(194,225,5,maxColorValue=255),
                   rgb(255,195,77,maxColorValue=255),
                   rgb(255,195,77,maxColorValue=255),
                   rgb(255,255,0,maxColorValue=255),
                   rgb(205,92,92,maxColorValue=255),
                   rgb(128,128,128,maxColorValue=255),
                   rgb(192,192,192,maxColorValue=255),
                   rgb(255,255,255,maxColorValue=255)))

d$cell1 <- paste(d$cell1, "annotations", sep=" ")
## d$cell2 <- paste(d$cell2, "chromatin states", sep=" ")
d$cell2 <- factor(d$cell2, levels = c("K562","HepG2","H1","GM12878"))
d$cell1 <- factor(d$cell1, levels = c("K562 annotations","HepG2 annotations","H1 annotations","GM12878 annotations"))
d <- renameAnnotations(d)

pdf(snakemake@output[[1]], height = 8, width = 7)
## makePlot(d)
makePlot(d)
dev.off()                               


