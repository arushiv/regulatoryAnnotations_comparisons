library(ggplot2)
library(scales)
library(gridExtra)
library(viridis)
args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T, sep='\t')

select = c('Alzheimer_disease_and_age_of_onset',
          'Autism_spectrum_disorder_or_schizophrenia',
          'Bipolar_disorder',
          'Bipolar_disorder_and_schizophrenia',
          'Cognitive_decline_rate_in_late_mild_cognitive_impairment',
          'Alcohol_dependence',
          'Alcoholic_chronic_pancreatitis',
          'Educational_attainment',
          'Major_depressive_disorder',
          'Migraine',
          'Neuroticism',
          'Night_sleep_phenotypes',
          'Parkinson_s_disease',
          'Schizophrenia',
          'Subjective_well-being',

          'Allergic_disease__asthma',
          'Asthma',
          'Celiac_disease',
          'Chronic_inflammatory_diseases__ankylosing_spondylitis',
          'Crohn_s_disease',
          'Inflammatory_bowel_disease',
          'Psoriasis',
          'Rheumatoid_arthritis',
          'Systemic_lupus_erythematosus',
          'Ulcerative_colitis',
          'Type_1_diabetes',

          'BMI',
          'Cholesterol',
          'LDL_cholesterol',
          'HDL_cholesterol',
          'Hip_circumference',
          'Obesity',
          'Waist_circumference',
          'Waist-to-hip_ratio',
          'Peripheral_arterial_disease',
          'Trans_fatty_acid_levels',
          'Triglycerides',
          'Type_2_diabetes',

          'Breast_cancer',
          'Lung_cancer',
          'Lung_cancer_in_ever_smokers',
          'Chronic_lymphocytic_leukemia',
          'Colorectal_cancer',
          'Prostate_cancer',
          'Squamous_cell_lung_carcinoma',
          'Testicular_germ_cell_tumor',


          'Blood_metabolite_levels',
          'Blood_protein_levels',
          'Diastolic_blood_pressure',
          'Eosinophil',
          'Fibrinogen_levels',
          'Granulocyte',
          'Hematocrit',
          'Hemoglobin_concentration',
          'High_light_scatter_reticulocyte',
          'Immature_fraction_of_reticulocytes',
          'Lymphocyte_counts',
          'Lymphocyte_percentage_of_white_cells',
          'Mean_corpuscular_hemoglobin',
          'Mean_corpuscular_volume',
          'Mean_platelet_volume',
          'Monocyte_count',
          'Monocyte_percentage_of_white_cells',
          'Myeloid_white_cell_count',
          'Neutrophil',
          'Platelet',
          'Post_bronchodilator_FEV1',
          'Pulse_pressure',
          'Red_blood_cell_count',
          'Red_cell_distribution_width',
          'Reticulocyte_count',
          'Reticulocyte_fraction_of_red_cells',
          'Sum_basophil_neutrophil_counts',
          'Sum_eosinophil_basophil_counts',
          'Sum_neutrophil_eosinophil_counts',
          'Systolic_blood_pressure',
          'White_blood_cell_count',

          'Menarche',
          'Height',
          'Mosquito_bite_size',
          'Multiple_sclerosis',
          'Bone_mineral_density',
          'Total_body_bone_mineral_density'
)
renameAnnotations <- function(d){
    d$annotation <- gsub("hotRegions","HOT Regions",d$annotation)
    d$annotation <- gsub("broadDomains","Broad Domains",d$annotation)
    d$annotation <- gsub("stretchEnhancer","Stretch Enhancers",d$annotation)
    d$annotation <- gsub("superEnhancers","Super Enhancers",d$annotation)
    d$annotation <- gsub("typicalEnhancers","Typical Enhancers",d$annotation)
    return(d)
}

makePlot <- function(d){
    p <- ggplot(d, aes(x=trait, y=Enrich)) +
        geom_point(aes(fill=annotation), size=2, shape=21, colour="black") +
        geom_line(aes(group=annotation, colour=annotation)) +
        labs(y="eQTL Fold enrichment", x="trait") +
        scale_fill_brewer(name="", palette="Set1") +
        scale_colour_brewer(name="", palette="Set1") +
        facet_wrap(~cell, ncol=1) +
        theme(axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5),
              strip.text.x = element_text(size = 9),
              panel.background = element_rect(fill = 'white', colour='black'),
              panel.grid=element_blank(),
              legend.position="bottom")
    return(p)
}

makeHeatmap = function(d, fillstring){
    maxE = max(d$enrichment, na.rm=TRUE)
    minE = min(d$enrichment, na.rm=TRUE)
    print(minE)
    colvector = c(0, 0.33*maxE, 0.66*maxE, maxE)
    print(colvector)
    p = ggplot(subset(d, pval <= bonfStringent), aes(y=trait, x=annotation)) +
        geom_tile(aes_string(fill=fillstring), colour="black") +
                                        #geom_tile(data = d[d$pval > 0.05,], aes(fill=enrichment), fill="grey") +
        facet_wrap(~cell, nrow=1) +
        theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=7),
              axis.text.y = element_text(size=6),
              legend.position="bottom",
              panel.grid = element_blank()) +
        scale_fill_viridis()
        ## geom_tile(data=d[d$pval>0.05,], fill="grey") +
        ## scale_fill_gradientn(colours=c("white","yellow", "red", "darkred"), values=rescale(colvector), breaks=c(colvector), labels=c(round(colvector,2)), name="Fold enrichment")
    return(p)
}

getcols <- function(d){
    p = ggplot(subset(d, pval <= bonfStringent), aes(y=trait, x=type)) +
        geom_tile(aes(fill=related)) +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=9),
              axis.ticks.y = element_blank()) +
        facet_wrap(~type, nrow=1)
    return(p)
}

d <- subset(d, ! annotation %in% c("ActiveTss","typicalEnhancerChromatin") )
d <- renameAnnotations(d)
d$trait <- factor(d$trait, levels=select)
d$type <- "Broad classification of trait"
d$minuslogpval <- -log10(d$pval)

bonfAnnot <- 0.05/5
bonfTrait <- 0.05/83
bonfStringent <- 0.05/(83*5)

p1 = makeHeatmap(d, "enrichment")
p2 = getcols(d)

p3 = makeHeatmap(d, "minuslogpval")
p4 = makeHeatmap(d, "overlap")

pdf(args[2], height=11, width=8)
grid.arrange(p1,p2, nrow=1, heights=c(1.0), widths=c(.80,.20))
grid.arrange(p3,p2, nrow=1, heights=c(1.0), widths=c(.80,.20))
grid.arrange(p4,p2, nrow=1, heights=c(1.0), widths=c(.80,.20))
dev.off()

