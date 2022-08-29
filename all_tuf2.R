#####----- Load all packages ----
BiocManager::install("metagenomeSeq", force= TRUE)
BiocManager::install("biomformat", force = TRUE)
BiocManager::install("glmnet", force = TRUE)
BiocManager::install("Matrix", force = TRUE)
BiocManager::install("ggplot2", force = TRUE)
BiocManager::install("miaViz", force = TRUE)
install.packages("remotes")
remotes::install_github("vmikk/metagMisc", force=TRUE)
library ("metagMisc")
library(metagenomeSeq); packageVersion("metagenomeSeq")
library(biomformat); packageVersion("biomformat")
library (phyloseq); packageVersion("phyloseq")
library (ggplot2)
library (miaViz)
library (Matrix)
library (glmnet)

####--- Create phyloseq object ---##

#Transform into matrixes otu and tax tables (sample table can be left as data frame)
#add taxa 
staph_taxa3 <- read.csv("~/Desktop/Rproject/updated_staph.txt", stringsAsFactors=FALSE,sep="\t", row.names = 1)
staph_OTU3 <- read.csv("~/Desktop/Rproject/all_tuf2_472.txt",stringsAsFactors=FALSE,sep="\t", row.names = 1)
staph_meta3 <- read.csv("~/Desktop/Rproject/meta_tuf2.txt",stringsAsFactors=FALSE,sep="\t", row.names = 1)
ncol (staph_OTU3)
View (staph_meta3)
nrow(staph_taxa3)
View (staph_taxa3)
#ZHlist<-read.table("~/Desktop/Rproject/trial_ZH.txt",sep="\t", header= TRUE)
##View(ZHlist)
#try6 <- colnames(staph_OTU2)
#write.table(try6, file = "~/Desktop/RProject/ZH_staph.txt", col.names=TRUE, sep = "\t", quote = FALSE)

#View (staph_OTU2)
##subset new metadata from the original 
#try <- which(rownames(final_MD) %in% colnames(staph_OTU2))
#try3 <- rownames (staph_meta2)
#staph_OTU.new <- staph_OTU2[ ,try3 ]


#staph_meta2<-final_MD[colnames(staph_OTU.new), ]
#View (staph_meta2)

#write.table(try3, file = "~/Desktop/RProject/ZH_updated_staph.txt", col.names=TRUE, sep = "\t", quote = FALSE)
#nrow(staph_taxa2)

#Transform to phyloseq objects
otu_mat3 <- as.matrix(staph_OTU3)
ncol (otu_mat3)
View (otu_mat3)
tax_mat3 <- as.matrix(staph_taxa3)
OTU3 = otu_table(otu_mat3, taxa_are_rows = TRUE)
TAX3 = tax_table(tax_mat3)
samples = sample_data(staph_meta3)
nrow (samples)
ncol (OTU3)
nrow(TAX3)
ps.staph3 <- phyloseq(OTU3, TAX3, samples)

#####---Filtering ------
#prevalence filter
#remove low prevalence ASV's (2o%)
#cleaned_ps.staph2<-filter_taxa(ps.staph2, function(x) sum(x) > .5, TRUE)
cleaned_ps.staph3 <- phyloseq_filter_prevalence(ps.staph3, prev.trh = 0.2, abund.trh = NULL, threshold_condition = "OR", abund.type = "total")
taxa_sums (cleaned_ps.staph3)

#prune ASVs with 4 or less.
cleaned_ps.staph3 <-prune_taxa(taxa_sums(cleaned_ps.staph3) > 4, cleaned_ps.staph3) 
tail (taxa_sums (cleaned_ps.staph2))
sort(taxa_sums(cleaned_ps.staph2), decreasing=FALSE)


staph_final_MD <- sample_data(cleaned_ps.staph3)
nrow (staph_final_MD)
## export 
write.table(staph_final_MD, file = "~/Desktop/Rproject/staph_MD_final.txt", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

##export_ps
export_staph_ps <- phyloseq_to_df(cleaned_ps.staph3)
write.table(export_staph_ps, file = "~/Desktop/Rproject/staph_all.txt", col.names=TRUE, row.names = FALSE, sep = "\t", quote = FALSE) 


####------Subsetting the data 8 ways------------####

input_table <- read.csv(file="~/Desktop/Rproject/staph_clr.txt", sep="\t", row.names = 1)
metadata <- read.csv(file="~/Desktop/Rproject/staphMD2.txt", sep="\t", row.names = 1)
nrow (metadata)

#####-----V1 VS V2 ------####
#Responders
R_V12 <- which(metadata$Response=="Improved")

subset <- rownames(metadata[R_V12,])
subset_MD <- metadata[R_V12,]
#ASV table 
ASV_input <- input_table %>% filter(row.names(input_table) %in% subset)
View (input_table)
nrow (metadata)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input$Staphylococcus_epidermidis,sample_data(subset_MD)$Visit, p.adjust.method= "fdr", paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input$Staphylococcus_aureus,sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input$Staphylococcus_hominis,sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input$Staphylococcus_capitis,sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input$Staphylococcus_lugdunensis,sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input$Staphylococcus_haemolyticus,sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input$Staphylococcus_pettenkoferi,sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input$Staphylococcus_warneri,sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
#Non-responders
NR_V12 <- which(metadata$Treatment_Response=="nbUVB_Worsened")

subset <- rownames(metadata[NR_V12,])
subset_MD <- metadata[NR_V12,]
#ASV table 
ASV_input <- input_table %>% filter(row.names(input_table) %in% subset)

# Next, run a t-test just like with the alpha diversity - in this case, paired
pairwise.wilcox.test(ASV_input$Staphylococcus_epidermidis, sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = TRUE, exact = FALSE)
pairwise.wilcox.test(ASV_input$Staphylococcus_aureus, sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = TRUE, exact = FALSE)
pairwise.wilcox.test(ASV_input$Staphylococcus_hominis, sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = TRUE, exact = FALSE)
pairwise.wilcox.test(ASV_input$Staphylococcus_capitis, sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = TRUE, exact = FALSE)
pairwise.wilcox.test(ASV_input$Staphylococcus_lugdunensis, sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = TRUE, exact = FALSE)
pairwise.wilcox.test(ASV_input$Staphylococcus_haemolyticus, sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = TRUE, exact = FALSE)
pairwise.wilcox.test(ASV_input$Staphylococcus_pettenkoferi, sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = TRUE, exact= FALSE)
pairwise.wilcox.test(ASV_input$Staphylococcus_warneri, sample_data(subset_MD)$Visit, p.adjust.method= "fdr",paired = TRUE, exact = FALSE)



##------- Les Responders V1 - V2-----
Les_R_V12 <- which(metadata$lesional_nonlesional=="lesional" & metadata$Response=="Improved")

subset1 <- rownames(metadata[Les_R_V12,])
subset1_MD <- metadata[Les_R_V12,]
#ASV table 
ASV_input1 <- input_table %>% filter(row.names(input_table) %in% subset1)
View (input_table)
nrow (metadata)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input1$Staphylococcus_epidermidis,sample_data(subset1_MD)$Visit, p.adjust.method= "fdr", paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input1$Staphylococcus_aureus,sample_data(subset1_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input1$Staphylococcus_hominis,sample_data(subset1_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input1$Staphylococcus_capitis,sample_data(subset1_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input1$Staphylococcus_lugdunensis,sample_data(subset1_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input1$Staphylococcus_haemolyticus,sample_data(subset1_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input1$Staphylococcus_pettenkoferi,sample_data(subset1_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
pairwise.wilcox.test(ASV_input1$Staphylococcus_warneri,sample_data(subset1_MD)$Visit, p.adjust.method= "fdr",paired = FALSE, exact = F)
####----Les Non Res V1 - V2-----#####

Les_NR_V12 <- which(metadata$lesional_nonlesional=="lesional" & metadata$Treatment_Response=="nbUVB_Worsened")

subset2 <- rownames(metadata[Les_NR_V12,])
subset2_MD <- metadata[Les_NR_V12,]
#ASV table 
ASV_input2 <- input_table %>% filter(row.names(input_table) %in% subset2)

# Next, run a t-test just like with the alpha diversity - in this case, paired
pairwise.wilcox.test(ASV_input2$Staphylococcus_epidermidis, sample_data(subset2_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input2$Staphylococcus_aureus, sample_data(subset2_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input2$Staphylococcus_hominis, sample_data(subset2_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input2$Staphylococcus_capitis, sample_data(subset2_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input2$Staphylococcus_lugdunensis, sample_data(subset2_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input2$Staphylococcus_haemolyticus, sample_data(subset2_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input2$Staphylococcus_pettenkoferi, sample_data(subset2_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input2$Staphylococcus_warneri, sample_data(subset2_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
####----- Non-les Responders V1V2-----####
NonLes_R_V12 <- which(metadata$lesional_nonlesional=="non_lesional" & metadata$Treatment_Response=="nbUVB_Improved")

subset3 <- rownames(metadata[NonLes_R_V12,])
subset3_MD <- metadata[NonLes_R_V12,]
#ASV table 
ASV_input3 <- input_table %>% filter(row.names(input_table) %in% subset3)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input3$Staphylococcus_epidermidis, sample_data(subset3_MD)$Visit, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input3$Staphylococcus_aureus, sample_data(subset3_MD)$Visit, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input3$Staphylococcus_hominis, sample_data(subset3_MD)$Visit, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input3$Staphylococcus_capitis, sample_data(subset3_MD)$Visit, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input3$Staphylococcus_lugdunensis, sample_data(subset3_MD)$Visit, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input3$Staphylococcus_haemolyticus, sample_data(subset3_MD)$Visit, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input3$Staphylococcus_pettenkoferi, sample_data(subset3_MD)$Visit, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input3$Staphylococcus_warneri, sample_data(subset3_MD)$Visit, p.adjust.method= "fdr",paired = FALSE)
####--- NonLes Non- Responders  V1V2------#####

NonLes_NR_V12 <- which(metadata$lesional_nonlesional=="non_lesional" & metadata$Treatment_Response=="nbUVB_Worsened")

subset4 <- rownames(metadata[NonLes_NR_V12,])
subset4_MD <- metadata[NonLes_NR_V12,]
#ASV table 
ASV_input4 <- input_table %>% filter(row.names(input_table) %in% subset4)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input4$Staphylococcus_epidermidis, sample_data(subset4_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)

pairwise.wilcox.test(ASV_input4$Staphylococcus_aureus, sample_data(subset4_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input4$Staphylococcus_hominis, sample_data(subset4_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input4$Staphylococcus_aureus, sample_data(subset4_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input4$Staphylococcus_capitis, sample_data(subset4_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input4$Staphylococcus_lugdunensis, sample_data(subset4_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input4$Staphylococcus_haemolyticus, sample_data(subset4_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input4$Staphylococcus_pettenkoferi, sample_data(subset4_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
pairwise.wilcox.test(ASV_input4$Staphylococcus_warneri, sample_data(subset4_MD)$Visit, p.adjust.method= "fdr",paired = TRUE)
####----V1 les- R vs NR-----
V1Les <- which(metadata$lesional_nonlesional=="lesional" & metadata$Visit=="1" & metadata$Treatment=="nbUVB")

subset5 <- rownames(metadata[V1Les,])
subset5_MD <- metadata[V1Les,]
#ASV table 
ASV_input5 <- input_table %>% filter(row.names(input_table) %in% subset5)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input5$Staphylococcus_epidermidis, sample_data(subset5_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input5$Staphylococcus_aureus, sample_data(subset5_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input5$Staphylococcus_hominis, sample_data(subset5_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input5$Staphylococcus_capitis, sample_data(subset5_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input5$Staphylococcus_lugdunensis, sample_data(subset5_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input5$Staphylococcus_haemolyticus,sample_data(subset5_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input5$Staphylococcus_pettenkoferi,sample_data(subset5_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input5$Staphylococcus_warneri,sample_data(subset5_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
####----V1 Nonles- R vs NR-----
V1NonLes <- which(metadata$lesional_nonlesional=="non_lesional" & metadata$Visit=="1" & metadata$Treatment=="nbUVB")

subset6 <- rownames(metadata[V1NonLes,])
subset6_MD <- metadata[V1NonLes,]
#ASV table 
ASV_input6 <- input_table %>% filter(row.names(input_table) %in% subset6)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input6$Staphylococcus_epidermidis, sample_data(subset6_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input6$Staphylococcus_aureus, sample_data(subset6_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input6$Staphylococcus_hominis, sample_data(subset6_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input6$Staphylococcus_capitis, sample_data(subset6_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input6$Staphylococcus_lugdunensis, sample_data(subset6_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input6$Staphylococcus_haemolyticus,sample_data(subset6_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input6$Staphylococcus_pettenkoferi,sample_data(subset6_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input6$Staphylococcus_warneri,sample_data(subset6_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
######--------Visit 2-----------######
####----V2 les- R vs NR-----
V2Les <- which(metadata$lesional_nonlesional=="lesional" & metadata$Visit=="2" & metadata$Treatment=="nbUVB")

subset7 <- rownames(metadata[V2Les,])
subset7_MD <- metadata[V2Les,]
#ASV table 
ASV_input7 <- input_table %>% filter(row.names(input_table) %in% subset7)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input7$Staphylococcus_epidermidis, sample_data(subset7_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input7$Staphylococcus_aureus, sample_data(subset7_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input7$Staphylococcus_hominis, sample_data(subset7_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input7$Staphylococcus_capitis, sample_data(subset7_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input7$Staphylococcus_lugdunensis, sample_data(subset7_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input7$Staphylococcus_haemolyticus, sample_data(subset7_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input7$Staphylococcus_pettenkoferi, sample_data(subset7_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input7$Staphylococcus_warneri, sample_data(subset7_MD)$Response, p.adjust.method= "fdr",paired = FALSE)

####----V2 Nonles- R vs NR-----
V2NonLes <- which(metadata$lesional_nonlesional=="non_lesional" & metadata$Visit=="2" & metadata$Treatment=="nbUVB")

subset8 <- rownames(metadata[V2NonLes,])
subset8_MD <- metadata[V2NonLes,]
#ASV table 
ASV_input8 <- input_table %>% filter(row.names(input_table) %in% subset8)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input8$Staphylococcus_epidermidis, sample_data(subset8_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input8$Staphylococcus_aureus, sample_data(subset8_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input8$Staphylococcus_hominis, sample_data(subset8_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input8$Staphylococcus_capitis, sample_data(subset8_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input8$Staphylococcus_lugdunensis, sample_data(subset8_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input8$Staphylococcus_haemolyticus, sample_data(subset8_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input8$Staphylococcus_pettenkoferi, sample_data(subset8_MD)$Response, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input8$Staphylococcus_warneri, sample_data(subset8_MD)$Response, p.adjust.method= "fdr",paired = FALSE)

####----V1 Not treated vs Treated -----#####
V1NT_TR <- which(metadata$Treatment_Response=="nbUVB_Improved" & metadata$Visit=="1" |metadata$Visit=="1" & metadata$Treatment=="No_nbUVB")

subset9 <- rownames(metadata[V1NT_TR,])
subset9_MD <- metadata[V1NT_TR,]
#ASV table 
ASV_input9 <- input_table %>% filter(row.names(input_table) %in% subset9)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input9$Staphylococcus_epidermidis, sample_data(subset9_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input9$Staphylococcus_aureus, sample_data(subset9_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input9$Staphylococcus_hominis, sample_data(subset9_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input9$Staphylococcus_capitis, sample_data(subset9_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input9$Staphylococcus_lugdunensis, sample_data(subset9_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input9$Staphylococcus_haemolyticus, sample_data(subset9_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input9$Staphylococcus_pettenkoferi, sample_data(subset9_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input9$Staphylococcus_warneri, sample_data(subset9_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
## Non responder 
V1NT_TNR <- which(metadata$Treatment_Response=="nbUVB_Worsened" & metadata$Visit=="1" |metadata$Visit=="1" & metadata$Treatment=="No_nbUVB")

subset10 <- rownames(metadata[V1NT_TNR,])
subset10_MD <- metadata[V1NT_TNR,]
#ASV table 
ASV_input10 <- input_table %>% filter(row.names(input_table) %in% subset10)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input10$Staphylococcus_epidermidis, sample_data(subset10_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input10$Staphylococcus_aureus, sample_data(subset10_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input10$Staphylococcus_hominis, sample_data(subset10_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input10$Staphylococcus_capitis, sample_data(subset10_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input10$Staphylococcus_lugdunensis, sample_data(subset10_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input10$Staphylococcus_haemolyticus, sample_data(subset10_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input10$Staphylococcus_pettenkoferi, sample_data(subset10_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input10$Staphylococcus_warneri, sample_data(subset10_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)

####----V2 Not treated vs Treated -----#####
V2NT_TR <- which(metadata$Treatment_Response=="nbUVB_Improved" & metadata$Visit=="2" |metadata$Visit=="2" & metadata$Treatment=="No_nbUVB")

subset11 <- rownames(metadata[V2NT_TR,])
subset11_MD <- metadata[V2NT_TR,]
#ASV table 
ASV_input11 <- input_table %>% filter(row.names(input_table) %in% subset11)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input11$Staphylococcus_epidermidis, sample_data(subset11_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input11$Staphylococcus_aureus, sample_data(subset11_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input11$Staphylococcus_hominis, sample_data(subset11_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input11$Staphylococcus_capitis, sample_data(subset11_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input11$Staphylococcus_lugdunensis, sample_data(subset11_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input11$Staphylococcus_haemolyticus, sample_data(subset11_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input11$Staphylococcus_pettenkoferi, sample_data(subset11_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input11$Staphylococcus_warneri, sample_data(subset11_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
## Non responder 
V2NT_TNR <- which(metadata$Treatment_Response=="nbUVB_Worsened" & metadata$Visit=="2" |metadata$Visit=="2" & metadata$Treatment=="No_nbUVB")

subset12 <- rownames(metadata[V2NT_TNR,])
subset12_MD <- metadata[V2NT_TNR,]
#ASV table 
ASV_input12 <- input_table %>% filter(row.names(input_table) %in% subset12)

# Next, run a t-test just like with the alpha diversity - in this case, paired

pairwise.wilcox.test(ASV_input12$Staphylococcus_epidermidis, sample_data(subset12_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input12$Staphylococcus_aureus, sample_data(subset12_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input12$Staphylococcus_hominis, sample_data(subset12_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input12$Staphylococcus_capitis, sample_data(subset12_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input12$Staphylococcus_lugdunensis, sample_data(subset12_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input12$Staphylococcus_haemolyticus, sample_data(subset12_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input12$Staphylococcus_pettenkoferi, sample_data(subset12_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)
pairwise.wilcox.test(ASV_input12$Staphylococcus_warneri, sample_data(subset12_MD)$Treatment, p.adjust.method= "fdr",paired = FALSE)

delta_mswat <-read.csv("~/Desktop/Rproject/try_mswat.txt", stringsAsFactors=FALSE,sep="\t", header= TRUE)
plot (delta_mswat)
dotchart(delta_mswat$delta_mSWAT, labels = delta_mswat$Patient_ID, xlab= "Patient_ID", ylab = "delta_mswat", pch = 21, bg = "green", pt.cex = 1.5)

hist(delta_mswat$delta_mSWAT,xlab = "Patient_ID",col = "yellow",border = "blue")

