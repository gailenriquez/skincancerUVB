load("/Users/gailhillary/Desktop/WorspaceforGail.RData")## Installing dada2 31 Aug 2021
 if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
##install and load packages
# Update bioC 
BiocManager::install(version = 3.15)
install.packages("tidyverse", force= TRUE)
library("tidyverse"); packageVersion("tidyverse")
BiocManager::install("phyloseq", force= TRUE)
library("phyloseq"); packageVersion("phyloseq")
library("data.table");packageVersion("data.table")
install.packages("ggplot2", force= TRUE)
library("ggplot2"); packageVersion("ggplot2")
install.packages("data.table")
install.packages("tibble", force=TRUE)
install.packages("dplyr", force= TRUE)
library("dplyr")
BiocManager::install("ggpubr", force= TRUE)
library ("ggpubr"); packageVersion ("ggpubr")
# Actually install
BiocManager::install("dada2", version = "3.15")

# The fastq files are somewhat piecemeal and need to be concateneated together and then fished out of their directory structure.
# I'll do this in bash using the following, non-R code. Just putting it here for safekeeping (and so I have a record of how it was done).

# The folders that contain the fastq sequences from BaseSpace contain a - at the start of the folder names. This fucks things up pretty badly. 
# I can fix this with the following bash loop, run from the terminal in the master folder containing the sample folders

# for dir in -476*
# do
# mv "./$dir" "folder$dir"
# done 

# In the loop above, the -476 is just something that all the folders in the directory had at the start of their names. Update as needed.

# Next, I want to uncompress the fastq.gz files to make them human-readable, then merge the relevant ones together.
# Run this code from the Linux or Mac terminal within the master folder downloaded from BaseSpace

# find . -maxdepth 3 -type fastq.gz -execdir gunzip {} \;

# These command have worked to generate the correct output by concatenating the forward reads
# find -s . -maxdepth 3 -name '*R1_00[0-9].fastq' -execdir ls {} \; -execdir bash -c 'cat {} >> forward.merged' \; -execdir bash -c 'i=$(ls {} | cut -d_ -f1,4);j=".fastq"; echo "${i}${j}" > namefile.txt' \;
# find -s . -maxdepth 3 -name forward.merged -execdir bash -c 'x=$(head -n1 namefile.txt); mv {} $x' \;

# Then with reverse reads
# find -s . -maxdepth 3 -name '*R2_00[0-9].fastq' -execdir ls {} \; -execdir bash -c 'cat {} >> reverse.merged' \; -execdir bash -c 'i=$(ls {} | cut -d_ -f1,4);j=".fastq"; echo "${i}${j}" > namefile.txt' \;
# find -s . -maxdepth 3 -name reverse.merged -execdir bash -c 'x=$(head -n1 namefile.txt); mv {} $x' \;

# Then pull these files into the relevant folder
# find -s . -name '*_R[1-2].fastq' -exec cp {} ~/Desktop/Sorted_Fastqs/ \;
# cd ~/Desktop/Sorted_Fastqs

# Next, get rid of the Alan-### prefix so the files all have the same formattting
# rename 's/Alan...-//' *.fastq

# Now the fils should be all ready to be used for adapter trimming.

## Primers used:
# 515F-modified: GTGYCAGCMGCCGCGGTAA
# 806R-modified: GGACTACNVGGGTWTCTAAT
# >E.coli_v4_292_nt_with_primers
# GTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCC
# >E.coli_v4_253_nt_without_primers
# TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGG
# >515F-modified
# GTGYCAGCMGCCGCGGTAA
# >806R-modified
# GGACTACNVGGGTWTCTAAT
# >806R-modified_rev_comp
# ATTAGAWACCCBNGTAGTCC

# I have cutadapt installed already (it's a command line application you can just call, once it's installed)
# The following is a bash for loop to run cutadapt over all the file. It should get run in the folder that contains all the concatenated fastq files.

# for read1 in *R1.fastq
# do
# read2=$(echo $read1 | sed 's/R1.fastq/R2.fastq/')
# out1=$(echo $read1 | sed 's/R1/R1.trimmed/')
# out2=$(echo $read2 |sed 's/R2/R2.trimmed/')
# cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT -o $out1 -p $out2 $read1 $read2 | tee >> cutadapt.log.txt
# done

# That worked, so now the trimmed reads are good to get fed into dada2 

# Figure out where we are
getwd()
inputDir <- getwd()


# Load dada2 and set the path variable
library(dada2)
packageVersion("dada2")
path <- "~/Desktop/Skin_Cancer/Sorted_Fastqs/Trimmed/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.trimmed.fastq and SAMPLENAME_R2.trimmed.fastq
fnFs <- sort(list.files(path, pattern="_R1.trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.trimmed.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Take a gander at the quality profiles of the first 6 forward and reverse fastq files
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnFs[1:6])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
length(fnFs)
length(fnRs)

# Got rid of the trunclen used by default and added a minQ funtion to reduce potential noise. Also adjusted the maxEE to a more stringent c(1,1) rather than the default c(2,2).
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minQ = 10,
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# I'm using way more of the dataset than normal to calculate errors (5e+08 instead of 1e+08). Also, I'm running the randomize approach to avoid potential batch effects.
errF <- learnErrors(filtFs, nbases = 5e+08, randomize = TRUE, multithread=TRUE)
errR <- learnErrors(filtRs, nbases = 5e+08, randomize = TRUE, multithread=TRUE)

# Now we see what the error model looks like.
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# Now merge for F and R reads. I had to reduce the minoverlap as there isn't much to work with. If you don't do this, dada2 will not run beyond this point.
# I suspect this the step that the UofI team got stuck on a used an outside tool.
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, maxMismatch = 0, minOverlap = 8, )

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Now we make our initial ASV table
seqtab <- makeSequenceTable(mergers)

# Count the samples and the ASVs
dim(seqtab)

# Take a look at the length distributions of the ASVs in the table. Based on the primer and V4 region work above, they shold be in the ballpark of 252ish.
table(nchar(getSequences(seqtab)))

# They mostly are, but there's also some garbage that resulted from clearly poor merging that needs to get tossed.
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 252:254]
dim(seqtab2)

# Set to pooled to improve PCR chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="pooled", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Eyeball how much of the dataset was removed as chimeras
sum(seqtab.nochim)/sum(seqtab)

# Start tracking the reads all the way through these first steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
barplot(track)

summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                          filtered=out[,2], dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN), merged=sapply(mergers, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/out[,1]*100, 1))
summary_tab
mean(summary_tab$final_perc_reads_retained)
write.table(summary_tab, "~/Desktop/read-count-tracking.txt", quote=FALSE, sep="\t", col.names=NA)

# Now assign taxonomy using the alternative k-mer appraoch with DECIPHER

library(DECIPHER); packageVersion("DECIPHER")

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("/Volumes/GoogleDrive/My Drive/Loyola/RStudioProjects/dadabases/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=TRUE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

taxa <- taxid

# There are two samples that are unaccounted for in the metadata table that need to be removed
row_names_to_remove<-c("ZH91","ZH90")

# Store the original in case I need it again
seqtab.nochim2 <- seqtab.nochim

# Delete the problematic columns
seqtab.nochim <- seqtab.nochim[!(row.names(seqtab.nochim) %in% row_names_to_remove),]

# Load the things needed for phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

samdf <- read.csv(file="~/Desktop/Analysis_Files/Metadata_v7.txt", sep="\t", header=TRUE, row.names = 1)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Now to run decontam
library(decontam); packageVersion("decontam")

head(sample_data(ps))

# Inspect the library sizes
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point(size =0.5)
#changed this to default geom point
# Now to use the prevalence method to identify potential contaminants


sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

# How many of the ASVs are contaminants?
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)

# Now to re-run this with a more stringent approach:
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev05$contaminant)
which(contamdf.prev05$contaminant)

# Look at these results as phyloseq opbjects
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point(size = 0.05) +
  xlab("Prevalence (Controls)") + ylab("Prevalence (Samples)")

# Make a vector to store the contaminant ASVs
contaminant_asvs <- row.names(contamdf.prev05[contamdf.prev05$contaminant == TRUE, ])

# Now to remove the identified contaminants:
badTaxa = contaminant_asvs
allTaxa = taxa_names(ps)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ex1 = prune_taxa(allTaxa, ps)
# new phyloseq object with just the taxa you kept.
ex1

# Now to filer out the rest of the unassigned and the archaea
ex2 = subset_taxa(ps, domain=="Bacteria")

# Now to do some filtering and pruning

ps <- ex2

# Prune out the air samples
ps_no_env = subset_samples(ps, lesional_nonlesional != "Air")

ps <- ps_no_env
str(ps)
BiocManager::install("microbiome", force= TRUE) 
#not helpful they deprecated the thing I wanted
library(microbiome); packageVersion("microbiome")
# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")
library(microbiome); packageVersion("microbiome")
install.packages("remotes")
library ("remotes")
remotes::install_github("vmikk/metagMisc", force = TRUE)
library(metagMisc); packageVersion("metagMisc")

export_ps <- phyloseq_to_df(ps)
str(export_ps)
ncol(export_ps[170,])
write.table(export_ps, file = "~/Desktop/ps_contaminantsremoved_version_10_jan102022.txt", col.names=TRUE, row.names = FALSE, sep = "\t", quote = FALSE) 
?read.table
#rename and load into the new r-script 

export_MD <- sample_data(ps) 
head(export_MD)
write.table(export_MD, file = "~/Desktop/metadata_v8_contaminantsremoved_jan102022.txt", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
##load into new rscipt, rename 

# Make a simple alpha diversity plot
plot_richness(ps, x="Body_Location", measures=c("Shannon", "Simpson"), color="Status_From_V1_V2")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray", max)

ordinate("list")

plot_ordination(ps.prop, ord.nmds.bray, color="Status_From_V1_V2", title="Bray NMDS")

# Plot the top 20 taxa. 
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Patient_ID", fill="family") + geom_bar(stat="identity")
ps@sam_data

##add libraries needed
install.packages("tibble", force= TRUE)
library ("tibble")
BiocManager::install("DESeq2", force= TRUE)
BiocManager::install("edgeR", force= TRUE)
BiocManager::install("limma", force= TRUE)
BiocManager::install("metagenomeSeq", force= TRUE)
BiocManager::install( "microshades")
BiocManager::install("dplyr", force= TRUE)
library("metagenomeSeq"); packageVersion("metagenomeSeq")
library("edgeR")
library("limma")
library(BiocManager)
library (vegan); packageVersion("vegan")
BiocManager::valid("metagenomeSeq")
BiocManager::valid("edgeR")
BiocManager::valid("limma")
BiocManager::valid("DESeq2")
install.packages("vegan", force= TRUE)
install.packages("ggplot2",force =TRUE)
library (ggplot2)
install.packages("dplyr")
library (dplyr); packageVersion ("dplyr")
install.packages("data.table", force= TRUE)
library (data.table)

install.packages("metacoder")
library(metacoder)
library(vegan)
BiocManager::install("remotes", force= TRUE)
remotes::install_github("microbiome/microbiome")
BiocManager::valid("remotes")


###----------------------  Filter data----------####

#create a histogram to see the read count frequency 
hist(sample_sums(ps), main="Histogram: Read Counts", xlab="Total Reads", 
     border="steelblue", col="skyblue", las=1, breaks=40)


#now look at ASVs, then sort to see the lowest count
head(taxa_sums(ps))
tail (taxa_sums(ps))
sort(taxa_sums(ps), decreasing=FALSE)

#prune ASVs with 4 or less.
ps.f <-prune_taxa(taxa_sums(ps) >= 4, ps) 
tail (taxa_sums (ps.f))
sort(taxa_sums(ps.f), decreasing=FALSE)


#prevalence filter
#remove low prevalence ASV's (1o%)
filtered_ps.f <-filter_taxa(ps.f, function(x) sum(x) > .1, TRUE)
taxa_sums (filtered_ps.f)

#remove samples with >1000 reads
sample_sums (filtered_ps.f)
sort (sample_sums(filtered_ps.f), decreasing = TRUE)
fil_ps.f<-prune_samples(sample_sums(filtered_ps.f) >=1000, filtered_ps.f)
tail(sample_sums(fil_ps.f))

# create histogram to check 
hist(sample_sums(fil_ps.f), main="Histogram: Read Counts", xlab="Total Reads after filtering", 
     border="steelblue", col="skyblue", las=1, breaks=40)

##export new metadata
clean_ps <- fil_ps.f
clean_Meta <- sample_data(clean_ps)
head (clean_Meta)
##write csv 
write.table(clean_Meta, file = "~/Desktop/Rproject/OG_MD_1000.txt", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#export metadata
all_MD <-sample_data (clean_ps)
write.table(all_MD, file = "~/Desktop/Rproject/all_MD.1000.txt", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
###updated Metadata to include other interventions 
latestMD <-read.csv(file="~/Desktop/Rproject/edited_last_MD.txt",sep="\t", header=TRUE, row.names= 1)
sampledata <- sample_data(latestMD)

clean_ps <- merge_phyloseq(clean_ps, sampledata)

clean_Meta <- sample_data(clean_ps)
view(sample_data(clean_ps))
dim(sample_data(clean_ps))

###-----------------Sort metadata (paired or unpaired)-----------------####

##-------prune samples for v1 and v2--------###

v1 <- prune_samples(sample_data(clean_Meta)$Visit=="1", clean_Meta)
v2 <- prune_samples(sample_data(clean_Meta)$Visit=="2", clean_Meta)

##add column that has the lesional_nonlesional values plus patient ID
v1loc <- base::transform(v1, location = paste(Patient_ID,Body_Location, sep = "."))
v1loc.t <-tibble::rownames_to_column(v1loc, "Sample_ID") ## add the ZH columns 
v2loc <- base::transform(v2, location = paste(Patient_ID,Body_Location,sep = "."))
v2loc.t <-tibble::rownames_to_column(v2loc, "Sample_ID") ## add the ZH columns 

#add frequency column to our data frame
v1_merge <- merge(v1loc.t,data.frame(table(location = v1loc.t$location)),by = c("location")) 
                  

v2_merge <- merge(v2loc.t, data.frame(table(location = v2loc$location)), by = c("location"))

##keep samples with a frequency of at least 2
v1_merge_filtered <-filter(v1_merge, Freq >=2)
v2_merge_filtered <-filter(v2_merge, Freq >=2)

##subset lesional and nonlesional from v1
lv1 <- filter(v1_merge_filtered [ , -15], lesional_nonlesional== "lesional")
nv1 <- filter (v1_merge_filtered [ , -15],lesional_nonlesional== "non_lesional")

##subset lesional and nonlesional from v2
lv2 <- filter(v2_merge_filtered[ , -15], lesional_nonlesional== "lesional")

nv2 <- filter (v2_merge_filtered [ , -15],lesional_nonlesional== "non_lesional" )

## removing duplicates in body location
lv1.t <-lv1[!duplicated(lv1$location), ]
nv1.t = nv1[!duplicated(nv1$location), ]
lv2.t = lv2[!duplicated(lv2$location), ]
nv2.t = nv2[!duplicated(nv2$location), ]

## join unique locations
##combine with selected columns; all the functions in the codes below in one line

merge_all_v1 <- inner_join(lv1.t,
          nv1.t %>% dplyr::select(location,lesional_nonlesional, Sample_ID),
          by = "location", all=FALSE)

merge_all_v2 <- inner_join(lv2.t,
                                 nv2.t %>% dplyr::select(location,lesional_nonlesional, Sample_ID),
                                 by = "location", all=FALSE)

##write csvs
write.table(merge_all_v1, file = "~/Desktop/Rproject/v1.1000.txt", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(merge_all_v2, file = "~/Desktop/Rproject/v2.1000.txt", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

###-----------------Compare v1 and v2 -------#### 
##Join v1 and v2
location_v1<-merge_all_v1$location
location_v2 <-merge_all_v2$location

##check which Patient ID+Body locations are present on both Data Frames

index <- merge_all_v1$location %in% location_v2 
m1<- merge_all_v1[index, ]
m1.index <-  m1$location #output of the matched v1 DF 
index2 <- merge_all_v2$location %in% m1.index

m2<-merge_all_v2 [index2, ]

#merge two DFs, remove repeat column
t1 <-inner_join(m1,
           m2,
           by = "location", all=FALSE)

all_merged <- t1[!duplicated(as.list(t1))]# remove duplicated columns
nrow (all_merged)
##write txt file
write.table(all_merged, file = "~/Desktop/Rproject/all_1000.txt", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
dim(all_merged)

#####---- Update PS objects------####

##convert combined v1v2 DF to a PS object; only has all the matching samples (Nonlesional and lesional pairs from V1 and V2)

ps_v1v2 <- clean_ps

##prune samples from ps that match our final table
ps_all_samples <-prune_samples(sample_names(ps_v1v2) %in% c(all_merged$Sample_ID.x.x, all_merged$Sample_ID.y.x, all_merged$Sample_ID.x.y, all_merged$Sample_ID.y.y), ps_v1v2)

final_MD <-sample_data (ps_all_samples)
nrow (final_MD)
write.table(final_MD, file = "~/Desktop/Rproject/final_MD.1000.txt", col.names=TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

##export_ps
export_ps2 <- phyloseq_to_df(ps_all_samples)
str(ps_all_samples)
write.table(export_ps2, file = "~/Desktop/Rproject/ps_all_samples.txt", col.names=TRUE, row.names = FALSE, sep = "\t", quote = FALSE) 

####-----Subset samples -----####

##subset all lesional and non lesional
ps_all_les<- subset_samples(ps_all_samples, lesional_nonlesional=="lesional" )
ps_all_non<-subset_samples(ps_all_samples, lesional_nonlesional=="non_lesional")

##lesional Non responders and responders
ps_all_lesional_R <- subset_samples(ps_all_samples, lesional_nonlesional=="lesional" & Status_From_V1_V2 == "Disease_Improved")
ps_all_lesional_NR <- subset_samples(ps_all_samples, lesional_nonlesional=="lesional" & Status_From_V1_V2 == "Disease_Worsened")
#

ps_all_nonles_R <- subset_samples(ps_all_samples, lesional_nonlesional=="non_lesional" & Status_From_V1_V2 == "Disease_Improved")
ps_all_nonles_NR <- subset_samples(ps_all_samples, lesional_nonlesional=="non_lesional" & Status_From_V1_V2 == "Disease_Worsened")
## subset nbUVB samples from non nbUVB

ps_all_nbUVB <- subset_samples(ps_all_samples, Treatment == "nbUVB")
ps_all_No_nbUVB <- subset_samples(ps_all_samples, Treatment == "No_nbUVB")
  
#Subset nbUVB and non nbUVB - lesional/ non lesional , repsonders/nonresponders/no_change
  
ps_all_les_nbUVB<- subset_samples(ps_all_samples, lesional_nonlesional=="lesional" & Treatment == "nbUVB")
ps_all_les_No_nbUVB<- subset_samples(ps_all_samples, lesional_nonlesional=="lesional" & Treatment== "No_nbUVB")

ps_nbUVB_les_R <- subset_samples(ps_all_samples, lesional_nonlesional=="lesional" & Treatment_Response== "nbUVB_Improved")
ps_nbUVB_les_NR <- subset_samples(ps_all_samples, lesional_nonlesional=="lesional" & Treatment_Response== "nbUVB_Worsened")

ps_nbUVB_non_R <- subset_samples(ps_all_samples, lesional_nonlesional=="non_lesional" & Treatment_Response== "nbUVB_Improved")
ps_nbUVB_non_NR <- subset_samples(ps_all_samples, lesional_nonlesional=="non_lesional" & Treatment_Response== "nbUVB_Worsened")



ps_all_non_nbUVB<-subset_samples(ps_all_samples, lesional_nonlesional=="non_lesional" & Treatment == "nbUVB")
ps_all_non_no_nbUVB<-subset_samples(ps_all_samples, lesional_nonlesional=="non_lesional" & Treatment == "No_nbUVB")
 
ps_all_nbUVB_responders <-subset_samples(ps_all_samples, Treatment_Response == "nbUVB_Improved")
ps_all_nbUVB_nonresponders <-subset_samples(ps_all_samples, Treatment_Response == "nbUVB_Worsened")

ps_all_no_nbUVB_nonresponders <- subset_samples(ps_all_samples, Treatment_Response == "No_nbUVB_Worsened")
ps_all_no_nbUVB_No_Change <- subset_samples(ps_all_samples, Treatment_Response == "No_nbUVB_No_Change")

##all not treated -lesional, nonlesional v1v2
ps_NT_les <- subset_samples(ps_all_No_nbUVB, lesional_nonlesional=="lesional")
ps_NT_nonles <- subset_samples(ps_all_No_nbUVB, lesional_nonlesional=="non_lesional")


##extra

ps_NT_vs_R <-subset_samples(ps_all_samples, Treatment_Response == "nbUVB_Improved" & Visit ==1|Treatment == "No_nbUVB" & Visit ==1)
ps_2NT_vs_R <-subset_samples(ps_all_samples, Treatment_Response == "nbUVB_Improved" & Visit ==2|Treatment == "No_nbUVB" & Visit ==2)


ps_NT_vs_NR <- subset_samples(ps_all_samples, Treatment_Response == "nbUVB_Worsened" & Visit ==1|Treatment == "No_nbUVB" & Visit ==1)
ps_2NT_vs_NR <- subset_samples(ps_all_samples, Treatment_Response == "nbUVB_Worsened" & Visit ==2|Treatment == "No_nbUVB" & Visit ==2)


####------------Pruning Visit 1 from clean ps object----------####

ps_v1 <- clean_ps
## clean up v1 based on v1v2 dataframe 
v1.clean <- subset_samples(ps_all_samples, Visit=="1")
v1.clean_nbUVB <-subset_samples(ps_all_samples, Visit=="1" & Treatment == "nbUVB" )

v1.clean_No_nbUVB <-subset_samples(ps_all_samples, Visit=="1" & Treatment == "No_nbUVB" )

##lesional and non lesional
ps_v1_samples.les <- subset_samples(v1.clean, lesional_nonlesional=="lesional")
ps_v1_lesional_nbUVB<- subset_samples(v1.clean, lesional_nonlesional=="lesional" & Treatment == "nbUVB")
ps_v1_lesional_No_nbUVB<- subset_samples(v1.clean, lesional_nonlesional=="lesional" & Treatment == "No_nbUVB")

ps_v1_samples.nonles <- subset_samples(v1.clean, lesional_nonlesional=="non_lesional")
ps_v1_non_nbUVB<- subset_samples(v1.clean, lesional_nonlesional=="non_lesional" & Treatment == "nbUVB")
ps_v1_non_No_nbUVB<- subset_samples(v1.clean, lesional_nonlesional=="non_lesional" & Treatment == "No_nbUVB")


##Responders and Non Responders of PTX
ps_v1_R <- subset_samples(v1.clean, Response == "Improved")
ps_v1_NR <- subset_samples(v1.clean, Response == "Worsened")

#Responders and Non Responders (Lesional)
ps_v1_L_R <- subset_samples(ps_v1_samples.les, Treatment_Response == "nbUVB_Improved")
ps_v1_L_NR <- subset_samples(ps_v1_samples.les, Treatment_Response == "nbUVB_Worsened")

####---------Pruning Visit 2 from clean ps object-------------####
ps_v2 <- clean_ps

##cleaned up v2 based on combined v1v2 data frame
v2.clean <-subset_samples(ps_all_samples, Visit=="2")

##subset nbUVB and no_nbUVB 
v2.clean_nbUVB <-subset_samples(ps_all_samples, Visit=="2" & Treatment == "nbUVB" )
v2.clean_No_nbUVB <-subset_samples(ps_all_samples, Visit=="2" & Treatment == "No_nbUVB" )

##lesional 
ps_v2_samples.les <- subset_samples(v2.clean,lesional_nonlesional=="lesional")

ps_v2_lesional_nbUVB<- subset_samples(v2.clean, lesional_nonlesional=="lesional" & Treatment == "nbUVB")
ps_v2_lesional_No_nbUVB<- subset_samples(v2.clean, lesional_nonlesional=="lesional" & Treatment == "No_nbUVB")

##non lesional
ps_v2_samples.nonles <- subset_samples(v2.clean, lesional_nonlesional=="non_lesional")

ps_v2_non_nbUVB<- subset_samples(v2.clean, lesional_nonlesional=="non_lesional" & Treatment == "nbUVB")
ps_v2_non_No_nbUVB<- subset_samples(v2.clean, lesional_nonlesional=="non_lesional" & Treatment == "No_nbUVB")

##lesional -non responders/ responders 
ps_v2_Les_R<- subset_samples(ps_v2_samples.les, Treatment_Response == "nbUVB_Improved" )
ps_v2_Les_NR <- subset_samples(ps_v2_samples.les,Treatment_Response == "nbUVB_Worsened")

ps_v2_non_R <- subset_samples(ps_v2_samples.nonles, Treatment_Response == "nbUVB_Improved")
ps_v2_non_NR <- subset_samples(ps_v2_samples.nonles,Treatment_Response == "nbUVB_Worsened")

##treated with nbUVB vs no nbUVB
ps_all_worsened <-  subset_samples(ps_all_samples, Response == "Worsened" )
ps_worsened_nbUVB <- subset_samples(ps_all_samples, Treatment_Response == "nbUVB_Worsened" )
ps_worsened_no_nbUVB <- subset_samples(ps_all_samples, Treatment_Response == "No_nbUVB_Worsened" )

ps_v1_worsened <-  subset_samples(ps_all_samples, Visit=="1" & Response == "Worsened" )
ps_v2_worsened <-  subset_samples(ps_all_samples, Visit=="2" & Response == "Worsened" )

R_vs_NT_NR <- subset_samples(ps_v2_samples.les, Treatment_Response == "No_nbUVB_Worsened" |Treatment_Response=="nbUVB_Improved")

 
####------Rarefy PS objects ------####
set.seed(500)
ps_all_rarefied <- rarefy_even_depth(ps_all_samples, sample.size = min(sample_sums(ps_all_samples)),
                                     rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_all_nbUVB.rarefied <- rarefy_even_depth(ps_all_nbUVB, sample.size = min(sample_sums(ps_all_nbUVB)),
                                      rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_all_No_nbUVB.rarefied <- rarefy_even_depth(ps_all_No_nbUVB, sample.size = min(sample_sums(ps_all_No_nbUVB)),
                                              rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_all_les.rarefied <-rarefy_even_depth(ps_all_les, sample.size = min(sample_sums(ps_all_les)),
                                             rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_all_non.rarefied <-rarefy_even_depth(ps_all_non, sample.size = min(sample_sums(ps_all_non)),
                                             rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_all_less_R.rarefied <- rarefy_even_depth(ps_all_lesional_R, sample.size = min(sample_sums(ps_all_lesional_R)),
                                        rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_all_les_NR.rarefied<- rarefy_even_depth(ps_all_lesional_NR, sample.size = min(sample_sums(ps_all_lesional_NR)),
                                        rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_nbUVB_les_R.rarefied <- rarefy_even_depth(ps_nbUVB_les_R, sample.size = min(sample_sums(ps_nbUVB_les_R)),
                                    rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_nbUVB_les_NR.rarefied <- rarefy_even_depth(ps_nbUVB_les_NR, sample.size = min(sample_sums(ps_nbUVB_les_NR)),
                                              rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_all_non_R.rarefied<-rarefy_even_depth (ps_all_nonles_R,sample.size = min(sample_sums(ps_all_nonles_R)),
                                    rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_all_non_NR.rarefied <- rarefy_even_depth (ps_all_nonles_NR,sample.size = min(sample_sums(ps_all_nonles_NR)),
                                       rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_all_nbUVB_responders.rarefied <- rarefy_even_depth (ps_all_nbUVB_responders,sample.size = min(sample_sums(ps_all_nbUVB_responders)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_all_nbUVB_nonresponders.rarefied <- rarefy_even_depth (ps_all_nbUVB_nonresponders,sample.size = min(sample_sums(ps_all_nbUVB_nonresponders)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_all_worsened.rarefied <- rarefy_even_depth (ps_all_worsened,sample.size = min(sample_sums(ps_all_worsened)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_worsened_nbUVB.rarefied <-rarefy_even_depth (ps_worsened_nbUVB,sample.size = min(sample_sums(ps_worsened_nbUVB)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_worsened_no_nbUVB.rarefied <-rarefy_even_depth (ps_worsened_no_nbUVB,sample.size = min(sample_sums(ps_worsened_no_nbUVB)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


ps_v1_worsened.rarefied<- rarefy_even_depth (ps_v1_worsened,sample.size = min(sample_sums(ps_v1_worsened)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_v2_worsened.rarefied <- rarefy_even_depth (ps_v2_worsened,sample.size = min(sample_sums(ps_v2_worsened)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_NT_les.rarefied<- rarefy_even_depth (ps_NT_les,sample.size = min(sample_sums(ps_NT_les)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_NT_nonles.rarefied <- rarefy_even_depth (ps_NT_nonles,sample.size = min(sample_sums(ps_NT_nonles)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

##v1 and v2 ; lesional/nonlesional


#Visit 1
set.seed(500)
ps_v1.rarefied <- rarefy_even_depth(v1.clean, sample.size = min(sample_sums(v1.clean)),
                                    rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_v1.nbUVB.rarefied <-rarefy_even_depth(v1.clean_nbUVB, sample.size = min(sample_sums(v1.clean_nbUVB)),
                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_v1.No_nbUVB<- rarefy_even_depth(v1.clean_No_nbUVB, sample.size = min(sample_sums(v1.clean_No_nbUVB)),
                                       rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_v1_les.nbUVB.rarefied <-rarefy_even_depth(ps_v1_lesional_nbUVB, sample.size = min(sample_sums(ps_v1_lesional_nbUVB)),
                                     rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_v1_non.nbUVB.rarefied <-rarefy_even_depth(ps_v1_non_nbUVB, sample.size = min(sample_sums(ps_v1_non_nbUVB)),
                                             rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_v1_les_No_nbUVB.rarefied <- rarefy_even_depth(ps_v1_lesional_No_nbUVB, sample.size = min(sample_sums(ps_v1_lesional_No_nbUVB)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_v1_non.rarefied <-rarefy_even_depth(ps_v1_samples.nonles, sample.size = min(sample_sums(ps_v1_samples.nonles )),
                                       rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#Visit 2
set.seed(500)
ps_v2.rarefied <- rarefy_even_depth(v2.clean, sample.size = min(sample_sums(v2.clean)),
                                    rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_v2.nbUVB.rarefied <- rarefy_even_depth(v2.clean_nbUVB, sample.size = min(sample_sums(v2.clean_nbUVB)),
                                         rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_v2.No_nbUVB.rarefied<- rarefy_even_depth(v2.clean_No_nbUVB, sample.size = min(sample_sums(v2.clean_No_nbUVB)),
                                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

##lesional vs nonlesional
ps_L.rarefied <- rarefy_even_depth(ps_v1_samples.les, sample.size = min(sample_sums(ps_v1_samples.les)),
                                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_NL.rarefied <- rarefy_even_depth(ps_v1_samples.nonles, sample.size = min(sample_sums(ps_v1_samples.nonles)),
                                    rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#lesional/nonlesional with nbUVB
ps_v2_lesional_nbUVB.rarefied <- rarefy_even_depth(ps_v2_lesional_nbUVB, sample.size = min(sample_sums(ps_v2_lesional_nbUVB)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_v2_non_nbUVB.rarefied <- rarefy_even_depth(ps_v2_non_nbUVB, sample.size = min(sample_sums(ps_v2_non_nbUVB)),
                                             rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#Responders/non responders 
 
ps_L.R.rarefied <- rarefy_even_depth(ps_all_lesional_R, sample.size = min(sample_sums(ps_all_lesional_R)),
                                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  
ps_L.NR.rarefied <- rarefy_even_depth(ps_all_lesional_NR, sample.size = min(sample_sums(ps_all_lesional_NR)),
                                      rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#Visit 1 
ps_v1L.rarefied <- rarefy_even_depth(ps_v1_samples.les, sample.size = min(sample_sums(ps_v1_samples.les)),
                                    rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_v1NL.rarefied <- rarefy_even_depth(ps_v1_samples.nonles, sample.size = min(sample_sums(ps_v1_samples.nonles)),
                                     rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_v1R.rarefied <- rarefy_even_depth(ps_v1_R, sample.size = min(sample_sums(ps_v1_R)),
                                     rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_v1NR.rarefied <- rarefy_even_depth(ps_v1_NR, sample.size = min(sample_sums(ps_v1_NR)),
                                      rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_v1LR.rarefied <-  rarefy_even_depth(ps_v1_L_R, sample.size = min(sample_sums(ps_v1_L_R)),
                                       rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_v1LNR.rarefied <-  rarefy_even_depth(ps_v1_L_NR, sample.size = min(sample_sums(ps_v1_L_NR)),
                                       rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
##V2
ps_v2L.rarefied <- rarefy_even_depth(ps_v2_samples.les, sample.size = min(sample_sums(ps_v2_samples.les)),
                                      rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_v2NL.rarefied <- rarefy_even_depth(ps_v2_samples.nonles, sample.size = min(sample_sums(ps_v2_samples.nonles)),
                                     rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
##
R_vs_NT_NR.rarefied <- rarefy_even_depth(R_vs_NT_NR, sample.size = min(sample_sums(R_vs_NT_NR)),
                                         rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
##
ps_NT_vs_R.rarefied <- rarefy_even_depth(ps_NT_vs_R, sample.size = min(sample_sums(ps_NT_vs_R)),
                  rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_NT_vs_NR.rarefied <- rarefy_even_depth(ps_NT_vs_NR, sample.size = min(sample_sums(ps_NT_vs_NR)),
                                         rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps_2NT_vs_R.rarefied <- rarefy_even_depth(ps_2NT_vs_R, sample.size = min(sample_sums(ps_2NT_vs_R)),
                                         rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps_2NT_vs_NR.rarefied <- rarefy_even_depth(ps_2NT_vs_NR, sample.size = min(sample_sums(ps_2NT_vs_NR)),
                                          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)



####-------- Alpha Diversity analysis -------####

#alpha diversity measures
library("phyloseq")
library("ggplot2")
library("vegan")
library("DESeq2")
library ("ggpubr")
##all v1 vs v2
plot_richness(ps_all_rarefied, x="sample", color="Variables", measures= c ("Shannon","Observed"),title = "Visit 1 vs Visit 2 Alpha Diversity")

plot_richness(ps_all_rarefied, x="sample", color="Visit", measures= "Shannon")

plot_richness(ps_all_nbUVB.rarefied, x="sample", measures=c("Observed", "Shannon")) + geom_boxplot(alpha=0.07)


GP <- prune_taxa(taxa_sums(ps_all_rarefied) > 0, ps_all_rarefied)
plot_richness(GP)


plot_richness(ps_all_rarefied, x="sample", measures= "Shannon", title = "Shannon Alpha Diversity Visit 1 vs Visit 2", color = "Variables")+geom_boxplot(alpha=0.06)+theme(legend.position="none") +facet_grid (~Variables) +theme(strip.text = element_text(size=8, color="black",lineheight=5.0,angle=90 ,hjust=1,vjust=1))

plot_richness(ps_all_rarefied, x="sample", measures="Shannon", title = "Shannon Alpha Diversity Across All groups", color = "Variables")


plot_richness(ps_all_less_R.rarefied,x="sample",color="Visit", measures= c("Shannon","Observed"), title= "Alpha Diversity Visit 1 vs Visit 2 All lesional responders")

plot_richness(ps_all_les_NR.rarefied, x="sample",color="Visit", measures= c("Shannon","Observed"), title= "Alpha Diversity V1 vs V2 All lesional Non-responders")

#nbUVB
plot_richness(ps_all_nbUVB.rarefied, measures= c("Observed", "Shannon"), color = "Visit",title = "Alpha Diversity nbUVB-Visit 1 vs Visit 2" )
#no_nbUVB
plot_richness(ps_all_No_nbUVB.rarefied, measures= c("Observed", "Shannon"), color = "Visit",title = "Alpha Diversity No_nbUVB-Visit 1 vs Visit 2" )

##Visit 1 nbUVB
plot_richness(ps_v1.nbUVB.rarefied, x="sample", measures= c("Observed", "Shannon"), title= "Alpha Diversity nbUVB responders vs non-responders at Visit 1",  color = "Treatment_Response")+
  geom_boxplot(alpha=0.6)+
  theme(legend.position="none", axis.text.x=element_text(angle=50,hjust=1,vjust=1,size=10))

##V1 Lesional - Status 
plot_richness(ps_v1_les.nbUVB.rarefied, x="Response", measures= c ("Shannon","Observed","Simpson"), title = "Visit 1-Lesional", color = "Treatment_Response")+
  geom_boxplot(alpha=0.6)+
  theme(legend.position="none", axis.text.x=element_text(angle=50,hjust=1,vjust=1,size=10))+ stat_summary(fun.data = function(x) data.frame(y=2, label = paste("Mean=",mean(x))), geom="text") + theme(legend.position="none")

plot_richness(ps_v1_les.nbUVB.rarefied, x="Response", measures= "Shannon", title = "Visit 1-Lesional nbUVB Samples: Responders vs Non Responders")+geom_boxplot(width = 0.08)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=", round(mean(x), digits = 2))), geom="text") + theme(legend.position="none")

plot_richness(ps_v1_non.nbUVB.rarefied, x="Response", measures= "Shannon", title = "Visit 1-Non Lesional nbUVB Samples: Responders vs Non Responders")+
  geom_boxplot(width = 0.08)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=", round(mean(x), digits = 2))), geom="text") + theme(legend.position="none")

plot_richness(ps_all_nbUVB.rarefied, x="Visit", measures= "Observed", title = "All nbUVB samples: Visit 1 vs Visit 2")+ theme(plot.title=element_text(hjust=0.5)) + geom_boxplot(width = 0.15)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=8))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.5)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=",mean(x)), round.off), geom="text") + theme(legend.position="none")

plot_richness(R_vs_NT_NR.rarefied, x="Treatment_Response", measures= "Shannon", title = "Visit 2 Lesional: Responders vs Not Treated-Non Responders")+
  geom_boxplot(fill='#A4A4A4', color="black",alpha= 1, width = 0.08)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=",mean(x))), geom="text") + theme(legend.position="none")


plot_richness(ps_v2L.rarefied , x="Treatment_Response", measures= "Shannon", title = "Visit 2-Lesional: Response")+geom_boxplot(fill='#A4A4A4', color="black",alpha= 1, width = 0.3)+theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=",mean(x))), geom="text") + theme(legend.position="none")


plot_richness(ps_v2_lesional_nbUVB.rarefied, x="Response", measures= "Shannon", title = "Visit 2-Lesional: Responders vs Non Responders")+
  geom_boxplot(fill='#A4A4A4', color="black",alpha= 0.6, width = 0.1)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=",mean(x))), geom="text") + theme(legend.position="none")

plot_richness(ps_all_non_R.rarefied, x="Visit", measures= "Shannon", title = "Non Lesional Responders V1 vs V2")+
  geom_boxplot(fill='#A4A4A4', color="black",alpha= 0.6, width = 0.1)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=",mean(x))), geom="text") + theme(legend.position="none")



##V1 nonlesional - status
plot_richness(ps_v1NL.rarefied, x="Status_From_V1_V2", measures= c("Shannon","Observed","Simpson"), title = "Visit 1-NonLesional", color = "Status_From_V1_V2")+
  geom_boxplot(alpha=0.6)+
  theme(legend.position="none", axis.text.x=element_text(angle=50,hjust=1,vjust=1,size=10))

plot_richness(ps_all_les.rarefied, measures= c("Shannon","Observed"), title = "Lesional-all samples V1 vs V2", color = " Visit")

##V2 Lesional -Status
plot_richness(ps_v2L.rarefied, x="Status_From_V1_V2", measures=c("Shannon","Observed", "Simpson"), title = "Visit 2-Lesional", color = "Status_From_V1_V2")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10))

##V2 Non-Lesional - Status 
plot_richness(ps_v2NL.rarefied, x="Status_From_V1_V2", measures=c("Shannon","Observed", "Simpson"), title = "Visit 2-NonLesional", color = "Status_From_V1_V2")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10))

####------Abundance plot----####

## V1-V2
top.P <- names(sort(taxa_sums(ps_all_rarefied), decreasing=TRUE))[1:100]
ps.topP <- transform_sample_counts(ps_all_rarefied, function(OTU) OTU/sum(OTU))
ps.topP <- prune_taxa(top.P, ps.topP)
new_labels <- c("lesional" = "Tumor", "non_lesional" = "Normal")
plot_bar(ps.topP, x="Sample", fill="family")+geom_bar(position = "stack", stat="identity") +facet_wrap(Visit+lesional_nonlesional~Treatment_Response, scales = "free_x", labeller = labeller(Visit = label_both,  lesional_nonlesional=new_labels)) + theme(strip.text = element_text(size=10, face="bold", color="black",lineheight=5.0,angle=0 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black")) 

plot_bar(ps.topP, x="Sample", fill="family")+geom_bar(position = "stack", stat="identity") +facet_grid(~ Variables, scales = "free")
##
plot_bar(ps.topP, x="Sample", fill="genus")+geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()+ labs(x = "Variab")


plot_bar(ps.topP, x="Sample", fill="family")+geom_bar(position = "stack", stat="identity") +facet_grid(~Variables, scales = "free",switch = "both", labeller= label_both)+theme(strip.text = element_text(size=10,lineheight=5.0,angle=0 ,hjust=1,vjust=1))


plot_bar(ps.topP, x="Sample", fill="family")+geom_bar(position = "stack", stat="identity") +facet_grid(Visit~Treatment_Response+lesional_nonlesional, switch = "both", scales = "free", labeller = labeller(Visit = label_both, lesional_nonlesional=label_context)) + theme(strip.text = element_text(size=10, face="italic", color="darkblue", lineheight=5.0,angle=0 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black"))


### Visit 1
top.V1 <- names(sort(taxa_sums(ps_v1.rarefied), decreasing=TRUE))[1:100]
ps.topV1 <- transform_sample_counts(ps_v1.rarefied, function(OTU) OTU/sum(OTU))
ps.topV1 <- prune_taxa(top.V1, ps.topV1)
plot_bar(ps.topV1, x="Sample", fill="family") + geom_bar(stat="identity")

##to see lesional and non lesional differences
plot_bar(ps.topV1, x="Sample", fill="family") + geom_bar(stat="identity")+ facet_wrap(~lesional_nonlesional)

## V1 Responders
top.R1 <- names(sort(taxa_sums(ps_v1R.rarefied), decreasing=TRUE))[1:100]
ps.topR1 <- transform_sample_counts(ps_v1R.rarefied, function(OTU) OTU/sum(OTU))
ps.topR1 <- prune_taxa(top.R1, ps.topR1)
plot_bar(ps.topR1, x="Sample", fill="family", title= "Visit 1_Responders") + geom_bar(stat="identity")

##V1 NonResponders
top.NR1 <- names(sort(taxa_sums(ps_v1NR.rarefied), decreasing=TRUE))[1:100]
ps.topNR1 <- transform_sample_counts(ps_v1NR.rarefied, function(OTU) OTU/sum(OTU))
ps.topNR1 <- prune_taxa(top.NR1, ps.topNR1)
plot_bar(ps.topNR1, x="Sample", fill="family", title= "Visit 1_Non_Responders") + geom_bar(stat="identity")

##V1 lesional-responders
V1_R <- names(sort(taxa_sums(ps_v1R.rarefied), decreasing=TRUE))[1:100]
ps.topPv1R <- transform_sample_counts(ps_v1R.rarefied, function(OTU) OTU/sum(OTU))
ps.topPv1R <- prune_taxa(V1_R, ps.topPv1R)
plot_bar(ps.topPv1R, x="Sample", fill="family", title= "Visit 1_Lesional_Responders") + geom_bar(stat="identity")

## V1- Lesional-NonResponders 

V1_NR <- names(sort(taxa_sums(ps_v1LNR.rarefied), decreasing=TRUE))[1:100]
ps.topPv1NR <- transform_sample_counts(ps_v1LNR.rarefied, function(OTU) OTU/sum(OTU))
ps.topPv1NR <- prune_taxa(V1_NR, ps.topPv1NR)
plot_bar(ps.topPv1NR, x="Sample", fill="family", title= "Visit 1_Non_Responders") + geom_bar(stat="identity")


## V2- lesional
top.P2 <- names(sort(taxa_sums(ps_v2L.rarefied), decreasing=TRUE))[1:100]
ps.topP2 <- transform_sample_counts(ps_v2L.rarefied, function(OTU) OTU/sum(OTU))
ps.topP2 <- prune_taxa(top.P2, ps.topP2)
plot_bar(ps.topP2, x="Sample", fill="family", title= "Visit 2- lesional") + geom_bar(stat="identity")

#v2 lesional Responders
V2_LR <- names(sort(taxa_sums(ps_v2LR.rarefied), decreasing=TRUE))[1:100]
ps.topV2R <- transform_sample_counts(ps_v2LR.rarefied, function(OTU) OTU/sum(OTU))
ps.topV2R <- prune_taxa(V2_LR, ps.topV2R)
plot_bar(ps.topV2R, x="Sample", fill="family", title= "Visit 2_Lesional_Responders") + geom_bar(stat="identity")


## V2- Lesional-NonResponders 

V2_NR <- names(sort(taxa_sums(ps_v2LNR.rarefied), decreasing=TRUE))[1:100]
ps.topPv2NR <- transform_sample_counts(ps_v2LNR.rarefied, function(OTU) OTU/sum(OTU))
ps.topPv2NR <- prune_taxa(V2_NR, ps.topPv2NR)
plot_bar(ps.topPv2NR, x="Sample", fill="family", title= "Visit 2_Lesional_Non_Responders") + geom_bar(stat="identity")

## V2 Non-lesional 
V2_NL <- names(sort(taxa_sums(ps_v2NL.rarefied), decreasing=TRUE))[1:100]
ps.top2NL <- transform_sample_counts(ps_v2NL.rarefied, function(OTU) OTU/sum(OTU))
ps.top2NL <- prune_taxa(V2_NL, ps.top2NL)
plot_bar(ps.top2NL, x="Sample", fill="family", title= "Visit 2- Nonlesional") + geom_bar(stat="identity")


####----Estimate richness ------####

##All samples
rich_all = estimate_richness(ps_all_rarefied)


rich_all_nbUVB<-estimate_richness(ps_all_nbUVB.rarefied)
rich_all_No_nbUVB<-estimate_richness(ps_all_No_nbUVB.rarefied)

rich_all_les=estimate_richness(ps_all_les.rarefied)
rich_all_nonles = estimate_richness(ps_all_non.rarefied)

rich_all_les_NR=estimate_richness(ps_all_les_NR.rarefied)
rich_all_les_R = estimate_richness(ps_all_less_R.rarefied)

rich_all_nbUVB_nonresponders=estimate_richness(ps_all_nbUVB_nonresponders.rarefied)
rich_all_nbUVB_responders=estimate_richness(ps_all_nbUVB_responders.rarefied)

rich_nbUVB_les_R = estimate_richness(ps_nbUVB_les_R.rarefied)
rich_nbUVB_les_NR = estimate_richness(ps_nbUVB_les_NR.rarefied)

rich_NT_les =estimate_richness(ps_NT_les.rarefied)
rich_NT_nonles =estimate_richness(ps_NT_nonles.rarefied)

rich_all_non_R =estimate_richness(ps_all_non_R.rarefied)
rich_all_non_NR=estimate_richness(ps_all_non_NR.rarefied)
##V1
rich_v1 = estimate_richness (ps_v1.rarefied)
rich_v1.nbUVB <-estimate_richness(ps_v1.nbUVB.rarefied)

rich_v1_les= estimate_richness (ps_v1L.rarefied)
rich_v1_nbUVB_les = estimate_richness(ps_v1_les.nbUVB.rarefied)
rich_v1_nbUVB_non = estimate_richness(ps_v1_non.nbUVB.rarefied)
rich_v1_No_nbUVB_les = estimate_richness(ps_v1_les_No_nbUVB.rarefied)

rich_v1_nonles =estimate_richness(ps_v1NL.rarefied)
rich_v1_les_R = estimate_richness(ps_v1LR.rarefied)
rich_v1_les_NR = estimate_richness(ps_v1LNR.rarefied)
rich_v1_R =estimate_richness(ps_v1R.rarefied)
rich_v1_NR= estimate_richness(ps_v1NR.rarefied)

##V2
rich_v2 = estimate_richness(ps_v2.rarefied)
rich_v2_nbUVB= estimate_richness(ps_v2.nbUVB.rarefied)
rich_v2_les = estimate_richness(ps_v2L.rarefied)
rich_v2_nbUVB_les= estimate_richness(ps_v2_lesional_nbUVB.rarefied)
rich_v2_nonles = estimate_richness(ps_v2NL.rarefied)
rich_v2_non_nbUVB <-estimate_richness(ps_v2_non_nbUVB.rarefied)

rich_v2_R_les = estimate_richness(ps_v2LR.rarefied)
rich_v2_NR_les = estimate_richness(ps_v2LNR.rarefied)
rich_v2_R_nonles = estimate_richness(ps_v2NonR.rarefied)
rich_v2_NR_nonles= estimate_richness(ps_v2NonNR.rarefied)

rich_all_worsened =estimate_richness(ps_all_worsened.rarefied)
rich_no_nbUVB_worsened =estimate_richness(ps_worsened_no_nbUVB.rarefied)
rich_nbUVB_worsened =estimate_richness(ps_worsened_nbUVB.rarefied)
rich_v1_worsened =estimate_richness(ps_v1_worsened.rarefied)
rich_v2_worsened = estimate_richness(ps_v2_worsened.rarefied)

rich_NT_vs_R = estimate_richness(ps_NT_vs_R.rarefied)
rich_NT_vs_NR = estimate_richness(ps_NT_vs_NR.rarefied)
rich_2NT_vs_R = estimate_richness(ps_2NT_vs_R.rarefied)
rich_2NT_vs_NR = estimate_richness(ps_2NT_vs_NR.rarefied)

####----------Pairwise Wilcox test------#####
#all samples v1 vs v2

#0.011 Sp value
pairwise.wilcox.test(rich_all_nbUVB$Simpson, sample_data(ps_all_nbUVB.rarefied)$Visit, p.adjust.method= "none", paired = TRUE)
#0.08 for all samples 
pairwise.wilcox.test(rich_all$Simpson,sample_data(ps_all_rarefied)$Visit,p.adjust.method= "none", paired = TRUE)


#0.52 
pairwise.wilcox.test(rich_all_No_nbUVB$Simpson, sample_data(ps_all_No_nbUVB.rarefied)$Visit,  p.adjust.method= "none", paired = TRUE)
##all lesional/nonlesional
#0.8
pairwise.wilcox.test(rich_all_les$Shannon, sample_data(ps_all_les.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)

#Non/responders lesional
pairwise.wilcox.test(rich_all_les_NR$Shannon, sample_data(ps_all_les_NR.rarefied)$Visit,p.adjust.method= "none", paired = TRUE)
#0.57
pairwise.wilcox.test(rich_all_les_R$Observed, sample_data(ps_all_less_R.rarefied)$Visit, p.adjust.method= "none",paired = TRUE, exact=FALSE)
#0.87
pairwise.wilcox.test(rich_all_les_NR$Observed, sample_data(ps_all_les_NR.rarefied)$Visit, p.adjust.method= "none",paired = TRUE, exact=FALSE)
#0.44


pairwise.wilcox.test(rich_all_les_NR$Observed, sample_data(ps_all_les_NR.rarefied)$Visit, p.adjust.method= "none",paired = TRUE, exact=FALSE)



#0.16
pairwise.wilcox.test(rich_all_non_R$Observed,sample_data(ps_all_non_R.rarefied)$Visit, p.adjust.method= "none",paired = TRUE, exact=FALSE)
pairwise.wilcox.test(rich_all_non_NR$Simpson,sample_data(ps_all_non_NR.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)
# just nbUVB -lesion R v1v2 0.31
pairwise.wilcox.test(rich_nbUVB_les_R$Chao1, sample_data(ps_nbUVB_les_R.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)

#
pairwise.wilcox.test(rich_nbUVB_les_NR$Simpson, sample_data(ps_nbUVB_les_NR.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)


#not treated lesional, nonlesional v1v2 0.82 and 0.48
pairwise.wilcox.test(rich_NT_les$Shannon, sample_data(ps_NT_les.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)
pairwise.wilcox.test(rich_NT_nonles$Shannon, sample_data(ps_NT_nonles.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)


#V1

#0.13
pairwise.wilcox.test(rich_v1.nbUVB$Simpson, sample_data(ps_v1.nbUVB.rarefied)$Response,p.adjust.method= "none", paired= FALSE)
#.077

pairwise.wilcox.test(rich_v1_nbUVB_les$Simpson, sample_data(ps_v1_les.nbUVB.rarefied)$Response,p.adjust.method= "none",paired= FALSE)#0.78
#0.47
pairwise.wilcox.test(rich_v1_nbUVB_non$Shannon, sample_data(ps_v1_non.nbUVB.rarefied)$Response,p.adjust.method= "none",paired= FALSE)
#0.022
pairwise.wilcox.test(rich_v2_nbUVB_les$Simpson, sample_data(ps_v2_lesional_nbUVB.rarefied)$Response,p.adjust.method= "none",paired= FALSE)
#0.036
pairwise.wilcox.test(rich_v2_non_nbUVB$Simpson, sample_data(ps_v2_non_nbUVB.rarefied)$Response,p.adjust.method= "none", paired= FALSE)

#across all samples V1 vs V2 = p 0.44
pairwise.wilcox.test(rich_all$Shannon, sample_data(ps_all_rarefied)$Visit, p.adjust.method= "none", paired = TRUE)

pairwise.wilcox.test(rich_all_nbUVB$Shannon, sample_data(ps_all_nbUVB.rarefied)$Visit, p.adjust.method= "none", paired = TRUE)

pairwise.wilcox.test(rich_all_nbUVB$Shannon, sample_data(ps_all_nbUVB.rarefied)$Response, p.adjust.method= "none", paired = FALSE) #0.028
pairwise.wilcox.test(rich_v1_nbUVB_les$Shannon, sample_data (ps_v1_les.nbUVB.rarefied)$Response,p.adjust.method= "none", paired= FALSE) #0.53

pairwise.wilcox.test(rich_all_les$Simpson, sample_data(ps_L.R.rarefied)$Visit, p.adjust.method= "none", paired = TRUE) 
#p value 0.19

#0.47
pairwise.wilcox.test(rich_v1_nbUVB_non$Shannon, sample_data(ps_v1_non.nbUVB.rarefied)$Response, p.adjust.method= "none", paired = FALSE) 

## all Visit 1- R and NR #0.97
pairwise.wilcox.test(rich_v1$Shannon, sample_data(ps_v1.rarefied)$Treatment_Response,p.adjust.method= "none", paired = FALSE, exact= FALSE)

#V1-lesional (R vs NR) #0.33 #0.57 for just nbUVB 
pairwise.wilcox.test(rich_v1_les$Simpson, sample_data(ps_v1L.rarefied)$Treatment_Response, p.adjust.method= "none",paired = FALSE)


wilcox.test (rich_v1_nbUVB_les$Shannon, rich_v1_No_nbUVB_les$Shannon, paired = FALSE) #0.24

#V2 lesional Responders / Non Responders

##R/Les post-UV richness > NR/Les post-UV richness (Shannon p=0.016) #we get p value = 0.037
pairwise.wilcox.test(rich_v2_les$Observed, sample_data(ps_v2L.rarefied)$Treatment_Response, paired = FALSE, p.adjust.method= "none", exact= FALSE)

##R/Non post-UV richness no different from NR/Non post-UV richness (Shannon p=0.497)
pairwise.wilcox.test(rich_v2_nonles$Shannon, sample_data(ps_v2NL.rarefied)$Treatment_Response, p.adjust.method= "none", paired = FALSE) # p value-0.026
##

pairwise.wilcox.test(rich_v2_nbUVB$Simpson, sample_data(ps_v2.nbUVB.rarefied)$Response, p.adjust.method= "none", paired = FALSE)
wilcox.test(rich_v2_R_les$Shannon, rich_v2_NR_les$Shannon, paired = FALSE) #0.02119

pairwise.wilcox.test(rich_v2_nonles$Shannon, sample_data(ps_v2_samples.nonles)$Response, p.adjust.method= "none", paired= FALSE)  #0.026

wilcox.test(rich_v2_les$Shannon, rich_v2_nonles$Shannon, paired = TRUE) #0.0048
 

pairwise.wilcox.test(rich_v1_R$Simpson, sample_data(ps_v1R.rarefied)$Visit, paired = FALSE)

pairwise.wilcox.test(rich_no_nbUVB_worsened$Shannon, sample_data(ps_worsened_no_nbUVB.rarefied)$Visit, paired = TRUE)


###nbUVB vs no nbUVB

#0.56
pairwise.wilcox.test(rich_v1$Shannon, sample_data(ps_v1.rarefied)$Treatment,p.adjust.method= "none", paired = FALSE)
#0.06
pairwise.wilcox.test(rich_v2$Simpson, sample_data(ps_v2.rarefied)$Treatment_Response,p.adjust.method= "none", paired = FALSE, exact= FALSE)
#0.79
pairwise.wilcox.test(rich_all_worsened$Shannon, sample_data(ps_all_worsened.rarefied)$Treatment, p.adjust.method = "none", paired = FALSE)
#0.69
pairwise.wilcox.test(rich_v1_worsened$Shannon, sample_data(ps_v1_worsened.rarefied)$Treatment, p.adjust.method = "none", paired = FALSE)
#0.78
pairwise.wilcox.test(rich_v2_worsened$Shannon, sample_data(ps_v2_worsened.rarefied)$Treatment, p.adjust.method = "none", paired = FALSE)


pairwise.wilcox.test(rich_all_nbUVB_responders$Shannon, sample_data(ps_all_nbUVB_responders.rarefied)$Visit, p.adjust.method = "none", paired = TRUE)
pairwise.wilcox.test(rich_all_nbUVB_nonresponders$Shannon, sample_data(ps_all_nbUVB_nonresponders.rarefied)$Visit, p.adjust.method = "none", paired = TRUE)

pairwise.wilcox.test(rich_all_No_nbUVB$Shannon, sample_data(ps_all_No_nbUVB.rarefied)$Visit, p.adjust.method = "none", paired = TRUE)

rich_all_No_nbUVB

 ####--------- Beta Diversity --------####

# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(ps_all_rarefied, method="bray", weighted=F)
ordination = ordinate(ps_all_rarefied, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps_all_rarefied, ordination, color="Visit", title= "Beta Diversity-All samples Visit 1 vs Visit 2") + theme(aspect.ratio=1)


#v1-: responders vs non responders

distv1 = phyloseq::distance(ps_all_nbUVB.rarefied, method="bray")
ordination.v1 = ordinate(ps_all_nbUVB.rarefied, method="PCoA", distance=distv1)
plot_ordination(ps_all_nbUVB.rarefied, ordination.v1, color="Visit", title= "Beta Diversity-nbUVB Visit 1 vs Visit 2") + 
  theme_classic() +
  theme(strip.background = element_blank())+stat_ellipse(aes(group = Visit), linetype = 7)
#
dist= phyloseq::distance(ps_all_rarefied, method="bray", weighted = T)
ordination = ordinate(ps_all_rarefied, method="PCoA", distance=dist)
plot_ordination(ps_all_rarefied,ordination,color="Visit", title= "Beta Diversity-All samples Visit 1 vs Visit 2") +
  theme_classic() +
  theme(strip.background = element_blank())+stat_ellipse(aes(group = Visit), linetype = 1)
#not treated 
distNT = phyloseq::distance(ps_all_No_nbUVB.rarefied, method="bray")
ordination=ordinate(ps_all_No_nbUVB.rarefied, method="PCoA", distance=distNT)
plot_ordination(ps_all_No_nbUVB.rarefied,ordination,color="Visit", title= "Beta Diversity-All Not Treated samples Visit 1 vs Visit 2") +
  theme_classic() +
  theme(strip.background = element_blank())+stat_ellipse(aes(group = Visit), linetype = 2)

distNTles = phyloseq::distance(ps_NT_les.rarefied, method="bray")
ordination=ordinate(ps_NT_les.rarefied, method="PCoA", distance=distNTles)
plot_ordination(ps_NT_les.rarefied,ordination,color="Visit", title= "Beta Diversity-Not Treated lesional samples Visit 1 vs Visit 2") +
  theme_classic() +
  theme(strip.background = element_blank())

distNTnon = phyloseq::distance(ps_NT_nonles.rarefied, method="bray")
ordination=ordinate(ps_NT_nonles.rarefied, method="PCoA", distance=distNTnon)
plot_ordination(ps_NT_nonles.rarefied,ordination,color="Visit", title= "Beta Diversity-Not Treated lesional samples Visit 1 vs Visit 2") +
  theme_classic() +
  theme(strip.background = element_blank())


distR= phyloseq::distance(ps_all_nbUVB_responders.rarefied, method="bray")
ordinationR = ordinate(ps_all_nbUVB_responders.rarefied, method="PCoA", distance=distR)
plot_ordination(ps_all_nbUVB_responders.rarefied,ordinationR,color="Visit", title= "Beta Diversity-All responders Visit 1 vs Visit 2") +
  theme_classic() +
  theme(strip.background = element_blank())



distNR= phyloseq::distance(ps_all_nbUVB_nonresponders.rarefied, method="bray")
ordinationNR = ordinate(ps_all_nbUVB_nonresponders.rarefied, method="PCoA", distance=distNR)
plot_ordination(ps_all_nbUVB_nonresponders.rarefied,ordinationNR,color="Visit", title= "Beta Diversity-All nonresponders Visit 1 vs Visit 2") +
  theme_classic() +
  theme(strip.background = element_blank())

##worsened 
distW1= phyloseq::distance(ps_v2_worsened.rarefied, method="bray")
ordinationW1 = ordinate(ps_v2_worsened.rarefied, method="PCoA", distance=distW1)
plot_ordination(ps_v2_worsened.rarefied, ordinationW1,color="Treatment", title= "Beta Diversity-All Worsened at Visit 2") +
  theme_classic() +
  theme(strip.background = element_blank())


#v1- lesional: responders vs non responders

distv1L.t = phyloseq::distance(ps_L.rarefied, method="bray")
ordination.v1L.t = ordinate(ps_L.rarefied, method="PCoA", distance=distv1L.t)
plot_ordination(ps_L.rarefied, ordination.v1L.t, color="Treatment_Response", title= "Visit 1-All Samples Lesional: Treatment Response ") + 
  theme_classic() +
  theme(strip.background = element_blank())
##	R/Les vs NR/Les clustering after UV (p=0.045) but not before (p=0.181)
distv1L = phyloseq::distance(ps_v1_les.nbUVB.rarefied, method="bray")
ordination.v1L = ordinate(ps_v1_les.nbUVB.rarefied, method="PCoA", distance=distv1L)
plot_ordination(ps_v1_les.nbUVB.rarefied, ordination.v1L, color="Treatment_Response", title= "Visit 1- nbUVB Lesional Samples: Responders vs Non Responders") + 
  theme_classic() +
  theme(strip.background = element_blank())+stat_ellipse(aes(group = Treatment_Response), linetype = 2)


#V1 Non lesional: Status
distv1N = phyloseq::distance(ps_v1_non.nbUVB.rarefied, method="bray")
ordination.v1N = ordinate(ps_v1_non.nbUVB.rarefied, method="PCoA", distance=distv1N)
plot_ordination(ps_v1_non.nbUVB.rarefied, ordination.v1N, color="Response", title= "Visit 1 nbUVB-Nonlesional") + 
  theme_classic() +
  theme(strip.background = element_blank())

##V1 All Non lesional: Status
distv1N.t = phyloseq::distance(ps_v1_non.rarefied, method="bray")
ordination.v1N.t= ordinate(ps_v1_non.rarefied, method="PCoA", distance=distv1N.t)
plot_ordination(ps_v1_non.rarefied, ordination.v1N.t, color="Treatment_Response", title= "Visit 1-Nonlesional All Samples: Treatment Response") + theme_classic() +
  theme(strip.background = element_blank())  

##V2 

distv2 = phyloseq::distance(ps_v2.nbUVB.rarefied, method="bray")
ordination.v2= ordinate(ps_v2.rarefied, method="PCoA", distance=distv2)
plot_ordination(ps_v2.rarefied, ordination.v2, color="lesional_nonlesional", title= "Visit 2-All samples lesional vs non-lesional") + theme_classic() +theme(strip.background = element_blank())  

##V2 lesional
distv2L.t = phyloseq::distance(ps_v2L.rarefied, method="bray")
ordination.v2L.t = ordinate(ps_v2L.rarefied, method="PCoA", distance=distv2L.t)
plot_ordination(ps_v2L.rarefied, ordination.v2L.t, color="Treatment_Response", title= "Visit 2-All Lesional Samples: Treatment Response") + theme_classic() + theme(strip.background = element_blank())


distv2.les = phyloseq::distance(ps_v2_lesional_nbUVB.rarefied, method="bray")
ordination.v2.les = ordinate(ps_v2_lesional_nbUVB.rarefied, method="PCoA", distance=distv2.les)
plot_ordination(ps_v2_lesional_nbUVB.rarefied, ordination.v2.les, color="Response", title= "Visit 2- nbUVB Lesional Samples : Responders vs Non Responders") + 
  theme_classic() +
  theme(strip.background = element_blank()) 

#V2 Nonlesional
distv2.non.t = phyloseq::distance(ps_v2NL.rarefied, method="bray")
ordination.v2.non.t = ordinate(ps_v2NL.rarefied, method="PCoA", distance=distv2.non.t)
plot_ordination(ps_v2NL.rarefied, ordination.v2.non.t, color="Treatment_Response", title= "Visit 2- All Samples Non-Lesional: Treatment Response") + theme_classic() +theme(strip.background = element_blank()) 

distv2.non = phyloseq::distance(ps_v2_non_nbUVB.rarefied, method="bray")
ordination.v2.non = ordinate(ps_v2_non_nbUVB.rarefied, method="PCoA", distance=distv2.non)
plot_ordination(ps_v2_non_nbUVB.rarefied, ordination.v2.non, color="Response", title= "Visit 2 nbUVB-NonLesional") + 
  theme_classic() +
  theme(strip.background = element_blank())

###------Adonis tests-----####
perm <- how(nperm = 500)
#not treated 
adonis2(distNT~sample_data(ps_all_No_nbUVB.rarefied)$Visit, method = "bray", permutations= perm) #r2=0.00659, p=0.4251 
adonis2(distNTles~sample_data(ps_NT_les.rarefied)$Visit, method = "bray", permutations= perm)#0.0075; 0.8862
adonis2(distNTnon~sample_data(ps_NT_nonles.rarefied)$Visit, method = "bray", permutations= perm)#0.00836; 0.8423

adonis2(distv1L.t~sample_data(ps_v1L.rarefied)$Treatment_Response, method = "bray", permutations= perm) #0.06857; 0.005988
adonis2(distv1L~sample_data(ps_v1_les.nbUVB.rarefied)$Response, method = "bray", permutations= perm) #0.04918; 0.01397

adonis2(distv2L.t~sample_data(ps_v2L.rarefied)$Treatment_Response, method = "bray", permutations= perm) #0.06821;0.003992
vegan::adonis2(distv2.les~sample_data(ps_v2_lesional_nbUVB.rarefied)$Response, method = "bray", permutations= perm)#0.03166; 0.08982

vegan::adonis2(distv1~sample_data(ps_all_nbUVB.rarefied)$Visit, method = "bray", permutations= perm)#0.005; 0.36
vegan::adonis2(distv2~sample_data(ps_v2.nbUVB.rarefied)$Response, method = "bray", permutations= perm) #0.03358; 0.001

adonis2(distv1~sample_data(ps_all_nbUVB_nonresponders.rarefied)$Treatment_Response, method = "bray", permutations= perm)

adonis2(distv~sample_data(ps_v1L.rarefied)$Treatment_Response, method = "bray", permutations= perm)

adonis2(distv1L~sample_data(ps_v1_les.nbUVB.rarefied)$Response, method = "bray", permutations= perm)#0.04918; 0.01796

adonis2(distv1N~sample_data(ps_v1_non.nbUVB.rarefied)$Response, method = "bray", permutations= perm) #0.07685; 0.005988

adonis2(distR~sample_data(ps_all_nbUVB_responders.rarefied)$Visit, method = "bray", permutations= perm) #0.0138; 0.08184
adonis2(distNR~sample_data(ps_all_nbUVB_nonresponders.rarefied)$Visit, method = "bray", permutations= perm)#0.02247; 0.07984 

adonis2(distv2.non~sample_data(ps_v2_non_nbUVB.rarefied,)$Response, method = "bray", permutations= perm) #0.03; 0.045
##nonlesional
perm <- how(nperm = 500) 
adonis2(distv1N.t~sample_data(ps_v1_non.rarefied)$Treatment_Response, method = "bray", permutations= perm) #0.09238; 0.001996
adonis2(distv1N~sample_data(ps_v1_non.nbUVB.rarefied)$Response, method = "bray", permutations= perm) #0.07685; 0.003992

adonis2(distv2.non.t~sample_data(ps_v2NL.rarefied)$Treatment_Response, method = "bray", permutations= perm)

adonis2(distW1~sample_data(ps_all_worsened.rarefied)$Response, method = "bray", permutations= perm)

####---------Skeleton questions------####

#1.)
#across all samples V1 vs V2 = p 0.11 ; 0.079
pairwise.wilcox.test(rich_all$Observed, sample_data(ps_all_rarefied)$Visit, p.adjust.method= "none", paired = TRUE)
##all v1 vs v2
plot_richness(ps_all_rarefied, x="sample", color="Visit", measures= c ("Shannon","Observed"),title = "Alpha Diversity- All samples Visit 1 vs Visit 2")
#2.)

## all Visit 1- R and NR #0.25 #0.89
pairwise.wilcox.test(rich_v1$Observed, sample_data(ps_v1.rarefied)$Treatment,p.adjust.method= "none", paired = FALSE)

plot_richness(ps_v1.rarefied, x="sample", color="Response", measures= c ("Shannon","Observed"),title = "Alpha Diversity- All samples at Visit 1")

pairwise.wilcox.test(rich_v1.nbUVB$Shannon, sample_data(ps_v1.nbUVB.rarefied)$Response,p.adjust.method= "none", paired= FALSE)#0.077
##
#across all nbUVB samples V1 vs V2 p-value= 0.008
pairwise.wilcox.test(rich_all_nbUVB$Shannon, sample_data(ps_all_nbUVB.rarefied)$Visit, p.adjust.method= "none", paired = FALSE)

#3.)#0.0044
pairwise.wilcox.test(rich_v2_les$Shannon, sample_data(ps_v2L.rarefied)$Response, paired = FALSE, p.adjust.method= "none")
#0.022
pairwise.wilcox.test(rich_v2_nbUVB_les$Shannon, sample_data(ps_v2_lesional_nbUVB.rarefied)$Treatment_Response,p.adjust.method= "none",paired= FALSE)

plot_richness(ps_v2_lesional_nbUVB.rarefied, measures = c("Observed", "Shannon"), color= "Response", title = "At Visit 2 Lesional nbUVB responders vs non responders")
plot_richness(ps_v2L.rarefied, measures = c("Observed", "Shannon"), color= "Response", title = "At Visit 2 Lesional all samples")

#4.) #0.021 N vs NR; 0.026 N vs NR nbUVB
pairwise.wilcox.test(rich_v2_nonles$Shannon, sample_data(ps_v2NL.rarefied)$Treatment_Response, p.adjust.method= "none", paired = FALSE)
pairwise.wilcox.test(rich_v2_non_nbUVB$Shannon, sample_data(ps_v2_non_nbUVB.rarefied)$Response,p.adjust.method= "none", paired= FALSE)

plot_richness(ps_v2NL.rarefied, measures = c("Observed", "Shannon"), color= "Response", title = "At Visit 2- nonlesional all samples")
plot_richness(ps_v2_non_nbUVB.rarefied, measures = c("Observed", "Shannon"), color= "Response", title = "At Visit 2- nonlesional nbUVB responders vs non responders")

#5.) n vs NR 0.35
pairwise.wilcox.test(rich_v1_les$Shannon, sample_data(ps_v1L.rarefied)$Response, paired = FALSE, p.adjust.method= "none")


plot_richness(ps_v1L.rarefied, measures= c ("Shannon","Observed"), title = "At Visit 1-Lesional all samples", color = "Response")
plot_richness(ps_v1_les.nbUVB.rarefied, measures = c("Observed", "Shannon"), color= "Response", title = "At Visit 1 Lesional nbUVB responders vs non responders")

#6.) #0.47
pairwise.wilcox.test(rich_v1_nonles$Shannon, sample_data(ps_v1NL.rarefied)$Treatment_Response, p.adjust.method = "none", paired= FALSE)
pairwise.wilcox.test(rich_v1_nbUVB_non$Shannon, sample_data(ps_v1_non.nbUVB.rarefied)$Response, p.adjust.method = "none", paired= FALSE) #0.51

#0.070
pairwise.wilcox.test(rich_v2_nonles$Shannon, sample_data(ps_v2NL.rarefied)$Response, p.adjust.method= "none", paired = FALSE)
#0.026
pairwise.wilcox.test(rich_v2_non_nbUVB$Shannon, sample_data(ps_v2_non_nbUVB.rarefied)$Response,p.adjust.method= "none", paired= FALSE)

#7.) #0.49
pairwise.wilcox.test(rich_all_les$Simpson, sample_data(ps_all_les.rarefied)$Visit,p.adjust.method= "none", paired= TRUE)


plot_richness(ps_v2NL.rarefied, measures = c("Observed", "Shannon"), color= "Response", title = "At Visit 2- nonlesional all samples")
plot_richness(ps_v2_non_nbUVB.rarefied, measures = c("Observed", "Shannon"), color= "Response", title = "At Visit 2- nonlesional nbUVB responders vs non responders")

#8.) #0.22
pairwise.wilcox.test(rich_all_les_R$Shannon, sample_data(ps_all_less_R.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)
plot_richness (ps_all_less_R.rarefied, color= "Visit" , measures = c("Observed", "Shannon"), title= "All lesional Responders")
#.57
pairwise.wilcox.test(rich_all_les_NR$Simpson, sample_data(ps_all_les_NR.rarefied)$Visit,p.adjust.method= "none", paired = TRUE)
plot_richness (ps_all_les_NR.rarefied, color= "Visit" , measures = c("Observed", "Shannon"), title= "All lesional Non-Responders")


distLNR = phyloseq::distance(ps_all_les_NR.rarefied, method="bray")
ordination.LNR = ordinate(ps_all_les_NR.rarefied, method="PCoA", distance=distLNR)
plot_ordination(ps_all_les_NR.rarefied, ordination.LNR, color="Visit", title= "Beta Diversity-All Lesional Nonresponders Visit 1 vs Visit 2") +theme_classic() + theme(strip.background = element_blank())

distLR = phyloseq::distance(ps_all_less_R.rarefied, method="bray")
ordination.LR = ordinate(ps_all_less_R.rarefied, method="PCoA", distance=distLR)
plot_ordination(ps_all_less_R.rarefied, ordination.LR, color="Visit", title= "Beta Diversity-All Lesional responders Visit 1 vs Visit 2") +theme_classic() + theme(strip.background = element_blank())

#9.)
#0.22
pairwise.wilcox.test(rich_all_les_R$Shannon, sample_data(ps_all_less_R.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)
plot_richness (ps_all_less_R.rarefied, color= "Visit" , measures = c("Observed", "Shannon"), title= "All lesional Responders")

#0.57
pairwise.wilcox.test(rich_all_les_NR$Shannon, sample_data(ps_all_les_NR.rarefied)$Visit,p.adjust.method= "none", paired = TRUE) 
plot_richness (ps_all_les_NR.rarefied, color= "Visit" , measures = c("Observed", "Shannon"), title= "All lesional Non-Responders")

#0.016
pairwise.wilcox.test(rich_all_non_R$Shannon, sample_data(ps_all_non_R.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)
plot_richness (ps_all_non_R.rarefied, color= "Visit" , measures = c("Observed", "Shannon"), title= "All lesional Responders")
pairwise.wilcox.test(rich_all_non_NR$Shannon, sample_data(ps_all_non_NR.rarefied)$Visit, p.adjust.method= "none",paired = TRUE)
plot_richness (ps_all_les_NR.rarefied, color= "Visit" , measures = c("Observed", "Shannon"), title= "All lesional Non-Responders")

#additional 
adonis2(distv1L~sample_data(ps_v1_les.nbUVB.rarefied)$Response, method = "bray", permutations= perm)

adonis2(distv1N~sample_data(ps_v1_non.nbUVB.rarefied)$Response, method = "bray", permutations= perm)

##
pairwise.wilcox.test(rich_all_nbUVB_responders$Shannon, sample_data(ps_all_nbUVB_responders.rarefied)$Visit, p.adjust.method = "none", paired = TRUE)
pairwise.wilcox.test(rich_all_nbUVB_nonresponders$Shannon, sample_data(ps_all_nbUVB_nonresponders.rarefied)$Visit, p.adjust.method = "none", paired = TRUE)
pairwise.wilcox.test(rich_all_No_nbUVB$Shannon, sample_data(ps_all_No_nbUVB.rarefied)$Visit, p.adjust.method = "none", paired = TRUE)


#0.58 shannon p
pairwise.wilcox.test(rich_NT_vs_R$Shannon, sample_data(ps_NT_vs_R.rarefied)$Treatment, p.adjust.method = "none", paired = FALSE)
#0.76 shannon
pairwise.wilcox.test(rich_NT_vs_NR$Shannon, sample_data(ps_NT_vs_NR.rarefied)$Treatment, p.adjust.method = "none", paired = FALSE)

##0.0024 Shannon p
pairwise.wilcox.test(rich_2NT_vs_R$Shannon, sample_data(ps_2NT_vs_R.rarefied)$Treatment, p.adjust.method = "none", paired = FALSE)
##0.88 Shannon p
pairwise.wilcox.test(rich_2NT_vs_NR$Shannon, sample_data(ps_2NT_vs_NR.rarefied)$Treatment, p.adjust.method = "none", paired = FALSE)
#####-------Plots---------####
BiocManager::install("HMP16SData", force= TRUE)
BiocManager::install("microshades")
library(HMP16SData)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(microshades)
#Abundance plot
top.P <- names(sort(taxa_sums(ps_all_rarefied), decreasing=TRUE))[1:100]
ps.topP <- transform_sample_counts(ps_all_rarefied, function(OTU) OTU/sum(OTU))
ps.topP <- prune_taxa(top.P, ps.topP)
new_labels <- c("lesional" = "Tumor", "non_lesional" = "Normal")
plot_bar(ps.topP, x="Sample", fill="family")+geom_bar(position = "stack", stat="identity") +facet_grid(~Variables, switch = "both", scales = "free_x", labeller = labeller(Visit = label_both, Response= label_both, Treatment = label_both, lesional_nonlesional=new_labels)) + theme(strip.text = element_text(size=10, face="bold", color="black",lineheight=5.0,angle=90 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black")) 

#
install.packages("remotes")
remotes::install_github("vmikk/metagMisc", force=TRUE)
library ("metagMisc")
#
top.P <- names(sort(taxa_sums(ps_all_rarefied), decreasing=TRUE))[1:100]
ps.topP <- phyloseq_standardize_otu_abundance(ps_all_rarefied, method= "total")
ps.topP <- prune_taxa(top.P, ps.topP)
new_labels <- c("lesional" = "Tumor", "non_lesional" = "Normal")
plot_bar(ps.topP, x="Sample", fill="genus")+geom_bar(position = "stack", stat="identity") +facet_grid(~Variables, scales = "free")
## do 

#Alpha Diversity
install.packages("ggpubr")
library ("ggpubr")
plot_richness(ps_all_rarefied, x="sample", color="Variables", measures= "Shannon") +facet_grid (~Variables,switch = "both", scales = "free")+geom_boxplot(alpha= 0.6)+theme(legend.position="none", axis.text.x=element_text(angle=0 ,hjust=1,vjust=1,size=10))+theme(strip.text = element_text(size=8, color="black",lineheight=5.0,angle=0 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black")) 

plot_richness(ps_v1_les.nbUVB.rarefied, x="Treatment_Response", measures=c("Shannon", "Chao1", "Observed"), title= "Alpha Diversity",  color = "Treatment_Response")+
  geom_boxplot(alpha=3)+
  theme(legend.position="none", axis.text.x=element_text(angle=0
                                                         ,hjust=1,vjust=1,size=9))


plot_richness(ps_v1_les.nbUVB.rarefied, x="Response", measures=c("Chao1", "Observed", "Shannon"), title= "Alpha Diversity",  color = "Response")+ geom_boxplot(alpha=0.15)+
  theme(legend.position="none", axis.text.x=element_text(angle=0)) 

                                                    




install.packages("remotes")
remotes::install_github("vallenderlab/MicrobiomeR", force=  TRUE)

install.packages ("devtools")
install.packages("MicrobiomeR")
library (devtools)
devtools::install_github("wsteenhu/microbiomer", force= TRUE)
library ("microbiomer")
library(microbiomer)
library(ggpubr)
library(knitr)
library(dplyr)

plot_richness(ps_all_rarefied, x="sample", color="Variables", measures= "Shannon")+facet_grid (~Variables,switch = "both", scales = "free")+ theme(strip.text = element_text(size=10,face="bold",color="black",lineheight=5.0,angle=90 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black"))

##

ps_m <- prune_taxa(taxa_sums(ps_all_rarefied) > 0, ps_all_rarefied)
tab <- microbiome::alpha(ps_m, index = "all")
kable(head(tab))
ps1.meta <- meta(ps_m)
kable(head(ps1.meta))
ps1.meta$Shannon <- tab$diversity_shannon 
ps1.meta$InverseSimpson <- tab$diversity_inverse_simpson
Category <- levels(ps1.meta$Variables)
?levels


p1 <- ggviolin(ps1.meta, x = "Treatment_Response", y = "Shannon",
               add = c("boxplot", "mean_sd"), fill = "Variables", add.params = list( width= 0.1), facet.by=c("Visit", "lesional_nonlesional"),panel.labs = list(Visit = c("Pre-nbUVB", "Post-nbUVB"), lesional_nonlesional = c("Normal", "Tumor")), palette = c("#5D3A9B", "#E69F00", "#56B4E9", "#009E73", "#66CC99", "#0072B2", "#000099","#CC79A7", "#5D3A9B", "#E69F00", "#56B4E9", "#009E73", "#66CC99", "#0072B2", "#000099","#CC79A7"))+theme(strip.text = element_text(size=10, face="bold", color="black",lineheight=5.0,angle=0 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black"))

print(p1)

library(ggpubr)
library(rstatix)

ggpaired(ps1.meta, x = "Treatment_Response", y = "Shannon", line.color = "gray", line.size = 0.4,palette = c("#5D3A9B", "#E69F00", "#56B4E9", "#009E73", "#66CC99", "#0072B2", "#000099","#CC79A7", "#5D3A9B", "#E69F00", "#56B4E9", "#009E73", "#66CC99", "#0072B2", "#000099","#CC79A7"), add = c("boxplot", "mean_sd"), fill = "Variables", add.params = list(fill = "white"), facet.by=c("Visit", "lesional_nonlesional"),panel.labs = list(Visit = c("Pre-nbUVB", "Post-nbUVB"), lesional_nonlesional = c("Normal", "Tumor")))+theme(strip.text = element_text(size=10, face="bold", color="black",lineheight=5.0,angle=0 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black"))


ps_les <- prune_taxa(taxa_sums(ps_v1_les.nbUVB.rarefied) > 0, ps_v1_les.nbUVB.rarefied)
tab_les <- microbiome::alpha(ps_les, index = "all")
kable(head(tab_les))
ps_les.meta <- meta(ps_les)
kable(head(ps_les.meta))
ps_les.meta$Shannon <- tab_les$diversity_shannon 
ps_les.meta$InverseSimpson <- tab_les$diversity_inverse_simpson

ggpaired(ps_les.meta, x = "Visit", y = "Shannon", line.color = "gray", line.size = 0.4,palette = c("#5D3A9B", "#E69F00", "#56B4E9", "#009E73", "#66CC99", "#0072B2", "#000099","#CC79A7", "#5D3A9B", "#E69F00", "#56B4E9", "#009E73", "#66CC99", "#0072B2", "#000099","#CC79A7"), add = c("boxplot", "mean_sd"), fill = "Variables", add.params = list(fill = "white"), facet.by=c("Treatment_Response"),panel.labs = list(Visit = c("Pre-nbUVB", "Post-nbUVB")))+theme(strip.text = element_text(size=10, face="bold", color="black",lineheight=5.0,angle=0 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black"))


ggpaired(ps1.meta, x = "lesional_nonlesional", y = "Shannon", line.color = "black", line.size = 0.2, add = "mean", fill = "lesional_nonlesional", facet.by=c("Visit"), panel.labs = list(Visit = c("Pre-nbUVB", "Post-nbUVB")))+theme(strip.text = element_text(size=15, face="bold", color="black",lineheight=5.0,angle=0 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black"))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=6, label = paste("Mean=",mean(x))), geom="text") + theme(legend.position="none")

ggpaired(ps_les.meta, x = "Response", y = "Shannon", line.color = "black", line.size = 0.2, add = "mean", fill = "Response")+theme(strip.text = element_text(size=15, face="bold", color="black",lineheight=5.0,angle=0 ,hjust=1,vjust=1),strip.background = element_rect(fill="lightblue", colour="black"))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=6, label = paste("Mean=",round(mean(x), digits = 2))), geom="text") + theme(legend.position="none")

## alpha diversity plots 

##V1 nbUVB-lesional : R vs NR 0.78 p
plot_richness(ps_v1_les.nbUVB.rarefied, x="Response", measures= c("Shannon", "Observed"), title = "Visit 1-Lesional nbUVB Samples: Responders vs Non Responders", color = "Response")+geom_boxplot(width = 0.08)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=", round(mean(x), digits = 2))), geom="text") + theme(legend.position="none")
##V1 nbUVB nonlesional : R vs NR 0.51 p
plot_richness(ps_v1_non.nbUVB.rarefied, x="Response", measures= "Shannon", title = "Visit 1-Non Lesional nbUVB Samples: Responders vs Non Responders")+
  geom_boxplot(width = 0.08)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=", round(mean(x), digits = 2))), geom="text") + theme(legend.position="none")
##V2 nbUVB-lesional : R vs NR 0.018 p
plot_richness(ps_v2_lesional_nbUVB.rarefied, x="Response", measures= "Shannon", title = "Visit 2-Lesional nbUVB Samples: Responders vs Non Responders")+
  geom_boxplot(width = 0.08)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=", round(mean(x), digits = 2))), geom="text") + theme(legend.position="none")
##V2 nbUVB nonlesional : R vs NR 0.026p
plot_richness(ps_v2_non_nbUVB.rarefied, x="Response", measures= "Shannon", title = "Visit 2-Non Lesional nbUVB Samples: Responders vs Non Responders")+
  geom_boxplot(width = 0.08)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=", round(mean(x), digits = 2))), geom="text") + theme(legend.position="none")
## V2 - T_R vs_Nt_NR 0.0085p
plot_richness(R_vs_NT_NR.rarefied, x="Treatment_Response", measures= "Shannon", title = "Visit 2 Lesional: Responders vs Not Treated-Non Responders")+
  geom_boxplot(width = 0.08)+
  theme(legend.position="none",axis.text.x=element_text(angle=0,hjust=1,vjust=1,size=10))+ stat_compare_means(method = "wilcox.test", paired = F,  label = "p", label.x = 0.8)+ stat_summary(fun.data = function(x) data.frame(y=5, label = paste("Mean=", round(mean(x), digits = 2))), geom="text") + theme(legend.position="none")

#beta diversity
##all samples
dist= phyloseq::distance(ps_all_rarefied, method="bray", type= "samples")
ordination=ordinate(ps_all_rarefied, method="PCoA", distance=dist)
plot_ordination(ps_all_rarefied,ordination,color="Visit", title= "Beta Diversity-All samples Visit 1 vs Visit 2") +
  theme_classic() +
  theme(strip.background = element_blank())
#
#v1_nbUVB_les
distv1L = phyloseq::distance(ps_v1_les.nbUVB.rarefied, method="bray")
ordination.v1L = ordinate(ps_v1_les.nbUVB.rarefied, method="PCoA", distance=distv1L)
plot_ordination(ps_v1_les.nbUVB.rarefied, ordination.v1L, color="Treatment_Response", title= "Visit 1- nbUVB Lesional Samples: Responders vs Non Responders") + 
  theme_classic() +
  theme(strip.background = element_blank())
#v1 nbUVB Nonles
distv1N = phyloseq::distance(ps_v1_non.nbUVB.rarefied, method="bray")
ordination.v1N = ordinate(ps_v1_non.nbUVB.rarefied, method="PCoA", distance=distv1N)
plot_ordination(ps_v1_non.nbUVB.rarefied, ordination.v1N, color="Response", title= "Visit 1 nbUVB-Nonlesional Samples: Responders vs Non Responders") + 
  theme_classic() +
  theme(strip.background = element_blank())

# v2_lesional_nbUVB 
distv2.les = phyloseq::distance(ps_v2_lesional_nbUVB.rarefied, method="bray")
ordination.v2.les = ordinate(ps_v2_lesional_nbUVB.rarefied, method="PCoA", distance=distv2.les)
plot_ordination(ps_v2_lesional_nbUVB.rarefied, ordination.v2.les, color="Response", title= "Visit 2- nbUVB Lesional Samples: Responders vs Non Responders") + 
  theme_classic() +
  theme(strip.background = element_blank()) 

### updated beta 

#1st Beta Diversity plot
ps_ord <- ordinate(ps_all_nbUVB.rarefied, "PCoA", "bray", weighted = TRUE)
plot_ordination(ps_all_nbUVB.rarefied, ps_ord, type = "samples",color = "Visit", title= "All Treated Samples Visit 1 vs Visit 2") + theme(plot.title= element_text(family = "Mono",color="blue",size=10))+ geom_point(size = 2)+ theme_classic() +stat_ellipse(aes(group = Visit), linetype = 1)  +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
#2nd Beta Diversity plot        
ps_ord <- ordinate(ps_v1_les.nbUVB.rarefied, "PCoA", "bray", weighted = TRUE)
plot_ordination(ps_v1_les.nbUVB.rarefied, ps_ord, type = "samples",color = "Response", title= "Visit 1: Treated Lesional Responders vs Non Responders") + theme(plot.title= element_text(family = "Mono",color="blue",size=10))+ geom_point(size = 2)+ theme_classic() +stat_ellipse(aes(group = Response), linetype = 1)  +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
#3rd Beta Diversity plot
ps_ord <- ordinate(ps_v1_non.nbUVB.rarefied, "PCoA", "bray", weighted = TRUE)
plot_ordination(ps_v1_non.nbUVB.rarefied, ps_ord, type = "samples",color = "Response", title= "Visit 1: Treated Non-Lesional Responders vs Non Responders") + theme(plot.title= element_text(family = "Mono",color="blue",size=10))+ geom_point(size = 2)+ theme_classic() +stat_ellipse(aes(group = Response), linetype = 1)  +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#4th Beta Diversity plot
ps_ord <- ordinate(ps_v2_lesional_nbUVB.rarefied, "PCoA", "bray", weighted = TRUE)
plot_ordination(ps_v2_lesional_nbUVB.rarefied, ps_ord, type = "samples",color = "Response", title= "Visit 2: Treated Non- Lesional Responders vs Non Responders") + theme(plot.title= element_text(family = "Mono",color="blue",size=10))+ geom_point(size = 2)+ theme_classic() +stat_ellipse(aes(group = Response), linetype = 1)  +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


##
ps_ord <- ordinate(ps_v2_non_nbUVB.rarefied, "PCoA", "bray", weighted = TRUE)
plot_ordination(ps_v2_non_nbUVB.rarefied, ps_ord, type = "samples",color = "Response", title= "Visit 2: Treated Non-Lesional Responders vs Non Responders") + theme(plot.title= element_text(family = "Mono",color="blue",size=10))+ geom_point(size = 2)+ theme_classic() +stat_ellipse(aes(group = Response), linetype = 1)  +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))







####-----MetagenomeSeq analysis-----####
BiocManager::install("metagenomeSeq", force= TRUE)
BiocManager::install("biomformat", force = TRUE)
BiocManager::install("glmnet", force = TRUE)
BiocManager::install("Matrix", force = TRUE)
BiocManager::install("ggplot2", force = TRUE)
BiocManager::install("miaViz", force = TRUE)
library(metagenomeSeq); packageVersion("metagenomeSeq")
library(biomformat); packageVersion("biomformat")
library (phyloseq); packageVersion("phyloseq")
library (ggplot2)
library (miaViz)
library (Matrix)
library (glmnet)

## ----loadData---------------------------------------------
#read the OTU table as a data frame
biom <-read.csv(file="~/Desktop/Rproject/new_biom.txt",sep="\t", header= TRUE, row.names=1)
#convert to biom
ASV_biom<- make_biom(biom)

View (clin)
data.OTU = read.csv(file= "~/Desktop/Rproject/data.otu.txt", sep="\t", header= TRUE, row.names = 1)
dim(data.OTU)
## ----loadTaxa---------------------------------------------
taxa = read.delim("~/Desktop/Rproject/taxo.txt",stringsAsFactors=FALSE,sep="\t", row.names = 1)

## ----loadMeta---------------------------------------------
clin = read.csv("~/Desktop/Rproject/final_MD.1000.txt", sep="\t", header= TRUE, row.names=1)
ord = match(colnames(data.OTU), rownames(clin)) 
clin = clin[ord, ]
head(clin[1:2, ])
phenotypeData = AnnotatedDataFrame(clin)
phenotypeData
OTUdata = AnnotatedDataFrame(taxa)
OTUdata

#Creating a MRexperiment object
obj = newMRexperiment(data.OTU,phenoData=phenotypeData,featureData=OTUdata)

## ----pdata------------------------------------------------
phenoData(obj)
head(pData(obj),3)
dim (obj)

## ----fdata------------------------------------------------
featureData(obj)
head(fData(obj)[,-c(2,10)],3)

## ----MRcounts---------------------------------------------
head(metagenomeSeq::MRcounts(obj[,]))


## ------------------------Subsets------------------------------
head( obj)
####----nbUVB (worsened vs improved)----####
features.nbUVB = which(rowSums(obj) >= 4)
ToKeep = which(pData(obj)$Treatment=="nbUVB")
obj_nbUVB = obj[features.nbUVB,ToKeep]
#normalize 
obj_nbUVB<- cumNorm(obj_nbUVB, p = 0.5)
pd.nbUVB <- pData(obj_nbUVB)
#create a model
mod.nbUVB <- model.matrix(~1+ Response, data = pd.nbUVB) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.nbUVB = fitFeatureModel(obj_nbUVB, mod.nbUVB)

#get coefficients 
MRcoefs.nbUVB <- MRcoefs(objres.nbUVB, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
#export table with coefs
write.table(MRcoefs.nbUVB, file = "~/Desktop/RProject/nbUVB_all_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE)


# take only significant adj pvalues; subset OTU & metadata
significant.nbUVB <- which(MRcoefs.nbUVB$adjPvalues <"0.05")
nbUVB.sig = MRcoefs.nbUVB[significant.nbUVB, ]
nbUVB.ord <- which (rownames(OTUdata) %in% rownames(nbUVB.sig))
nbUVB_otu.sig <- obj_nbUVB[rownames(nbUVB.sig), ] 

# export taxa table to corresponding significant asvs
taxa.all_NBUVB <- taxa [ rownames (nbUVB_otu.sig), ]
write.table(taxa.all_NBUVB, file = "~/Desktop/RProject/nbUVB_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)



###---subset ps obj -based on ASVs
nbUVB_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(nbUVB_otu.sig), ps_all_samples)
taxatokeep <- rownames(nbUVB_otu.sig)
nbUVB_ps.sig <- prune_taxa(taxatokeep, nbUVB_ps)
taxa_names (nbUVB_ps.sig)
sample_data (nbUVB_ps)


#set a df with just adjusted pvalues 
nbUVB.annotate <- nbUVB.sig [rownames(nbUVB_otu.sig), 4]
nbUVB.annotate <- round (nbUVB.annotate, 5 )
dat_text <- data.frame(
  label = nbUVB.annotate,
  OTU = rownames(nbUVB.sig)
)
#Plot
phyloseq::taxa_names(nbUVB_ps.sig) 
phyloseq::psmelt(nbUVB_ps.sig) %>%
  ggplot(aes(x = Response, y = Abundance)) + ggtitle(label = "nbUVB: Worsened vs Improved")+
  geom_boxplot(outlier.shape = NA) +   theme(legend.position="none")+
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free") + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -0.03,
    vjust   = -6.0, size= 3.5
  )


####-----V1-nbUVB Improved vs no_nbUVB_Worsened-----

features.com = which(rowSums(obj) >= 4)
ToKeep.com = which(pData(obj)$Visit== "1" & pData(obj)$Treatment_Response=="nbUVB_Improved"|pData(obj)$Visit== "1" & pData(obj)$Treatment=="No_nbUVB")
obj_both = obj[features.com,ToKeep.com]
#normalize 
obj_both<- cumNorm(obj_both, p = 0.5)
pd.both <- pData(obj_both)
#create a model
mod.both <- model.matrix(~1+ Treatment, data = pd.both) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.both = fitFeatureModel(obj_both, mod.both)

#get coefficients 
MRcoefs.both <- MRcoefs(objres.both, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
#export table with coefs
write.table(MRcoefs.both, file = "~/Desktop/RProject/T_R_NT_allmtgs.txt", col.names=TRUE, sep = "\t", quote = FALSE)


# take only significant adj pvalues; subset OTU & metadata
significant.both <- which(MRcoefs.both$adjPvalues <"0.05")
both.sig = MRcoefs.both[significant.both, ]

# export taxa table to corresponding significant asvs
taxa.both <- taxa [ rownames (both.sig), ]
write.table(taxa.both, file = "~/Desktop/RProject/both_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)

#subset the edgeR table 
V1_NT_R.t<-which (rownames(both.sig) %in% rownames (v1_NT_vs_R_new))
V1_NT_R <- both.sig[V1_NT_R.t, ]

#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 
both.ord <- which (rownames(OTUdata) %in% rownames(V1_NT_R ))
both_otu.sig <- obj_both[rownames(V1_NT_R), ] 

###---subset phyloseq obj -based on ASVs
both_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(both_otu.sig), ps_all_samples)
taxatokeep <- rownames(both_otu.sig)
both_ps.sig <- prune_taxa(taxatokeep, both_ps)
taxa_names (both_ps.sig)
sample_data (both_ps)

#set a df with just adjusted pvalues 
annotate <- both.sig [rownames(V1_NT_R ), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(V1_NT_R)
)

#Plot
phyloseq::taxa_names(both_ps.sig) 
phyloseq::psmelt(both_ps.sig) %>%
  ggplot(aes(x = Treatment, y = Abundance)) + ggtitle(label = "Visit 1: Not Treated vs Treated Improved")+
  geom_boxplot(outlier.shape = NA) +   theme(legend.position="none")+
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+ annotate("text", x = -Inf, y = Inf, label = "FDR adj p-value= ", hjust = -0.15, vjust = 30) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.9,
    vjust   = -15, size= 3.8
  )

####-------V2- nbUVB Improved vs no_nbUVB_######
features.com = which(rowSums(obj) >= 4)
ToKeep.com = which(pData(obj)$Visit== "2" & pData(obj)$Treatment_Response=="nbUVB_Improved"|pData(obj)$Visit== "2" & pData(obj)$Treatment=="No_nbUVB")
obj_both2 = obj[features.com,ToKeep.com]
#normalize 
obj_both2<- metagenomeSeq::cumNorm(obj_both2, p = 0.5)
pd.both <- pData(obj_both2)
#create a model
mod.both <- model.matrix(~1+ Treatment, data = pd.both) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.both2 = metagenomeSeq::fitFeatureModel(obj_both2, mod.both)

#get coefficients 
MRcoefs.both2 <- metagenomeSeq::MRcoefs(objres.both2, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
#export table with coefs
write.table(MRcoefs.both2, file = "~/Desktop/RProject/V2_NT_vs_R_mtgs.txt", col.names=TRUE, sep = "\t", quote = FALSE)


# take only significant adj pvalues; subset OTU & metadata
significant.both2 <- which(MRcoefs.both2$adjPvalues <"0.05")
both2.sig = MRcoefs.both2[significant.both2, ]

#subset the edgeR table 
V2_NT_R.t <-which (rownames(both2.sig) %in% rownames (v2_NT_vs_R_new))
V2_NT_R <- both2.sig[V2_NT_R.t, ]

write.table(V2_NT_R, file = "~/Desktop/RProject/new_matched_v2_NT_R.txt", col.names=TRUE, sep = "\t", quote = FALSE)


#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 
both2.ord <- which (rownames(OTUdata) %in% rownames(both2.sig))
both2_otu.sig <- obj_both2[rownames(both2.sig), ] 

# export taxa table to corresponding significant asvs
taxa.both2 <- taxa [ rownames (both2_otu.sig), ]
write.table(taxa.both2, file = "~/Desktop/RProject/both2_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)

###---subset ps obj -based on ASVs
both2_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(both2_otu.sig), ps_all_samples)
taxatokeep <- rownames(both2_otu.sig)
both2_ps.sig <- prune_taxa(taxatokeep, both2_ps)
taxa_names (both2_ps.sig)
sample_data (both_ps)

#set a df with just adjusted pvalues 
annotate <- both2.sig [rownames(V2_NT_R), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(V2_NT_R)
)

#Plot
phyloseq::taxa_names(both2_ps.sig) 
phyloseq::psmelt(both2_ps.sig) %>%
  ggplot(aes(x = Treatment, y = Abundance)) + ggtitle(label = "Visit 2: Not Treated vs Treated Improved")+
  geom_boxplot(outlier.shape = NA) +   theme(legend.position="none")+
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+ annotate("text", x = -Inf, y = Inf, label = "FDR adj p-value= ", hjust = -0.16, vjust = 21) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.5,
    vjust   = -3.7, size= 3.5
  )
####-------Visit 1- nbUVB Improved vs no_nbUVB_no_change----###
features.com = which(rowSums(obj) >= 4)
ToKeep.com = which(pData(obj)$Visit== "1" & pData(obj)$Treatment_Response=="nbUVB_Improved"|pData(obj)$Visit== "1" & pData(obj)$Treatment_Response=="No_nbUVB_No_Change")
obj_both3 = obj[features.com,ToKeep.com]
#normalize 
obj_both3<- metagenomeSeq::cumNorm(obj_both3, p = 0.5)
pd.both <- pData(obj_both3)
#create a model
mod.both <- model.matrix(~1+ Treatment_Response, data = pd.both) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.both3 = metagenomeSeq::fitFeatureModel(obj_both3, mod.both)

#get coefficients 
MRcoefs.both3 <- metagenomeSeq::MRcoefs(objres.both3, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
#export table with coefs
write.table(MRcoefs.both3, file = "~/Desktop/RProject/both3_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE)


# take only significant adj pvalues; subset OTU & metadata
significant.both3 <- which(MRcoefs.both3$adjPvalues <"0.05")
both3.sig = MRcoefs.both3[significant.both3, ]
both3.ord <- which (rownames(OTUdata) %in% rownames(both3.sig))
both3_otu.sig <- obj_both3[rownames(both3.sig), ] 

# export taxa table to corresponding significant asvs
taxa.both3 <- taxa [ rownames (both3_otu.sig), ]
write.table(taxa.both3, file = "~/Desktop/RProject/both3_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)

###---subset ps obj -based on ASVs
both3_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(both3_otu.sig), ps_all_samples)
taxatokeep <- rownames(both3_otu.sig)
both3_ps.sig <- prune_taxa(taxatokeep, both3_ps)
taxa_names (both3_ps.sig)
sample_data (both3_ps)

#set a df with just adjusted pvalues 
annotate <- both3.sig [rownames(both3_otu.sig), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(both3.sig)
)

#Plot
phyloseq::taxa_names(both3_ps.sig) 
phyloseq::psmelt(both3_ps.sig) %>%
  ggplot(aes(x = Response, y = Abundance)) + ggtitle(label = "Visit 1: No nbUVB-No Change vs nbUVB-Improved")+
  geom_boxplot(outlier.shape = NA) +   theme(legend.position="none")+
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+ annotate("text", x = -Inf, y = Inf, size= 3.5,label = "FDR adj p-value= ", hjust = -0.16, vjust = 6.7) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.5,
    vjust   = -4.0, size= 3.5
  )
####---- V1: no_nbUVB_no_chnage vs _nbUVB_Improved----####
features.com = which(rowSums(obj) >= 4)
ToKeep.com = which(pData(obj)$Visit== "1" & pData(obj)$Treatment_Response=="nbUVB_Improved"|pData(obj)$Visit== "1" & pData(obj)$Treatment_Response=="No_nbUVB_No_Change")
obj_both3 = obj[features.com,ToKeep.com]
#normalize 
obj_both3<- metagenomeSeq::cumNorm(obj_both3, p = 0.5)
pd.both <- pData(obj_both3)
#create a model
mod.both <- model.matrix(~1+ Treatment_Response, data = pd.both) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.both3 = metagenomeSeq::fitFeatureModel(obj_both3, mod.both)

#get coefficients 
MRcoefs.both3 <- metagenomeSeq::MRcoefs(objres.both3, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
#export table with coefs
write.table(MRcoefs.both3, file = "~/Desktop/RProject/both3_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE)


# take only significant adj pvalues; subset OTU & metadata
significant.both3 <- which(MRcoefs.both3$adjPvalues <"0.05")
both3.sig = MRcoefs.both3[significant.both3, ]

# export taxa table to corresponding significant asvs
taxa.both3 <- taxa [ rownames (both3.sig), ]
write.table(taxa.both3, file = "~/Desktop/RProject/both3_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)
nrow (both3.sig)
#subset the edgeR table 
NT_NC_T_I.t <-which (rownames(both3.sig) %in% rownames (v1_NC_imp))
NT_NC_T_I <- both3.sig[NT_NC_T_I.t, ]
rownames (NT_NC_T_I)
#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 

both3.ord <- which (rownames(OTUdata) %in% rownames(NT_NC_T_I))
both3_otu.sig <- obj_both3[rownames(NT_NC_T_I), ] 

###---subset ps obj -based on ASVs
both3_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(both3_otu.sig), ps_all_samples)
taxatokeep <- rownames(both3_otu.sig)
both3_ps.sig <- prune_taxa(taxatokeep, both3_ps)
taxa_names (both3_ps.sig)

#set a df with just adjusted pvalues 
annotate <- both3.sig [rownames(NT_NC_T_I), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(NT_NC_T_I)
)

#Plot
phyloseq::taxa_names(both3_ps.sig) 
phyloseq::psmelt(both3_ps.sig) %>%
  ggplot(aes(x = Response, y = Abundance)) + ggtitle(label = "Visit 1: No nbUVB-No Change vs nbUVB-Improved")+
  geom_boxplot(outlier.shape = NA) +   theme(legend.position="none")+
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "",  y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+ annotate("text", x = -Inf, y = Inf, size= 3.5,label = "FDR adj p-value= ", hjust = -0.16, vjust = 10) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.5,
    vjust   = -3.9, size= 3.5
  )

####---- V2: no_nbUVB_no_chnage vs _nbUVB_Improved----####
features.com = which(rowSums(obj) >= 4)
ToKeep.com = which(pData(obj)$Visit== "2" & pData(obj)$Treatment_Response=="nbUVB_Improved"|pData(obj)$Visit== "2" & pData(obj)$Treatment_Response=="No_nbUVB_No_Change")
obj_both4 = obj[features.com,ToKeep.com]
#normalize 
obj_both4<- metagenomeSeq::cumNorm(obj_both4, p = 0.5)
pd.both <- pData(obj_both4)
#create a model
mod.both <- model.matrix(~1+ Treatment_Response, data = pd.both) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.both4 = metagenomeSeq::fitFeatureModel(obj_both4, mod.both)

#get coefficients 
MRcoefs.both4 <- metagenomeSeq::MRcoefs(objres.both4, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
#export table with coefs
write.table(MRcoefs.both4, file = "~/Desktop/RProject/both4_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE)


# take only significant adj pvalues; subset OTU & metadata
significant.both4 <- which(MRcoefs.both4$adjPvalues <"0.05")
both4.sig = MRcoefs.both4[significant.both4, ]

# export taxa table to corresponding significant asvs
taxa.both4 <- taxa [ rownames (both4_otu.sig), ]
write.table(taxa.both4, file = "~/Desktop/RProject/both4_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)
#subset the edgeR table 
V2_NT_NC_T_I.t <-which (rownames(both4.sig) %in% rownames (v2_NC_imp))
V2_NT_NC_T_I <- both4.sig[V2_NT_NC_T_I.t, ]
rownames (V2_NT_NC_T_I)
#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 

both4.ord <- which (rownames(OTUdata) %in% rownames(V2_NT_NC_T_I ))
both4_otu.sig <- obj_both4[rownames(V2_NT_NC_T_I), ] 

###---subset ps obj -based on ASVs
both4_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(both4_otu.sig), ps_all_samples)
taxatokeep <- rownames(both4_otu.sig)
both4_ps.sig <- prune_taxa(taxatokeep, both4_ps)
taxa_names (both4_ps.sig)
sample_data (both4_ps)

#set a df with just adjusted pvalues 
annotate <- both4.sig [rownames(V2_NT_NC_T_I), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(V2_NT_NC_T_I)
)

#Plot
phyloseq::taxa_names(both4_ps.sig) 
phyloseq::psmelt(both4_ps.sig) %>%
  ggplot(aes(x = Response, y = Abundance)) + ggtitle(label = "Visit 2: No nbUVB-No Change vs nbUVB-Improved")+
  geom_boxplot(outlier.shape = NA) +   theme(legend.position="none")+
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "",  y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+ annotate("text", x = -Inf, y = Inf, size= 3.5,label = "FDR adj p-value= ", hjust = -0.16, vjust = 10.5) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.5,
    vjust   = -4.2, size= 3.5
  )

####---- V1: NT vs _nbUVB_Worsened----####
features.com = which(rowSums(obj) >= 4)
ToKeep.com = which(pData(obj)$Visit== "1" & pData(obj)$Treatment_Response=="nbUVB_Worsened"|pData(obj)$Visit== "1" & pData(obj)$Treatment=="No_nbUVB")
obj_both5 = obj[features.com,ToKeep.com]

#normalize 
obj_both5<- metagenomeSeq::cumNorm(obj_both5, p = 0.5)
pd.both <- pData(obj_both5)
#create a model
mod.both <- model.matrix(~1+ Treatment, data = pd.both) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.both5 = metagenomeSeq::fitFeatureModel(obj_both5, mod.both)

#get coefficients 
MRcoefs.both5 <- metagenomeSeq::MRcoefs(objres.both5, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
#export table with coefs
write.table(MRcoefs.both5, file = "~/Desktop/RProject/new_v1_NT_vsNR.txt", col.names=TRUE, sep = "\t", quote = FALSE)


# take only significant adj pvalues; subset OTU & metadata
significant.both5 <- which(MRcoefs.both5$adjPvalues <"0.05")
both5.sig = MRcoefs.both5[significant.both5, ]
nrow (both5.sig)
# export taxa table to corresponding significant asvs
taxa.both5 <- taxa [ rownames (both5.sig), ]
write.table(taxa.both5, file = "~/Desktop/RProject/1NT_NC_T_W_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)
nrow (taxa.both5)
#subset the edgeR table 
NT_vs_T_W.t <-which (rownames(both5.sig) %in% rownames (v1_NT_vs_NR_new))
NT_vs_T_W <- both5.sig[NT_vs_T_W.t, ]
nrow (NT_vs_T_W)

write.table(NT_vs_T_W, file = "~/Desktop/RProject/new_matched_v1_NT_NR.txt", col.names=TRUE, sep = "\t", quote = FALSE)
#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 
both5.ord <- which (rownames(OTUdata) %in% rownames(NT_vs_T_W))
both5_otu.sig <- obj_both5[rownames(NT_vs_T_W), ] 

###---subset ps obj -based on ASVs
both5_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(both5_otu.sig), ps_all_samples)
taxatokeep <- rownames(both5_otu.sig)
both5_ps.sig <- prune_taxa(taxatokeep, both5_ps)
taxa_names (both5_ps.sig)


#set a df with just adjusted pvalues 
annotate <- both5.sig [rownames(NT_vs_T_W), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(NT_vs_T_W)
)

#Plot
phyloseq::taxa_names(both5_ps.sig) 
phyloseq::psmelt(both5_ps.sig) %>%
  ggplot(aes(x = Treatment, y = Abundance)) + ggtitle(label = "Visit 1: Not-Treated vs Treated NR")+
  geom_boxplot(outlier.shape = NA) +   theme(legend.position="none")+
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "",  y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+ annotate("text", x = -Inf, y = Inf, size= 3.3,label = "FDR adj p-value= ", hjust = -0.16, vjust = 15) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.5,
    vjust   = -4.4, size= 3.5
  )

####---- V2: NT vs _nbUVB_Worsened----####
features.com = which(rowSums(obj) >= 4)
ToKeep.com = which(pData(obj)$Visit== "2" & pData(obj)$Treatment_Response=="nbUVB_Worsened"|pData(obj)$Visit== "2" & pData(obj)$Treatment=="No_nbUVB")
obj_both6 = obj[features.com,ToKeep.com]
#normalize 
obj_both6<- metagenomeSeq::cumNorm(obj_both6, p = 0.5)
pd.both <- pData(obj_both6)
#create a model
mod.both <- model.matrix(~1+ Treatment, data = pd.both) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.both6 = metagenomeSeq::fitFeatureModel(obj_both6, mod.both)

#get coefficients 
MRcoefs.both6 <- metagenomeSeq::MRcoefs(objres.both6, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
#export table with coefs
write.table(MRcoefs.both6, file = "~/Desktop/RProject/new_v2_NT_ns_NR.txt", col.names=TRUE, sep = "\t", quote = FALSE)


# take only significant adj pvalues; subset OTU & metadata
significant.both6 <- which(MRcoefs.both6$adjPvalues <"0.05")
both6.sig = MRcoefs.both6[significant.both6, ]
nrow (V2_NT_NC_T_W )
# export taxa table to corresponding significant asvs
taxa.both6 <- taxa [ rownames (both6.sig), ]
write.table(taxa.both6, file = "~/Desktop/RProject/both6_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)

#subset the edgeR table 
V2_NT_vs_T_W.t <-which (rownames(both6.sig) %in% rownames (v2_NT_vs_NR_new))
V2_NT_vs_T_W <- both6.sig[V2_NT_vs_T_W.t, ]

write.table(V2_NT_vs_T_W, file = "~/Desktop/RProject/new_matched_v2_NT_NR.txt", col.names=TRUE, sep = "\t", quote = FALSE)

#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 
both6.ord <- which (rownames(OTUdata) %in% rownames(V2_NT_vs_T_W))
both6_otu.sig <- obj_both6[rownames(V2_NT_vs_T_W), ] 


###---subset ps obj -based on ASVs
both6_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(both6_otu.sig), ps_all_samples)
taxatokeep <- rownames(both6_otu.sig)
both6_ps.sig <- prune_taxa(taxatokeep, both6_ps)

#set a df with just adjusted pvalues 
annotate <- both6.sig [rownames(V2_NT_vs_T_W), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(V2_NT_vs_T_W)
)

#Plot
phyloseq::taxa_names(both6_ps.sig) 
phyloseq::psmelt(both6_ps.sig) %>%
  ggplot(aes(x = Treatment, y = Abundance)) + ggtitle(label = "Visit 2: Not Treated vs Treated-NR")+
  geom_boxplot(outlier.shape = NA) +   theme(legend.position="none")+
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "",  y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+ annotate("text", x = -Inf, y = Inf, size= 3.5,label = "FDR adj p-value= ", hjust = -0.16, vjust = 38.3) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.5,
    vjust   = -6.8, size= 3.5
  )



#####---V1 nbUVB worsened  vs Improved----####

features.1 = which(rowSums(obj_nbUVB) >= 4)
ToKeep.1 = which(pData(obj_nbUVB)$Visit=="1")
obj_1 = obj_nbUVB[features.1 ,ToKeep.1]

#normalize

obj_1<- cumNorm(obj_1, p = 0.5)
pd.1 <- pData(obj_1)
mod.1 <- model.matrix(~1+ Response, data = pd.1) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.1 = fitFeatureModel(obj_1, mod.1)
MRcoefs.1 <-MRcoefs(objres.1, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
write.table(MRcoefs.1, file = "~/Desktop/RProject/V1_nbUVB_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE)

## take only significant adj pvalues; subset OTU & metadata
significant.v1nbuvb <- which(MRcoefs.1$adjPvalues <"0.05")
v1nbUVB.sig = MRcoefs.1[significant.v1nbuvb, ]
# export taxa table to corresponding significant asvs
taxa.1_nbuvb <- taxa [ rownames (v1nbUVB.sig), ]
write.table(taxa.1_nbuvb, file = "~/Desktop/RProject/v1_nbUVB_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)
#subset the edgeR table 
v1nbUVB.t <-which (rownames(v1nbUVB.sig) %in% rownames (v1_nbUVB))
v1nbUVB <- v1nbUVB.sig[v1nbUVB.t, ]

#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 
v1nbUVB.ord<- which (rownames(OTUdata) %in% rownames(v1nbUVB))
v1nbUVB.otu.sig <- obj_1[rownames(v1nbUVB), ] 


###---subset ps obj -based on ASVs
v1_nbUVB_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(v1nbUVB.otu.sig), ps_all_samples)
taxatokeep <- rownames(v1nbUVB.otu.sig)
v1_nbUVB_ps.sig <- prune_taxa(taxatokeep, v1_nbUVB_ps)


#set a df with just adjusted pvalues 
annotate <- v1nbUVB.sig [rownames(v1nbUVB), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4 )

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(v1nbUVB)
)

#Plot
phyloseq::taxa_names(v1_nbUVB_ps.sig) 
phyloseq::psmelt(v1_nbUVB_ps.sig) %>%
  ggplot(aes(x = Response, y = Abundance)) + ggtitle(label = " Visit 1 nbUVB: Worsened vs Improved")+
  geom_boxplot(outlier.shape = NA) +   theme(legend.position="none")+
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free") +annotate("text", x = -Inf, y = Inf, size =3.5,label = "FDR adj. p-value= ", hjust = -0.10, vjust = 41.5) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.3,
    vjust   = -4, size= 3.5
  )


#####---V2 nbUVB  worsened vs Improved----####

features.2 = which(rowSums(obj_nbUVB) >= 4)
ToKeep.2 = which(pData(obj_nbUVB)$Visit=="2")
obj_2 = obj_nbUVB[features.2 ,ToKeep.2]

#normalize

obj_2<- cumNorm(obj_2, p = 0.5)
pd.2 <- pData(obj_2)
mod.2 <- model.matrix(~1+ Response, data = pd.2) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.2 = fitFeatureModel(obj_2, mod.2)
MRcoefs.2 <-MRcoefs(objres.2, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
write.table(MRcoefs.2, file = "~/Desktop/RProject/V2_nbUVB_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE)

## take only significant adj pvalues; subset OTU & metadata
significant.v2nbuvb <- which(MRcoefs.2$adjPvalues <"0.05")
v2nbUVB.sig = MRcoefs.2[significant.v2nbuvb, ]
# export taxa table to corresponding significant asvs
taxa.2_nbuvb <- taxa [ rownames (v2nbUVB.sig), ]
write.table(taxa.2_nbuvb, file = "~/Desktop/RProject/v2_nbUVB_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)
#subset the edgeR table 
v2nbUVB.t <-which (rownames(v2nbUVB.sig) %in% rownames (v2_nbUVB))
v2nbUVB <- v2nbUVB.sig[v2nbUVB.t, ]

##
v2nbuvb.ord <- which (rownames(OTUdata) %in% rownames(v2nbUVB))
v2nbUVB.otu.sig <- obj_2[rownames(v2nbUVB), ] 


###---subset ps obj -based on ASVs
v2_nbUVB_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(v2nbUVB.otu.sig), ps_all_samples)
taxatokeep <- rownames(v2nbUVB.otu.sig)
v2_nbUVB_ps.sig <- prune_taxa(taxatokeep, v2_nbUVB_ps)

#set a df with just adjusted pvalues 
annotate <- v2nbUVB.sig [rownames(v2nbUVB), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(v2nbUVB)
)
#Plot
phyloseq::taxa_names(v2_nbUVB_ps.sig) 
phyloseq::psmelt(v2_nbUVB_ps.sig) %>%
  ggplot(aes(x = Response, y = Abundance)) + ggtitle(label = "Visit 2: Worsened nbUVB vs Improved nbUVB")+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+annotate("text", x = -Inf, y = Inf, size =3.5,label = "FDR adj. p-value= ", hjust = -0.10, vjust = 13.5) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.3,
    vjust   = -6, size= 3.5
  )

####----- No nbUVB-----####
features.no_nbUVB = which(rowSums(obj) >= 4)
ToKeep = which(pData(obj)$Treatment=="No_nbUVB" )
obj_no_nbUVB = obj[features.no_nbUVB,ToKeep]
#normalize 
obj_no_nbUVB<- metagenomeSeq::cumNorm(obj_no_nbUVB, p = 0.5)
pd.no_nbUVB <- pData(obj_no_nbUVB)
#create a model
mod.no_nbUVB <- model.matrix(~1+ Response, data = pd.no_nbUVB) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.no_nbUVB = metagenomeSeq::fitFeatureModel(obj_no_nbUVB, mod.no_nbUVB)

#get coefficients 
MRcoefs.no_nbUVB <- metagenomeSeq::MRcoefs(objres.no_nbUVB, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
#export table with coefs
write.table(MRcoefs.no_nbUVB, file = "~/Desktop/RProject/no_nbUVB_all_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE) 

# take only significant adj pvalues; subset OTU & metadata
significant.nonbUVB <- which(MRcoefs.no_nbUVB$adjPvalues <"0.05")
no_nbUVB.sig = MRcoefs.no_nbUVB[significant.nonbUVB, ]
no_nbUVB.ord <- which (rownames(OTUdata) %in% rownames(no_nbUVB.sig))
no_nbUVB_otu.sig <- obj_no_nbUVB[rownames(no_nbUVB.sig), ] 
featureNames(no_nbUVB_otu.sig)
# export taxa table to corresponding significant asvs
taxa.all_noNBUVB <- taxa [ rownames (no_nbUVB_otu.sig), ]
write.table(taxa.all_noNBUVB, file = "~/Desktop/RProject/no_nbUVB_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)


###---subset ps obj -based on ASVs
no_nbUVB_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(no_nbUVB_otu.sig), ps_all_samples)
taxatokeep <- rownames(no_nbUVB_otu.sig)
no_nbUVB_ps.sig <- prune_taxa(taxatokeep, no_nbUVB_ps)
taxa_names (no_nbUVB_ps.sig)

#Plot
phyloseq::taxa_names(no_nbUVB_ps.sig) 
phyloseq::psmelt(no_nbUVB_ps.sig) %>%
  ggplot(aes(x = Response, y = Abundance)) + ggtitle(label = "No nbUVB: Worsened vs No Change")+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")

####----All worsened samples ----####
features.w = which(rowSums(obj) >= 4)
ToKeep.w = which(pData(obj)$Response=="Worsened")
obj_w = obj[features.w ,ToKeep.w]

#fitfeature model 
obj_w<- metagenomeSeq::cumNorm(obj_w, p = 0.5)
pd.w <- pData(obj_w)
mod.w <- model.matrix(~1+ Treatment, data = pd.w) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres = metagenomeSeq::fitFeatureModel(obj_w, mod.w)
MRcoefs.w <- metagenomeSeq::MRcoefs(objres, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
write.table(MRcoefs.w, file = "~/Desktop/RProject/worsened_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE) 

## take only significant adj pvalues; subset OTU & metadata
significant.worsened <- which(MRcoefs.w$adjPvalues <"0.05")
worsened.sig = MRcoefs.w[significant.worsened, ]
worsened.ord <- which (rownames(OTUdata) %in% rownames(worsened.sig))
worsened.otu.sig <- obj_w[rownames(worsened.sig), ] 

# export taxa table to corresponding significant asvs
taxa.all_worsened <- taxa [ rownames (worsened.sig), ]
write.table(taxa.all_worsened, file = "~/Desktop/RProject/all_worsened_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)

###---subset ps obj -based on ASVs
worsened_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(worsened.otu.sig), ps_all_samples)
taxatokeep <- rownames(worsened.otu.sig)
worsened_ps.sig <- prune_taxa(taxatokeep, worsened_ps)

#Plot


phyloseq::taxa_names(worsened_ps.sig) 
phyloseq::psmelt(worsened_ps.sig) %>%
  ggplot(aes(x = Treatment, y = Abundance)) + ggtitle(label = "Worsened: nbUVB vs No nbUVB")+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")

#####---V1 worsened nbUVB vs no nbUVB----####

features.v1 = which(rowSums(obj_w) >= 4)
ToKeep.v1 = which(pData(obj_w)$Visit=="1")
obj_v1 = obj_w[features.v1 ,ToKeep.v1]

#normalize

obj_v1<- metagenomeSeq::cumNorm(obj_v1, p = 0.5)
pd.v1 <- pData(obj_v1)
mod.v1 <- model.matrix(~1+ Treatment, data = pd.v1) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.v1 = metagenomeSeq::fitFeatureModel(obj_v1, mod.v1)
MRcoefs.v1 <-metagenomeSeq::MRcoefs(objres.v1, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
write.table(MRcoefs.v1, file = "~/Desktop/RProject/V1_all_worsened_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE)

## take only significant adj pvalues; subset OTU & metadata
significant.v1worsened <- which(MRcoefs.v1$adjPvalues <"0.05")
v1worsened.sig = MRcoefs.v1[significant.v1worsened, ]

# export taxa table to corresponding significant asvs
taxa.v1_worsened <- taxa [ rownames (v1worsened.sig), ]
write.table(taxa.v1_worsened, file = "~/Desktop/RProject/v1_worsened_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)
#subset the edgeR table 
v1worsened.t <-which (rownames(v1worsened.sig) %in% rownames (v1_w_nbUVB_vs_nnbUVB))
v1worsened <- v1worsened.sig[v1worsened.t, ]
nrow(v1worsened)

#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 
v1worsened.ord <- which (rownames(OTUdata) %in% rownames(v1worsened))
v1worsened.otu.sig <- obj_v1[rownames(v1worsened), ] 

###---subset ps obj -based on ASVs
v1_worsened_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(v1worsened.otu.sig), ps_all_samples)
taxatokeep <- rownames(v1worsened.otu.sig)
v1_worsened_ps.sig <- prune_taxa(taxatokeep, v1_worsened_ps)
taxa_names (v1_worsened_ps.sig)

#set a df with just adjusted pvalues 
annotate <-  v1worsened.sig[rownames(v1worsened), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(v1worsened)
)
#Plot
phyloseq::taxa_names(v1_worsened_ps.sig) 
phyloseq::psmelt(v1_worsened_ps.sig) %>%
  ggplot(aes(x = Treatment, y = Abundance)) + ggtitle(label = "Visit 1: Worsened nbUVB vs no nbUVB")+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+annotate("text", x = -Inf, y = Inf, size =3.5,label = "FDR adj. p-value= ", hjust = -.08, vjust = 31.5) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.3,
    vjust   = -15, size= 3.5
  )


####---V2 worsened nbUVB vs no_nbUVB-----
features.v2 = which(rowSums(obj_w) >= 4)
ToKeep.v2 = which(pData(obj_w)$Visit=="2")
obj_v2 = obj_w[features.v2 ,ToKeep.v2]

#create a model, normalize 

obj_v2<- metagenomeSeq::cumNorm(obj_v2, p = 0.5)
pd.v2 <- pData(obj_v2)
mod.v2 <- model.matrix(~1+ Treatment, data = pd.v2) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.v2 = metagenomeSeq::fitFeatureModel(obj_v2, mod.v2)
MRcoefs.v2 <- metagenomeSeq::MRcoefs(objres.v2, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
write.table(MRcoefs.v2, file = "~/Desktop/RProject/V2_all_worsened_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE) 
# take only significant adj pvalues; subset OTU & metadata
significant.v2worsened <- which(MRcoefs.v2$adjPvalues <"0.05")
v2worsened.sig = MRcoefs.v2[significant.v2worsened, ]

# export taxa table to corresponding significant asvs
taxa.v2_worsened <- taxa [ rownames (v2worsened.sig), ]
write.table(taxa.v2_worsened, file = "~/Desktop/RProject/v2_worsened_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)

#subset the edgeR table 
v2worsened.t <-which (rownames(v2worsened.sig) %in% rownames (v2_w_nbUVB_vs_nnbUVB))
v2worsened <- v2worsened.sig[v2worsened.t, ]
nrow(no_nbUVB_v1)

#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 

v2worsened.ord <- which (rownames(OTUdata) %in% rownames(v2worsened))
v2worsened.otu.sig <- obj_v2[rownames(v2worsened), ] 

###---subset ps obj -based on ASVs
v2_worsened_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(v2worsened.otu.sig), ps_all_samples)
taxatokeep2 <- rownames(v2worsened.otu.sig)
v2_worsened_ps.sig <- prune_taxa(taxatokeep2, v2_worsened_ps)
taxa_names (v2_worsened_ps.sig)

#set a df with just adjusted pvalues 
annotate <-  v2worsened.sig[rownames(v2worsened), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(v2worsened)
)
#Plot
phyloseq::taxa_names(v2_worsened_ps.sig) 
phyloseq::psmelt(v2_worsened_ps.sig) %>%
  ggplot(aes(x = Treatment, y = Abundance)) + ggtitle(label = "Visit 2: Worsened nbUVB vs no nbUVB") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+annotate("text", x = -Inf, y = Inf, size =4.5,label = "FDR adj. p-value= ", hjust = -0.10, vjust = 65) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.7,
    vjust   = -15, size= 3.5
  )


####----nbUVB- worsened_no_ v1 vs v2----##

#subset no_nbUVB- worsened 
feature_no= which(rowSums(obj_no_nbUVB) >= 4)
Tokeep_no = which(pData(obj_no_nbUVB)$Response=="Worsened")
obj_no = obj_no_nbUVB[ ,Tokeep_no]

#
obj_no<- metagenomeSeq::cumNorm(obj_no, p = 0.5)
pd.no <- pData(obj_no)
mod.no <- model.matrix(~1+ Visit, data = pd.no) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.no = metagenomeSeq::fitFeatureModel(obj_no, mod.no)
MRcoefs.no <- metagenomeSeq::MRcoefs(objres.no, by = 2, number=1000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
write.table(MRcoefs.no, file = "~/Desktop/RProject/no_nbUVB_w_coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE) 

# take only significant adj pvalues; subset OTU & metadata
significant.no_nbUVB_w <- which(MRcoefs.no$adjPvalues <"0.05")
no_nbUVB_w.sig = MRcoefs.no[significant.no_nbUVB_w, ]
no_nbUVB_w.ord <- which (rownames(OTUdata) %in% rownames(no_nbUVB_w.sig))
no_nbUVB_w.otu.sig <- obj_no_nbUVB[rownames(no_nbUVB_w.sig), ] 

# export taxa table to corresponding significant asvs
taxa_no_nbUVB_w<- taxa [ rownames (no_nbUVB_w.sig), ]
write.table(taxa_no_nbUVB_w, file = "~/Desktop/RProject/taxa_no_nbUVB_w_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)

#subset ps obj -based on ASVs
no_nbUVB_worsened_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(no_nbUVB_w.otu.sig), ps_all_samples)
taxatokeep <- rownames(no_nbUVB_w.otu.sig)
no_nbUVB_worsened_ps.sig <- prune_taxa(taxatokeep,no_nbUVB_worsened_ps )
taxa_names (v2_worsened_ps.sig)

##set a df with just adjusted pvalues 
annotate <-  no_nbUVB_w.sig[rownames(no_nbUVB_w.otu.sig), 4]

#round off pvalues  to 4 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(no_nbUVB_w.sig)
)
#Plot
phyloseq::taxa_names(no_nbUVB_worsened_ps.sig) 
phyloseq::psmelt(no_nbUVB_worsened_ps.sig) %>%
  ggplot(aes(x = Visit, y = Abundance)) + ggtitle(label = "Worsened-No nbUVB")+ 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Visit", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+annotate("text", x = -Inf, y = Inf, size =3.5,label = "FDR adj. p-value= ", hjust = -.08, vjust = 78.5) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -.5,
    vjust   = -20, size= 3.5
  )

####---V1 NT : No_change vs worsened----
features.v1_all = which(rowSums(obj_no_nbUVB) >= 4)
ToKeep.v1_all = which(pData(obj_no_nbUVB)$Visit=="1")
obj_v1_all = obj_no_nbUVB[features.v1_all ,ToKeep.v1_all]

#v1 - all  
obj_v1_all<- metagenomeSeq::cumNorm(obj_v1_all, p = 0.5)
pd.v1_all <- pData(obj_v1_all)
mod.v1_all <- model.matrix(~1+ Response, data = pd.v1_all) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.v1_all = metagenomeSeq::fitFeatureModel(obj_v1_all, mod.v1_all)
MRcoefs.v1_all <- metagenomeSeq::MRcoefs(objres.v1_all, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
write.table(MRcoefs.v1_all, file = "~/Desktop/RProject/V1_NC_vs_NR.coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE) 

# take only significant adj pvalues; subset OTU & metadata
significant.no_nbUVB_v1 <- which(MRcoefs.v1_all$adjPvalues <"0.05")
no_nbUVB_v1.sig = MRcoefs.v1_all[significant.no_nbUVB_v1, ]

# export taxa table to corresponding significant asvs
taxa_no_nbUVB_v1<- taxa [ rownames (no_nbUVB_v1.sig), ]
write.table(taxa_no_nbUVB_v1, file = "~/Desktop/RProject/taxa_no_nbUVB_v1_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)
#subset the edgeR table 
no_nbUVB_v1.t <-which (rownames(no_nbUVB_v1.sig) %in% rownames (v1_no_nc_no_w))
no_nbUVB_v1 <- no_nbUVB_v1.sig[no_nbUVB_v1.t, ]
nrow(no_nbUVB_v1)

#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 
no_nbUVB_v1.ord <- which (rownames(OTUdata) %in% rownames(no_nbUVB_v1))
no_nbUVB_v1.otu.sig <- obj_v1_all[rownames(no_nbUVB_v1), ] 

###---subset ps obj -based on ASVs
no_nbUVB_v1_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(no_nbUVB_v1.otu.sig), ps_all_samples)
taxatokeep <- rownames(no_nbUVB_v1.otu.sig)
no_nbUVB_v1_ps.sig <- prune_taxa(taxatokeep,no_nbUVB_v1_ps)

##set a df with just adjusted pvalues from metagenomeSeq
annotate <-  no_nbUVB_v1.sig[rownames(no_nbUVB_v1 ), 4]

#round off pvalues to 4 digits 
annotate <- round (annotate, 4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(no_nbUVB_v1)
)

#Plot
phyloseq::taxa_names(no_nbUVB_v1_ps.sig) 
phyloseq::psmelt(no_nbUVB_v1_ps.sig) %>%
  ggplot(aes(x = Response, y = Abundance)) + ggtitle(label = "Visit 1:No_nbUVB- No Change vs Worsened")+ 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Response", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+annotate("text", x = -Inf, y = Inf, size =3.0,label = "FDR adj. p-value= ", hjust = -0.03, vjust = 11) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -0.3,
    vjust   = -4.9, size= 3.7
  )

####------V2 NT: No_Change vs Worsened-----
ToKeep.v2_all = which(pData(obj_no_nbUVB)$Visit=="2")
obj_v2_all = obj_no_nbUVB[features.v1_all ,ToKeep.v2_all]

#v2 - all 
obj_v2_all<- metagenomeSeq::cumNorm(obj_v2_all, p = 0.5)
pd.v2_all <- pData(obj_v2_all)
mod.v2_all <- model.matrix(~1+ Response, data = pd.v2_all) #works for treatment, lesional or non lesional, and visit  but not for categories with more than 3 levels
objres.v2_all = metagenomeSeq::fitFeatureModel(obj_v2_all, mod.v2_all)
MRcoefs.v2_all <- metagenomeSeq::MRcoefs(objres.v2_all, by = 2, number=10000, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)
write.table(MRcoefs.v2_all, file = "~/Desktop/RProject/V2_Nc_vs_w.coefs.txt", col.names=TRUE, sep = "\t", quote = FALSE) 

# take only significant adj pvalues; subset OTU & metadata
significant.no_nbUVB_v2 <- which(MRcoefs.v2_all$adjPvalues <"0.05")
no_nbUVB_v2.sig = MRcoefs.v2_all[significant.no_nbUVB_v2, ]
# export taxa table to corresponding significant asvs
taxa_no_nbUVB_v2<- taxa [ rownames (no_nbUVB_v2.sig), ]
write.table(taxa_no_nbUVB_v2, file = "~/Desktop/RProject/taxa_no_nbUVB_v2_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)
#subset the ASVs that are in the significant ASVs from the edgeR table 
no_nbUVB_v2.t <-which (rownames(no_nbUVB_v2.sig) %in% rownames (v2_no_nc_no_w))
no_nbUVB_v2 <- no_nbUVB_v2.sig[no_nbUVB_v2.t, ]

#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 
no_nbUVB_v2.ord <- which (rownames(OTUdata) %in% rownames(no_nbUVB_v2))
no_nbUVB_v2.otu.sig <- obj_v2_all[rownames(no_nbUVB_v2), ] 

###---subset ps obj -based on ASVs
no_nbUVB_v2_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(no_nbUVB_v2.otu.sig), ps_all_samples)
taxatokeep <- rownames(no_nbUVB_v2.otu.sig)
no_nbUVB_v2_ps.sig <- prune_taxa(taxatokeep,no_nbUVB_v2_ps)
taxa_names (no_nbUVB_v2_ps.sig)

##set a df with just adjusted pvalues 
annotate <-  no_nbUVB_v2.sig[rownames(no_nbUVB_v2), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate,4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(no_nbUVB_v2)
)


#Plot
phyloseq::taxa_names(no_nbUVB_v2_ps.sig) 
phyloseq::psmelt(no_nbUVB_v2_ps.sig) %>%
  ggplot(aes(x = Response, y = Abundance)) + ggtitle(label = "Visit 2:No nbUVB- No Change vs Worsened")+ 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Response", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")+annotate("text", x = -Inf, y = Inf, size =3.0,label = "FDR adj. p-value= ", hjust = -0.03, vjust = 25) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -0.3,
    vjust   = -7.0, size= 3.7
  )

####---- v1v2 NT: No change ------

features.NC = which(rowSums(obj_no_nbUVB) >= 4)
Tokeep_NC = which(pData(obj_no_nbUVB)$Response=="No_Change")
obj_NC = obj_no_nbUVB[features.NC,Tokeep_NC]

# No_change- No_nbUVB
obj_NC<- metagenomeSeq::cumNorm(obj_NC, p = 0.5)
pd.NC <- pData(obj_NC)
mod.NC <- model.matrix(~1+ Visit, data = pd.NC) 
objres.NC= metagenomeSeq::fitFeatureModel(obj_NC, mod.NC)

MRcoefs.NC <- metagenomeSeq::MRcoefs(objres.NC, number= 5000, by = 2, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)

write.table(MRcoefs.NC, file = "~/Desktop/RProject/V1v2_no_nbUVB_NC.txt", col.names=TRUE, sep = "\t", quote = FALSE) 

# take only significant adj pvalues; subset OTU & metadata
significant.no_nbUVB_NC <- which(MRcoefs.NC$adjPvalues <"0.05")
no_nbUVB_NC.sig = MRcoefs.NC[significant.no_nbUVB_NC, ]
# export taxa table to corresponding significant asvs
taxa_no_nbUVB_NC<- taxa [ rownames (no_nbUVB_NC.sig), ]
write.table(taxa_no_nbUVB_NC, file = "~/Desktop/RProject/taxa_no_nbUVB_NC_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)

#subset the edgeR table 
no_nbUVB_NC.t <-which (rownames(no_nbUVB_NC.sig) %in% rownames (v1v2_nnbUVB_nc))
no_nbUVB_NC <- no_nbUVB_NC.sig[no_nbUVB_NC.t, ]

#subset the main metagenomeSeq object to only have the ASVs that are present from edgeR GLM analysis 
no_nbUVB_NC.ord <- which (rownames(OTUdata) %in% rownames(no_nbUVB_NC))
no_nbUVB_NC.otu.sig <- obj_NC[rownames(no_nbUVB_NC), ] 

###---subset ps obj -based on ASVs
no_nbUVB_NC_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(no_nbUVB_NC.otu.sig), ps_all_samples)
taxatokeep <- rownames(no_nbUVB_NC.otu.sig)
no_nbUVB_NC_ps.sig <- prune_taxa(taxatokeep,no_nbUVB_NC_ps)
taxa_names (no_nbUVB_NC_ps.sig)

##set a df with just adjusted pvalues 
annotate <-  no_nbUVB_NC.sig[rownames(no_nbUVB_NC), 4]

#round off pvalues  to 5 digits 
annotate <- round (annotate,4)

#make an object defining the text to annotate on our graph
dat_text <- data.frame(
  label = annotate,
  OTU = rownames(no_nbUVB_NC)
)

#Plot
phyloseq::taxa_names(no_nbUVB_NC_ps.sig) 
phyloseq::psmelt(no_nbUVB_NC_ps.sig) %>%
  ggplot(aes(x = Visit, y = Abundance)) + ggtitle(label = "No nbUVB-No Change: Visit 1 vs Visit 2")+ 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Visit", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free") +annotate("text", x = -Inf, y = Inf, size =3.5,label = "FDR adj. p-value= ", hjust = -0.05, vjust = 77.5) + geom_text(
    data    = dat_text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -0.6,
    vjust   = -19, size= 3.7
  )

#subset worsened- No nbUVB
features.nnb = which(rowSums(obj_w) >= 4)
ToKeep.nnb = which(pData(obj_w)$Treatment=="No_nbUVB")
obj_nnb = obj_w[ ,ToKeep.nnb]

#  Normalize no nbuvb worsened 
obj_nnb<- cumNorm(obj_nnb, p = 0.5)
pd.nnb <- pData(obj_nnb)
mod.nnb <- model.matrix(~1+ Visit, data = pd.nnb) 
objres.nnb= fitFeatureModel(obj_nnb, mod.nnb)

MRcoefs.nnb <- MRcoefs(objres.nnb, number= 10000, by = 2, coef = NULL, uniqueNames = FALSE, adjustMethod = "fdr", group = 0, eff = 0, numberEff = FALSE, counts = 0, file = NULL)

write.table(MRcoefs.nnb, file = "~/Desktop/RProject/worsened_nnb.txt", col.names=TRUE, sep = "\t", quote = FALSE) 

# take only significant adj pvalues; subset OTU & metadata
significant.nnb <- which(MRcoefs.nnb$adjPvalues <"0.05")
nnb.sig = MRcoefs.nnb[significant.nnb, ]
nnb.ord <- which (rownames(OTUdata) %in% rownames(nnb.sig))
nnb.otu.sig <- obj_nnb[rownames(nnb.sig), ] 

# export taxa table to corresponding significant asvs
taxa_nnb<- taxa [ rownames (nnb.sig), ]
write.table(taxa_nnb, file = "~/Desktop/RProject/nnb_taxa.txt", col.names=TRUE, sep = "\t", quote = FALSE)

###---subset ps obj -based on ASVs
nnb_ps <- prune_samples(sample_names(ps_all_samples) %in% colnames(nnb.otu.sig), ps_all_samples)
taxatokeep <- rownames(nnb.otu.sig)
nnb_ps.sig <- prune_taxa(taxatokeep,nnb_ps)
taxa_names (nnb_ps.sig)

#Plot
phyloseq::taxa_names(no_nbUVB_NC_ps.sig) 
phyloseq::psmelt(no_nbUVB_NC_ps.sig) %>%
  ggplot(aes(x = Visit, y = Abundance)) + ggtitle(label = "No Change: Visit 1 vs Visit 2")+ 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "Visit", y = "Abundance\n")+
  facet_wrap(~ OTU, scales = "free")


#### ----edgeR ----####
BiocManager::install("edgeR", force= TRUE)
library (edgeR); packageVersion ("edgeR")
citation ("edgeR")

BiocManager::install("PathoStat", force= TRUE)
library (PathoStat)
install.packages("statmod")
library(statmod)

#convert phyloseq objects to edgeR list 

#y<- PathoStat::phyloseq_to_edgeR(ps_all_samples, group="Treatment")

##load count (OTU) data
rawdata <-read.csv(file= "~/Desktop/Rproject/data.otu.txt", sep="\t", header= TRUE, row.names = 1)
dim(rawdata)
##load metadata 
metadata.edgeR <- read.csv("~/Desktop/Rproject/MD_for_edgeR.txt", sep="\t", header= TRUE, row.names=1)

group = metadata.edgeR$Visit_status
length(group)

#create the dge object
y <- DGEList(counts = rawdata, group = group)
dim (y)
y$samples$group
sum(y$all.zeros)
#make a copy of the original list
y_full <- y
head (y$counts)
##filtering
keep <- filterByExpr(y, min.prop= 0.1, group=group)
table (keep)
y <- y[keep, , keep.lib.sizes=FALSE]#keep only OTUs that pass the filter 
dim (y) 
##effective library sizes
y$samples$lib.size <- colSums (y$counts) 
#We now perform normalization steps, which is totally independent from our experimental design
y <- calcNormFactors(y)
#Now we can see the scaling factors: these should be "reasonably" similar among all samples
y$samples

y <- estimateDisp(y, design)
#Alternatively, one can use the following calling sequence to estimate them one by one. To estimate common dispersion:
y <- estimateGLMCommonDisp(y, design)

#To estimate trended dispersions:
y <- estimateGLMTrendedDisp(y, design)
#To estimate tagwise dispersions:
y <- estimateGLMTagwiseDisp(y, design)

#GLM method 
design <- model.matrix(~0+group, data=y$samples)
rownames(design) <- colnames (y)
design

fit <- glmQLFit(y, design)
fit$coefficients
qlf.2vs1 <- glmQLFTest(fit, contrast=c(0,0,0,0,-1,1))

v1_NT_vs_R_new <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v2_NT_vs_R_new <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v1_NT_vs_NR_new <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v2_NT_vs_NR_new <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)


write.table(v2_NT_vs_NR_new , file = "~/Desktop/RProject/NEW.v2_nbUVB_NR_NT.txt", col.names=TRUE, sep = "\t", quote = FALSE)

v1_no_nc_no_w<- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v2_no_nc_no_w<- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v1_w_nbUVB_vs_nnbUVB<- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v2_w_nbUVB_vs_nnbUVB<- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v1v2_w_nnbUVB<- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v1v2_nnbUVB_nc<- topTags(qlf.2vs1, n=100, sort.by="PValue", p.value=1)
v1_nbuvb_i_no_w <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v1_nbUVB <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v2_nbUVB <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v2_nbuvb_i_no_w <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v1_NC_imp  <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v2_NC_imp  <- topTags(qlf.2vs1, n=300, sort.by="PValue", p.value=1)
v1_NC_worsened  <- topTags(qlf.2vs1, n=500, sort.by="PValue", p.value=1)
v2_NC_worsened  <- topTags(qlf.2vs1, n=500, sort.by="PValue", p.value=1)
#export table with coefs
write.table(v1_no_nc_no_w, file = "~/Desktop/RProject/v1_no_nc_no_w.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v2_no_nc_no_w, file = "~/Desktop/RProject/v2_no_nc_no_w.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v1_w_nbUVB_vs_nnbUVB, file = "~/Desktop/RProject/v1_w_nbUVB_vs_nnbUVB.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v2_w_nbUVB_vs_nnbUVB, file = "~/Desktop/RProject/v2_w_nbUVB_vs_nnbUVB.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v1v2_w_nnbUVB, file = "~/Desktop/RProject/v1v2_w_nnbUVB.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v1_nbuvb_i_no_w, file = "~/Desktop/RProject/v1_nbuvb_i_no_w.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v2_nbuvb_i_no_w, file = "~/Desktop/RProject/v2_nbuvb_i_no_w.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v1_NC_imp, file = "~/Desktop/RProject/v1_NC_imp.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v2_NC_imp, file = "~/Desktop/RProject/v2_NC_imp.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v1_NC_worsened, file = "~/Desktop/RProject/v1_NC_worsened.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v2_NC_worsened, file = "~/Desktop/RProject/V2_NC_worsened.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v1_nbUVB, file = "~/Desktop/RProject/v1_nbUVB_R_NR.txt", col.names=TRUE, sep = "\t", quote = FALSE)
write.table(v2_nbUVB, file = "~/Desktop/RProject/v2_nbUVB_R_NR.txt", col.names=TRUE, sep = "\t", quote = FALSE)



