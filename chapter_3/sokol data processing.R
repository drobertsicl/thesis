# https://benjjneb.github.io/dada2/tutorial.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html

#source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")

library(dada2); packageVersion("dada2") ### Sequence processing 
library(DECIPHER) ### Performs sequence alignment
library(phangorn) ## Phylogenetic tree generation
library(ggplot2) ### Data visualisation and analysis
library(phyloseq) ### Data visualisation and analysis
library(pacman)
library(knitr)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(mia)
library(lme4)
library(ggpubr)
library(vegan)

setwd("C:/Users/eatmo/Desktop/CAC/Sokol/")
getwd()

path <- "C:/Users/eatmo/Desktop/CAC/Sokol/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##Inspect read quality profiles

#Visualize the quality profiles of the forward reads:
plotQualityProfile(fnFs[1:2])

#Visualize the quality profiles of the reverse reads:
plotQualityProfile(fnRs[1:2])

##Filter annd trim

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered2", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered2", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19, 19), truncLen=c(250,190),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

##Learn the error rates
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

#Visualize the estimated error rates:
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#save.image("16S_dada2.Rdata")

##Dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

save.image("16S_dada2.Rdata")

##Merge paired reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

save.image("16S_dada2.Rdata")

##Construct sequence table

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

hist(nchar(getSequences(seqtab)))

#You can remove non-target-length sequences with base R manipulations of the sequence table
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,390)]

# Inspect distribution of new sequence lengths
table(nchar(getSequences(seqtab2)))

hist(nchar(getSequences(seqtab2)))

# Overwrite original table

seqtab <- seqtab2

##Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

##Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

save.image("16S_dada2.Rdata")

##Assign taxonomy

chunk1 <- seqtab.nochim[c(1:21), c(1:500)]
chunk2 <- seqtab.nochim[c(1:21), c(501:1000)]
chunk3 <- seqtab.nochim[c(1:21), c(1001:1439)]

nochimchunk3 <- assignTaxonomy(chunk3, "D:/Sokol_study_data/Silva/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa_silva <- rbind(nochimchunk1, nochimchunk2, nochimchunk3)

taxa_silva <- assignTaxonomy(seqtab.nochim, "C:/Users/eatmo/Downloads/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa_silva_species <- addSpecies(taxa_silva, "C:/Users/eatmo/Desktop/16s/silva_species_assignment_v132.fa.gz")

save.image("16S_dada2.Rdata")

#Inspect the taxonomic assignments:
taxa.print <- taxa_silva # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#If your reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences 
#are being assigned as Eukaryota NA NA NA NA NA, your reads may be in the opposite orientation as the reference database. 
#Tell dada2 to try the reverse-complement orientation and see if this fixes the assignments:
#taxa <- assignTaxonomy(seqtab.nochim, "D:/JULIE_FILES/dada2/dada2_healthy_serial/silva_nr_v128_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)

write.csv(taxa_silva, "taxa_silva_idiot.csv")
write.csv(taxa_silva_species, "taxa_silva_species_idiot.csv")
write.csv(seqtab.nochim, "USV_counts_silva.csv")

save.image("16S_dada2.Rdata")

###FASTA Export
# giving our seq headers more manageable names (ASV_1, ASV_2...)
#asv_seqs <- colnames(seqtab.nochim)
#asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

#for (i in 1:dim(seqtab.nochim)[2]) {
#  asv_headers[i] <- paste(">ASV", i, sep="_")
#}

# making and writing out a fasta of our final ASV seqs:
#asv_fasta <- c(rbind(asv_headers, asv_seqs))
#write(asv_fasta, "ASVs.fa")

# count table:
#asv_tab <- t(seqtab.nochim)
#row.names(asv_tab) <- sub(">", "", asv_headers)
#write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
#asv_tax <- taxa_silva
#row.names(asv_tax) <- sub(">", "", asv_headers)
#write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

# change path to wherever your files are, format might be different on mac
path_to_files <- "C:/Users/eatmo/Desktop/CAC/Sokol"

setwd(path_to_files)

# load in the components to construct our Tree Summarised Experiment object

# There is a difference between read.csv (base R) and read_csv (tidyverse)
# This can cause issues with certain packages as its data frames are a 
# different format, you can get around this by changing them to data frames 
# after loading

taxtable <- read_csv("taxa_silva_species.csv")
seqtab <- read_csv("USV_counts_silva.csv")

# in this case we have a column to change to rownames
seqtab <- column_to_rownames(seqtab, "...1")
metadata <- read_tsv("metadatanojla.txt")
metadata <- data.frame(metadata)

# subset the metadata to only rows which it has in common with the sequences

metadata <- metadata[which(metadata[,1] %in% rownames(seqtab)),]
# and vice versa, to remove these two hangovers?
seqtab <- seqtab[metadata[,1],]
seqtab <- t(seqtab)

# in addition to the column_to_rownames function we can directly set rownames from a vector
rownames(metadata) <- metadata[,1]
# and then remove that column with negative addressing
metadata <- metadata[,-1]

# we want the patient id as there were errors with excel filling MRN numbers
# using string split, we can split the sample id by underscores and extract
# just the useful part for use as a patient id field
#test <- t(as.data.frame(strsplit(metadata$SampleID, "_")))
# and add it to our metadata frame
#metadata$ptid <- test[,2]

#metadata[which(metadata$Diagnosis == "normal"),6] <- "Normal"

# fix imported non standard characters
taxtable <- as.data.frame(taxtable)
rownames(taxtable) <- taxtable[,1]
taxtable <- taxtable[,-1]

# transpose the sequence table if necessary as this dumb function only works row-wise
# seqtab <- t(seqtab)

# create a copy of the taxonomy table with NA values replaced with best available taxonomy
#taxtable_nona <- as.matrix(taxtable)
#taxtable_nona <- turdpolisher(taxtable_nona, seqtab)
# rerun with outputseqtab
#taxtable_nona <- as.matrix(taxtable)
#taxtable_nona <- turdpolisher(as.matrix(outputtaxtable), as.matrix(outputseqtab))

# remove fully zero rows if any



# check that neither data frame is larger than the other or has duplicates
length(which(rownames(taxtable) %in% rownames(seqtab)))

# check that rownames are the same positionally
check <- data.frame(rownames(taxtable) == rownames(seqtab))

# subset by common rownames
#taxtable_nona <- taxtable_nona[which(rownames(taxtable_nona) %in% rownames(outputseqtab)),]
#outputseqtab <- outputseqtab[which(rownames(outputseqtab) %in% rownames(taxtable_nona)),]

# check that rownames are the same positionally
#check <- data.frame(rownames(taxtable_nona) == rownames(outputseqtab))

# remove zero sum rows

check <- data.frame(rowSums(seqtab))
seqtab <- seqtab[which(check$rowSums.seqtab. != 0),]
#taxtable_nona <- taxtable_nona[which(check$rowSums.outputseqtab. != 0),]

# remove eukaryotes and no kingdom assignment

seqtab <- seqtab[-which(taxtable$Kingdom == "Eukaryota"),]
taxtable <- taxtable[-which(taxtable$Kingdom == "Eukaryota"),]

seqtab <- seqtab[-which(is.na(taxtable$Kingdom) == TRUE),]
taxtable <- taxtable[-which(is.na(taxtable$Kingdom) == TRUE),]

# check that neither data frame is larger than the other or has duplicates
length(which(rownames(taxtable) %in% rownames(seqtab)))

# check that rownames are the same positionally
check <- data.frame(rownames(taxtable) == rownames(seqtab))

# overwrite original data with processed data
#seqtab <- outputseqtab
#taxtable <- taxtable_nona
#rownames(taxtable) <- taxtable[,1]
#taxtable <- taxtable[,-1]

# and tidy up after yourself
#rm(outputseqtab, outputtaxtable, taxtable_nona, test)

# store sequences for later use
#sequences <- as.data.frame(rownames(seqtab))

sequences <- data.frame(matrix(data = NA, nrow = nrow(taxtable), ncol = 2))
sequences[,2] <- rownames(taxtable)

# give sequences more manageable names
for(i in 1:nrow(seqtab)){
  rownames(seqtab)[i] <- paste0("ASV_", i)
}

# store corresponding ASV number for ease of use in subsetting
rownames(taxtable) <- rownames(seqtab)
sequences[,1] <- rownames(seqtab)
colnames(sequences) <- c("ASV", "Sequence")

#for (i in 1:dim(sequences)[1]) {
 # sequences[i,1] <- paste(">ASV", i, sep="_")
#}

# We won't be using them today but if you wanted to you could predict
# metagenomic function using these sequences

################################################
# IT IS UNLIKELY YOU NEED TO USE THE FOLLOWING #
################################################

# Metrics such as faith diversity and UniFrac distances utilise phylogenetic information.
# 16S phylogenetic relationships are falling out of favour because they usually
# only reflect an estimate of the whole genome phylogenetic relationship on a
# broad basis where you have lots of data, in a small study you may end up with
# strange, wrong, or unhelpful estimates

# read in phylogenetic tree if applicable
phylo <- ape::read.tree(file = "PhyloTree_Microbiomeanalyst.tre")

# replace tip label sequences with ASV number
for(i in 1:length(phylo[["tip.label"]])){
  phylo[["tip.label"]][i] <- sequences[which(sequences[,1] == phylo[["tip.label"]][i]),2]
}

###############################################

# tse structure:
# rowData - taxtable (samples as rows)
# colData - sample metadata (samples as rows)
# counts - sequence table (samples as columns)

tse <- TreeSummarizedExperiment(assays = list(counts = seqtab),
                                colData = metadata,
                                rowData = taxtable)
#rowTree = phylo)

# ensure mia doesn't fill its nappy, convert data frame to matrix
tse@assays@data@listData[["counts"]] <- as.matrix(tse@assays@data@listData[["counts"]])

# estimate shannon (measure of alpha diversity using number and distribution of taxa) and faith 
# (which uses phylogenetic relationships and will fail without a tree object)
# or the inverse Simpson metric, which takes into account distribution and relative
# abundance

indices <- c("shannon", "inverse_simpson")
names <- c("Shannon_diversity", "Inverse_Simpson")
tse <- estimateDiversity(tse, index = indices, name = names)

tse <- estimateRichness(tse, assay_name = abund_values, abund_values = "counts", index = "chao1")

# Alpha diversity is useful to say something quite broad about your data, each
# metric has strengths and weaknesses, and there are two fundamental approaches:
# 1. Assess alpha diversity in unfiltered data, which will be affected by
# rare taxa which must be kept in mind when interpreting results
# 2. Assess alpha diversity in fully processed, prevalence filtered data, which
# will likely only show you significant differences with a very stark signal in
# the data, ie IBD patients

# assess distribution, right skew is normal, adjust the bins to fit your data

shannon_hist <- ggplot(as.data.frame(colData(tse)), 
                       aes(x = Shannon_diversity)) + 
  geom_histogram(bins = 50, fill = "gray", color = "black") +
  labs(x = "Shannon index", y = "Sample frequency")

shannon_hist

tse <- transformCounts(tse, abund_values = "counts", method = "relabundance")

# assess alpha diversity by categorical metadata variables

shannon_box_disease <- ggplot(as.data.frame(colData(tse)),
                              aes(x = IBD_HS, 
                                  y = Shannon_diversity,
                                  fill = IBD_HS)) + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") +
  ylab(label = "Shannon Diversity Index") + xlab(label = NULL)


shannon_box_disease

kruskal.test(Shannon_diversity ~ IBD_HS, data = colData(tse))

# pval between IBD and HS unadjusted - 0.000001516

lme_data <- data.frame(tse@colData@listData)
rownames(lme_data) <- tse@colData@rownames
lme_data <- data.frame(lme_data)
lme_data$Subject <- as.factor(lme_data$Subject)
nullmod <- lmer(Shannon_diversity ~ (1|Subject), data=lme_data, REML = F)
fullmod <- lmer(Shannon_diversity ~ IBD_HS + (1|Subject), data=lme_data, REML = F)
tempres <- stats::anova(nullmod,fullmod)
tempres$`Pr(>Chisq)`[2]

# pval between IBD and HS by mixed effects model - 0.0000074372

is_box_disease <- ggplot(as.data.frame(colData(tse)),
                         aes(x = IBD_HS, 
                             y = Inverse_Simpson,
                             fill = IBD_HS)) + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") +
  ylab(label = "Inverse Simpson Index") + xlab(label = NULL)

is_box_disease

kruskal.test(Inverse_Simpson ~ IBD_HS, data = colData(tse))

# pval between IBD and HS unadjusted - 0.000006794

lme_data <- data.frame(tse@colData@listData)
rownames(lme_data) <- tse@colData@rownames
lme_data <- data.frame(lme_data)
lme_data$Subject <- as.factor(lme_data$Subject)
nullmod <- lmer(Inverse_Simpson ~ (1|Subject), data=lme_data, REML = F)
fullmod <- lmer(Inverse_Simpson ~ IBD_HS + (1|Subject), data=lme_data, REML = F)
tempres <- stats::anova(nullmod,fullmod)
tempres$`Pr(>Chisq)`[2]

# pval between IBD and HS by mixed effects model - 0.0003114062

# patients with cancer

shannon_box_cancer <- ggplot(as.data.frame(colData(tse)),
                                 aes(x = Cancer, 
                                     y = Shannon_diversity,
                                     fill = Cancer)) + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") +
  ylab(label = "Shannon Diversity Index") + xlab(label = NULL)

shannon_box_cancer

kruskal.test(Shannon_diversity ~ Cancer, data = colData(tse))

# pval between no and yes - 0.000216

lme_data <- data.frame(tse@colData@listData)
rownames(lme_data) <- tse@colData@rownames
lme_data <- data.frame(lme_data)
lme_data$Subject <- as.factor(lme_data$Subject)
nullmod <- lmer(Shannon_diversity ~ (1|Subject), data=lme_data, REML = F)
fullmod <- lmer(Shannon_diversity ~ Cancer + (1|Subject), data=lme_data, REML = F)
tempres <- stats::anova(nullmod,fullmod)
tempres$`Pr(>Chisq)`[2]

# adjusted pval 0.0261227

# cancer vs non, within ibd only

tempdf <- as.data.frame(colData(tse)[which(tse@colData@listData$IBD_HS == "IBD"),])
shannon_box_ibd_cancer <- ggplot(tempdf,
                             aes(x = Cancer, 
                                 y = Shannon_diversity,
                                 fill = Cancer)) + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") +
  ylab(label = "Shannon Diversity Index") + xlab(label = NULL)

shannon_box_ibd_cancer

kruskal.test(Shannon_diversity ~ Cancer, data = tempdf)

# pval between cancer pts and not within ibd subjects = 0.01545

nullmod <- lmer(Shannon_diversity ~ (1|Subject), data=tempdf, REML = F)
fullmod <- lmer(Shannon_diversity ~ Cancer + (1|Subject), data=tempdf, REML = F)
tempres <- stats::anova(nullmod,fullmod)
tempres$`Pr(>Chisq)`[2]

# pval between cancer patients and not within ibd subjects by mixed effects model = 0.1416008

# tests between cancers only

tempdf <- as.data.frame(colData(tse)[c(which(tse@colData@listData$Group == "Sporadic_Cancer") , which(tse@colData@listData$Group == "Colitis_Associated_Cancer")),])

#tse_cancersonly <- tse[,tse$Group %in% c("Sporadic_Cancer", "Colitis_Associated_Cancer")]

shannon_box_cacvs <- ggplot(tempdf,
                              aes(x = Group, 
                                  y = Shannon_diversity,
                                  fill = Group)) + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") +
  ylab(label = "Shannon Diversity Index") + xlab(label = NULL)

shannon_box_cacvs

kruskal.test(Shannon_diversity ~ Group, data = tempdf)

# pval between cancer samples - 0.7697, no mixed effects model as there are only singles

tempdf <- as.data.frame(colData(tse)[c(which(tse@colData@listData$Group == "Colitis_Associated_Cancer") , which(tse@colData@listData$Group == "CAC_normal_near_K")),])
tempdf$Group[which(tempdf$Group == "Colitis_Associated_Cancer")] <- "Colitis Associated Cancer"
tempdf$Group[which(tempdf$Group == "CAC_normal_near_K")] <- "Paired Normal Tissue"

#tse_cancerpts <- tse[,tse$Cancer %in% "yes"]
#tse_cacvsnear <- tse[,tse$Group %in% c("Colitis_Associated_Cancer", "CAC_normal_near_K")]
#tse_582 <- tse_genus[,tse_genus$]

shannon_box_cacvsnear <- ggplot(tempdf,
                            aes(x = Group, 
                                y = Shannon_diversity,
                                fill = Group)) + 
  geom_boxplot(lwd = 1, width = 0.5, outlier.shape = NA) + geom_jitter(width = 0.2, size = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1),
  axis.ticks = element_line(size = 1),
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 13),
  
  panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none") +
  ylab(label = "Shannon Diversity Index") + xlab(label = NULL)

shannon_box_cacvsnear

p <-ggpaired(data = tempdf,x = "Group", y = "Shannon_diversity", id = "Subject", fill = "Group", width = 0, point.size = 3, line.size = 1) + 
  geom_boxplot(data = tempdf, lwd = 1, width = 0.5, outlier.shape = NA, fill = c(scales::hue_pal()(2))) +
  xlab(label = NULL) + 
  ylab(label = "Shannon Diversity Index") +
  theme(axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1),
axis.ticks = element_line(size = 1),
axis.title = element_text(size = 13),
axis.text = element_text(size = 13))

p$layers[2:as.numeric(length(p$layers)+1)] <- p$layers[1:as.numeric(length(p$layers))]

p$layers[1] <- p$layers[5]
p$layers[5] <- NULL
p$layers[2] <- NULL

shannon_box_cacvsnear_paired <- p
shannon_box_cacvsnear_paired
rm(p)

kruskal.test(Shannon_diversity ~ Group, data = tempdf)

# pval between CAC and paired normal samples - 0.7494, no mixed effects model as there are only singles

#tse <- subsetByPrevalentTaxa(tse, detection = 0, prevalence = 0.15)
# check pcas with and without subsetting

tse <- transformCounts(tse, method = "rclr")

genus <- which(is.na(rowData(tse)$Genus))
species <- which(is.na(rowData(tse)$Species))
head(getTaxonomyLabels(tse[genus,]))

head(getTaxonomyLabels(tse[genus,], make_unique = FALSE))

head(getUniqueTaxa(tse, rank = "Genus"))

tse <- relAbundanceCounts(tse)
# if agglomeration is needed downstream, re-apply rclr transform
altExp(tse, "Family") <- agglomerateByRank(tse, rank = "Family")
altExp(tse, "Family") <- transformCounts(altExp(tse, "Family"), method = "rclr")

# if agglomeration is needed downstream, re-apply rclr transform
altExp(tse, "Genus") <- agglomerateByRank(tse, rank = "Genus")
altExp(tse, "Genus") <- transformCounts(altExp(tse, "Genus"), method = "rclr")

# if agglomeration is needed downstream, re-apply rclr transform
altExp(tse_cancersonly, "Genus") <- agglomerateByRank(tse_cancersonly, rank = "Genus")
altExp(tse_cancersonly, "Genus") <- transformCounts(altExp(tse_cancersonly, "Genus"), method = "rclr")

##### PERMANOVA

# Set seed for reproducibility
set.seed(12345)
# Perform dbRDA
dbrda <- dbrda(t(assay(tse,"rclr")) ~ Group, 
               data = colData(tse))

sppscores(dbrda) <- t(assay(tse,"rclr"))

# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term (here only 'Group') analyzed individually
                        method = "euclidean",
                        permutations = 9999)

tse_cancersonly <- tse[,tse$Group %in% c("Sporadic_Cancer", "Colitis_Associated_Cancer")]

# Set seed for reproducibility
set.seed(12345)
# Perform dbRDA
dbrda <- dbrda(t(assay(tse_cancersonly,"rclr")) ~ Group, 
               data = colData(tse_cancersonly))

sppscores(dbrda) <- t(assay(tse_cancersonly,"rclr"))

# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term (here only 'Group') analyzed individually
                        method = "euclidean",
                        permutations = 9999)
# extract coefficients
coef <- dbrda$CCA$v

# give nice names for plotting

rownames(coef) <- getTaxonomyLabels(tse)

abs_coef <- coef[rev(order(abs(coef))), , drop = FALSE]

# Get the taxa with biggest weights
top.coef <- head( coef[rev(order(abs(coef))), , drop = FALSE], 100)
# Sort weights in increasing order
top.coef <- top.coef[ order(top.coef), ]
# Get top names
top_names <- names(top.coef)[ order(abs(top.coef), decreasing = TRUE) ]

# 

library(scales)

# run aldex preliminarily to get significant taxa for highlighting

library(ALDEx2)
set.seed(123)
x <- aldex.clr(
  reads = assay(tse),
  conds = colData(tse)$IBD_HS,
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 1000,
  denom = "all",
  verbose = FALSE
)
x_tt <- aldex.ttest(
  x,
  paired.test = FALSE,
  verbose = FALSE)
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
aldex_out <- data.frame(x_tt, x_effect)

aldex_sig <- aldex_out
aldex_sig$ASV <- getTaxonomyLabels(tse)

aldex_sig <- aldex_sig %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than tt
  dplyr::select(ASV, we.eBH, wi.eBH, effect, overlap)

# ggplot visualisations

aldex_out$labels <- NA

# increased in non-ibd
aldex_out$labels[which(aldex_out$wi.eBH < 0.05 & aldex_out$effect > 0)] <- "Increased in non-IBD"

# increased in ibd
aldex_out$labels[which(aldex_out$wi.eBH < 0.05 & aldex_out$effect < 0)] <- "Increased in IBD"

aldex_out$labels[which(is.na(aldex_out$labels))] <- "Not differentially abundant"

aldex_out$size <- NA
aldex_out$size <- as.numeric(as.factor(aldex_out$labels))

aldex_out$size[which(aldex_out$size == 3)] <- 1
aldex_out$size[which(aldex_out$size != 1)] <- 2

  
# effect size plot

ggplot(data=aldex_out, aes(x=diff.win, y=diff.btw)) + 
  geom_segment(linetype = 2, size = 1.5, colour = "grey", 
               x = 0, y = 0, xend = max(aldex_out$diff.win), yend = max(aldex_out$diff.win)) + 
  geom_segment(linetype = 2, size = 1.5, colour = "grey", 
               x = 0, y = 0, xend = max(aldex_out$diff.win), yend = -max(aldex_out$diff.win)) + 
  ylim(-max(aldex_out$diff.btw), max(aldex_out$diff.btw)) +
  geom_point(aes(colour = labels, size = labels)) +
  scale_size_manual(values=c(2,2,1)) +
  xlab("Median Log Dispersion \n (Difference within groups)") +
  ylab("Median Difference between groups")

# ma plot

ggplot(data=aldex_out, aes(x=rab.all, y=diff.btw)) + 
  #geom_segment(linetype = 2, size = 1.5, colour = "grey", 
  #             x = 0, y = 0, xend = max(aldex_out$diff.win), yend = max(aldex_out$diff.win)) + 
  #geom_segment(linetype = 2, size = 1.5, colour = "grey", 
  #             x = 0, y = 0, xend = max(aldex_out$diff.win), yend = -max(aldex_out$diff.win)) + 
  #ylim(-max(aldex_out$diff.btw), max(aldex_out$diff.btw)) +
  geom_point(aes(colour = labels, size = labels)) +
  scale_size_manual(values=c(2,2,1)) +
  xlab("Median Log relative abundance") +
  ylab("Median Log Difference between groups")

# effect size vs q score volcano plot

ggplot(data=aldex_out, aes(x=effect, y=-log(we.eBH))) + 
  #geom_segment(linetype = 2, size = 1.5, colour = "grey", 
  #             x = 0, y = 0, xend = max(aldex_out$diff.win), yend = max(aldex_out$diff.win)) + 
  #geom_segment(linetype = 2, size = 1.5, colour = "grey", 
  #             x = 0, y = 0, xend = max(aldex_out$diff.win), yend = -max(aldex_out$diff.win)) + 
  xlim(-max(aldex_out$effect), max(aldex_out$effect)) +
  geom_point(aes(colour = labels, size = labels)) +
  scale_size_manual(values=c(2,2,1)) +
  xlab("Effect size") +
  ylab("Q score \n Negative log of p values")

# median difference between groups vs q score

ggplot(data=aldex_out, aes(x=diff.btw, y=-log(we.eBH))) + 
  #geom_segment(linetype = 2, size = 1.5, colour = "grey", 
  #             x = 0, y = 0, xend = max(aldex_out$diff.win), yend = max(aldex_out$diff.win)) + 
  #geom_segment(linetype = 2, size = 1.5, colour = "grey", 
  #             x = 0, y = 0, xend = max(aldex_out$diff.win), yend = -max(aldex_out$diff.win)) + 
  xlim(-max(aldex_out$diff.btw), max(aldex_out$diff.btw)) +
  geom_point(aes(colour = labels, size = labels)) +
  scale_size_manual(values=c(2,2,1)) +
  xlab("Median log difference between groups") +
  ylab("Q score \n Negative log of p values")

# to do in future - tease out which ASVs 95%ci does not cross 0
# loop for effect > 0, check min > 0 and colour accordingly, etc
# none in this dataset but maybe others?

library(scales)

fillvar <- names(top.coef) %in% aldex_sig[,1]
fillvar[which(fillvar != TRUE)] <- "gray"
fillvar[which(fillvar == TRUE)] <- hue_pal()(1)

ggplot(data.frame(x = top.coef,
                  y = factor(names(top.coef),
                             unique(names(top.coef)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity", fill = fillvar) +
  labs(x="",y="",title="Top ASVs associated with PERMANOVA model") +
  theme_bw()

#tse_cancersonly <- transformCounts(tse_cancersonly, method = "rclr")  

# Gets clr table
clr_assay <- assays(tse)$clr

# Transposes it to get taxa to columns
clr_assay <- t(clr_assay)

clr_pca <- prcomp(t(tse_cancersonly@assays@data@listData[["clr"]]))
df_out <- as.data.frame(clr_pca$x)

# construct df of all pc combinations, warning, if you have a lot this might die
allpcs <- as.data.frame(t(combn(colnames(df_out), 2)))

df_out$Disease <- as.factor(colData(tse_cancersonly)$Group)

pcx <- "PC1"
pcy <- "PC2"

p<-ggplot(df_out,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=Disease))
p<-p+geom_point()
p

# get convex hulls

hull <- df_out %>%
  group_by(Disease) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p + aes(fill = factor(Disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(x = paste(pcx, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcx, "PC")[[1]][2])][2])*100, "% variance explained"),
       y = paste(pcy, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcy, "PC")[[1]][2])][2])*100, "% variance explained"))

mult=min(max(clr_pca$x[,1])/max(clr_pca$rotation[,1]),max(clr_pca$x[,1+1])/max(clr_pca$rotation[,1+1]))
loadings <- data.frame(clr_pca$rotation)
#rownames(loadings) <- paste()

x <- sym(paste0("PC", "1"))
y <- sym(paste0("PC", "2"))

p <- p + aes(fill = factor(Disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(x = paste(pcx, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcx, "PC")[[1]][2])][2])*100, "% variance explained"),
       y = paste(pcy, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcy, "PC")[[1]][2])][2])*100, "% variance explained")) + 
  geom_segment(inherit.aes = FALSE, aes(x = 0, y = 0,xend=mult*loadings[6,1],yend=mult*loadings[6,2]),arrow=arrow(length=unit(.3,"lines")),color="gray60",size=.4)
p

# add known species to geom_label?

library(ggrepel)

tempNames <- data.frame(tse_cancersonly@rowRanges@elementMetadata@listData[["Genus"]][which(tse_cancersonly@rowRanges@elementMetadata@listData[["Genus"]] == "Fusobacterium")])
colnames(tempNames)[1] <- "names"
tempNames$names[which(is.na(tempNames$names))] <- " "

p + geom_text(aes(label = tempNames$names, x = mult*loadings[which(tse_cancersonly@rowRanges@elementMetadata@listData[["Genus"]] == "Fusobacterium"),"PC1"], y = mult*loadings[which(tse_cancersonly@rowRanges@elementMetadata@listData[["Genus"]] == "Fusobacterium"),"PC2"] +1))

tse_phylum <- agglomerateByRank(tse, rank ="Phylum", onRankOnly=TRUE)
top_taxa <- getTopTaxa(tse_phylum,top = 5, abund_values = "relabundance")

phylum_renamed <- lapply(rowData(tse)$Phylum,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})
#rowData(tse)$Phylum <- as.character(phylum_renamed)

boxplotdf <- data.frame("group" = colData(tse_cacvsnear_genus)$Group, relabund = getAbundanceFeature(tse_cacvsnear_genus, feature_id = "Fusobacterium", abund_values = "relabundance"))
boxplotdf$relabund <- boxplotdf$relabund * 100

box_fuso <- ggplot(boxplotdf,
                   aes(x = group, 
                       y = relabund,
                       fill = group)) + 
  geom_boxplot() + geom_point()

set.seed(123)
x <- aldex.clr(
  reads = assay(tse_genus),
  conds = colData(tse_genus)$IBD_HS, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 1000, 
  denom = "all",
  verbose = FALSE
)

x_tt <- aldex.ttest(
  x, 
  paired.test = FALSE, 
  verbose = FALSE)

x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)

aldex_out <- data.frame(x_tt, x_effect)

par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.05
)
aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.05
)

rownames_to_column(aldex_out, "genus") %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than tt
  select(genus, we.eBH, wi.eBH, effect, overlap)

# Identify structural zeros from ANCOM-II

# Identify structural zeros
get_struc_zero = function(tse, tax_level, assay_name,
                          alt = FALSE, group, neg_lb) {
  if (alt) {
    tse_alt = SingleCellExperiment::altExp(tse, tax_level)
    feature_table = SummarizedExperiment::assay(tse_alt, assay_name)
    meta_data = SummarizedExperiment::colData(tse_alt)
    tax_name = rownames(tse_alt)
  } else {
    feature_table = SummarizedExperiment::assay(tse, assay_name)
    meta_data = SummarizedExperiment::colData(tse)
    tax_name = rownames(tse)
  }
  group_data = factor(meta_data[, group])
  present_table = as.matrix(feature_table)
  present_table[is.na(present_table)] = 0
  present_table[present_table != 0] = 1
  n_tax = nrow(feature_table)
  n_group = nlevels(group_data)
  
  p_hat = matrix(NA, nrow = n_tax, ncol = n_group)
  rownames(p_hat) = rownames(feature_table)
  colnames(p_hat) = levels(group_data)
  for (i in seq_len(n_tax)) {
    p_hat[i, ] = tapply(present_table[i, ], group_data,
                        function(x) mean(x, na.rm = TRUE))
  }
  
  samp_size = matrix(NA, nrow = n_tax, ncol = n_group)
  rownames(samp_size) = rownames(feature_table)
  colnames(samp_size) = levels(group_data)
  for (i in seq_len(n_tax)) {
    samp_size[i, ] = tapply(as.matrix(feature_table)[i, ], group_data,
                            function(x) length(x[!is.na(x)]))
  }
  
  p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)
  
  output = (p_hat == 0)
  # Shall we classify a taxon as a structural zero by its negative lower bound?
  if (neg_lb) output[p_hat_lo <= 0] = TRUE
  
  output = cbind(tax_name, output)
  colnames(output) = c("taxon",
                       paste0("structural_zero (", group,
                              " = ", colnames(output)[-1], ")"))
  output = data.frame(output, check.names = FALSE, row.names = NULL)
  output[, -1] = apply(output[, -1], 2, as.logical)
  return(output)
}













seqtab_roberts <- read_csv(file = "C:/Users/eatmo/Desktop/CAC/Roberts/USV_counts_silva.csv")
seqtab_roberts <- data.frame(seqtab_roberts)
rownames(seqtab_roberts) <- seqtab_roberts[,1]
seqtab_roberts <- seqtab_roberts[,-1]
seqtab_roberts <- t(seqtab_roberts)
metadata_roberts <- read_csv(file = "C:/Users/eatmo/Desktop/CAC/Roberts/tempmetadata.csv")
metadata_roberts <- data.frame(metadata_roberts)
rownames(metadata_roberts) <- metadata_roberts[,1]
metadata_roberts <- metadata_roberts[,-1]
taxtable_roberts <- read_csv(file = "C:/Users/eatmo/Desktop/CAC/Roberts/taxa_silva_species.csv")

seqtab_roberts <- seqtab_roberts[,-146]
seqtab_roberts <- seqtab_roberts[,-146]

sequences_roberts <- taxtable[,1]
taxtable_roberts <- taxtable_roberts[,-1]

taxtable_nona <- as.matrix(taxtable_roberts)

taxtable_nona <- turdpolisher(taxtable_nona, seqtab_roberts)

which(rownames(metadata_roberts) %in% colnames(seqtab_roberts))

tse_roberts <- TreeSummarizedExperiment(assays = list(counts = seqtab_roberts),
                                colData = metadata_roberts,
                                rowData = taxtable_nona)

tse_roberts_genus <- agglomerateByRank(tse_roberts, rank ="Genus", onRankOnly=TRUE)
tse_roberts_genus <- transformCounts(tse_roberts_genus, abund_values = "counts", method = "relabundance")
tse_roberts_genus <- subsetByPrevalentTaxa(tse_roberts_genus, detection = 0, prevalence = 0.15)



tse_roberts_cacvscrc <- tse_roberts_genus[,tse_roberts_genus$allCACvsCRC %in% c("sCRC", "CAC")]

tse_roberts_onlyt <- tse_roberts_cacvscrc[,tse_roberts_cacvscrc$TvsN %in% "T"]

tse_roberts_cacvscrc <- transformCounts(tse_roberts_cacvscrc, method = "clr", pseudocount = 1)
tse_roberts_onlyt <- transformCounts(tse_roberts_onlyt, method = "clr", pseudocount = 1)

clr_pca <- prcomp(t(tse_roberts_onlyt@assays@data@listData[["clr"]]))
df_out <- as.data.frame(clr_pca$x)

# construct df of all pc combinations, warning, if you have a lot this might die
allpcs <- as.data.frame(t(combn(colnames(df_out), 2)))

df_out$Group <- as.factor(colData(tse_roberts_onlyt)$Type)

pcx <- "PC1"
pcy <- "PC2"

p<-ggplot(df_out,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=Group))
p<-p+geom_point()
#p

# get convex hulls

hull <- df_out %>%
  group_by(Group) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p + aes(fill = Group) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(x = paste(pcx, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcx, "PC")[[1]][2])][2])*100, "% variance explained"),
       y = paste(pcy, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcy, "PC")[[1]][2])][2])*100, "% variance explained"))

mult=min(max(clr_pca$x[,1])/max(clr_pca$rotation[,1]),max(clr_pca$x[,1+1])/max(clr_pca$rotation[,1+1]))
loadings <- data.frame(clr_pca$rotation)
#rownames(loadings) <- paste()

x <- sym(paste0("PC", "1"))
y <- sym(paste0("PC", "2"))

p + aes(fill = factor(Disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(x = paste(pcx, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcx, "PC")[[1]][2])][2])*100, "% variance explained"),
       y = paste(pcy, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcy, "PC")[[1]][2])][2])*100, "% variance explained")) + 
  geom_segment(data=loadings[which(tse_cancersonly@rowRanges@elementMetadata@listData[["Genus"]] == "Fusobacterium"),], inherit.aes = FALSE,aes(x = 0, y = 0,xend=mult*!!x,yend=mult*!!y),arrow=arrow(length=unit(.3,"lines")),color="gray60",size=.4)
p

clr_pca <- prcomp(t(tse_roberts_onlyt@assays@data@listData[["clr"]]))
df_out <- as.data.frame(clr_pca$x)

# construct df of all pc combinations, warning, if you have a lot this might die
allpcs <- as.data.frame(t(combn(colnames(df_out), 2)))

df_out$Group <- as.factor(colData(tse_roberts_onlyt)$allCACvsCRC)

pcx <- "PC2"
pcy <- "PC3"

p<-ggplot(df_out,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=Group))
p<-p+geom_point()
#p

# get convex hulls

hull <- df_out %>%
  group_by(Group) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p + aes(fill = factor(Group)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(x = paste(pcx, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcx, "PC")[[1]][2])][2])*100, "% variance explained"),
       y = paste(pcy, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcy, "PC")[[1]][2])][2])*100, "% variance explained"))

mult=min(max(clr_pca$x[,1])/max(clr_pca$rotation[,1]),max(clr_pca$x[,1+1])/max(clr_pca$rotation[,1+1]))
loadings <- data.frame(clr_pca$rotation)
#rownames(loadings) <- paste()

x <- sym(paste0("PC", "1"))
y <- sym(paste0("PC", "2"))

p + aes(fill = factor(Disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(x = paste(pcx, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcx, "PC")[[1]][2])][2])*100, "% variance explained"),
       y = paste(pcy, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcy, "PC")[[1]][2])][2])*100, "% variance explained")) + 
  geom_segment(data=loadings[which(tse_cancersonly@rowRanges@elementMetadata@listData[["Genus"]] == "Fusobacterium"),], inherit.aes = FALSE,aes(x = 0, y = 0,xend=mult*!!x,yend=mult*!!y),arrow=arrow(length=unit(.3,"lines")),color="gray60",size=.4)
p

tse_roberts_578 <- tse_roberts_genus[,tse_roberts_genus$JLA %in% "578"]

c("Terminal Ileum","Terminal Ileum", "Caecum","Caecum", "Ascending Colon","Ascending Colon", "Hepatic Flexure","Hepatic Flexure", "Mid Transverse Colon", "Mid Transverse Colon", "Splenic Flexure","Splenic Flexure", "Proximal Descending Colon","Proximal Descending Colon", "Distal Descending Colon","Distal Descending Colon", "Sigmoid Colon","Sigmoid Colon", "Recto-Sigmoid","Recto-Sigmoid")

tse_roberts_578@colData@listData[["locations"]] <- c("Ascending Colon", "Caecum", "Distal Descending Colon", "Hepatic Flexure", "Mid Transverse Colon", "Proximal Descending Colon", "Recto-Sigmoid", "Sigmoid Colon", "Splenic Flexure", "Terminal Ileum")

boxplotdf <- data.frame("group" = colData(tse_roberts_578)$locations, relabund = getAbundanceFeature(tse_roberts_578, feature_id = "Fusobacterium", abund_values = "relabundance"))
boxplotdf$relabund <- boxplotdf$relabund * 100

boxplotdf <- boxplotdf[match(unique(c("Terminal Ileum","Terminal Ileum", "Caecum","Caecum", "Ascending Colon","Ascending Colon", "Hepatic Flexure","Hepatic Flexure", "Mid Transverse Colon", "Mid Transverse Colon", "Splenic Flexure","Splenic Flexure", "Proximal Descending Colon","Proximal Descending Colon", "Distal Descending Colon","Distal Descending Colon", "Sigmoid Colon","Sigmoid Colon", "Recto-Sigmoid","Recto-Sigmoid")), boxplotdf$group),]
boxplotdf2 <- data.frame(matrix(data = 0, nrow = nrow(boxplotdf)*2, ncol = 1))
boxplotdf2$group <- c("Terminal Ileum","Terminal Ileum", "Caecum","Caecum", "Ascending Colon","Ascending Colon", "Hepatic Flexure","Hepatic Flexure", "Mid Transverse Colon", "Mid Transverse Colon", "Splenic Flexure","Splenic Flexure", "Proximal Descending Colon","Proximal Descending Colon", "Distal Descending Colon","Distal Descending Colon", "Sigmoid Colon","Sigmoid Colon", "Recto-Sigmoid","Recto-Sigmoid")
boxplotdf2$Bacteria <- rep(c("Fusobacterium", "Other"), 10)
boxplotdf2$sample <- rep(1:10, each = 2)
colnames(boxplotdf2)[1] <- "value"
boxplotdf2$value <- rep(getAbundanceFeature(tse_roberts_578, feature_id = "Fusobacterium", abund_values = "relabundance"), each = 2)
boxplotdf2$value <- boxplotdf2$value*100

for(i in 1:(nrow(boxplotdf2)/2)){
  boxplotdf2[i*2,1] <- 100 - boxplotdf2[(i*2)-1,1]
}

ggplot(boxplotdf2, aes(x=as.numeric(factor(group)), y=value, fill = Bacteria)) + scale_x_continuous(breaks = unique(boxplotdf2$sample), labels = unique(boxplotdf2$group)) + xlab("Anatomical Region") + ylab("Relative Abundance (%)") +
  geom_area() + geom_vline(xintercept = 2, linetype="dotted", color = "red", size=1.5)

box_fuso <- ggplot(boxplotdf,
                   aes(x = group, 
                       y = relabund,
                       fill = group)) + 
  geom_boxplot() + geom_point()



tse_roberts_582 <- tse_roberts_genus[,tse_roberts_genus$JLA %in% "582"]

check <- matrix(nrow = 2, ncol = 14, data = unlist(strsplit(tse_roberts_582@colData@rownames, "582")))
check <- c("Ascending Colon", "Caecum", "Descending Colon", "Distal Descending Colon", "Hepatic Flexure", "Distal Sigmoid Colon", "Sigmoid Colon", "Proximal Sigmoid Colon", "Sigmoid Proximal to Tumour", "Splenic Flexure", "Terminal Ileum", "Transverse Colon", "Distal Transverse Colon", "Tumour")
tse_roberts_582@colData@listData[["locations"]] <- check

boxplotdf <- data.frame("group" = colData(tse_roberts_582)$locations, relabund = getAbundanceFeature(tse_roberts_582, feature_id = "Fusobacterium", abund_values = "relabundance"))

boxplotdf <- boxplotdf[match(c("Terminal Ileum", "Caecum", "Ascending Colon", "Hepatic Flexure", "Transverse Colon", "Distal Transverse Colon", "Splenic Flexure", "Descending Colon", "Distal Descending Colon", "Proximal Sigmoid Colon", "Sigmoid Proximal to Tumour", "Tumour", "Sigmoid Colon", "Distal Sigmoid Colon"), boxplotdf$group),]
boxplotdf2 <- data.frame(matrix(data = 0, nrow = 28, ncol = 1))
boxplotdf2$group <- rep(c("Terminal Ileum", "Caecum", "Ascending Colon", "Hepatic Flexure", "Transverse Colon", "Distal Transverse Colon", "Splenic Flexure", "Descending Colon", "Distal Descending Colon", "Proximal Sigmoid Colon", "Sigmoid Proximal to Tumour", "Tumour", "Sigmoid Colon", "Distal Sigmoid Colon"), each=2)
boxplotdf2$Bacteria <- rep(c("Fusobacterium", "Other"), 14)
boxplotdf2$sample <- rep(1:14, each = 2)
colnames(boxplotdf2)[1] <- "value"
boxplotdf2$value <- rep(getAbundanceFeature(tse_roberts_582, feature_id = "Streptococcus", abund_values = "relabundance"), each = 2)
boxplotdf2$value <- boxplotdf2$value*100

for(i in 1:(nrow(boxplotdf2)/2)){
  boxplotdf2[i*2,1] <- 100 - boxplotdf2[(i*2)-1,1]
}

ggplot(boxplotdf2, aes(x=as.numeric(factor(group)), y=value, fill = Bacteria)) + scale_x_continuous(breaks = unique(boxplotdf2$sample), labels = unique(boxplotdf2$group)) + xlab("Anatomical Region") + ylab("Relative Abundance (%)") +
  geom_area() + geom_vline(xintercept = 12, linetype="dotted", color = "red", size=1.5)


tse_roberts_636 <- tse_roberts_genus[,tse_roberts_genus$JLA %in% "636"]
tse_roberts_636 <- tse_roberts_636[,tse_roberts_636$SampleType %in% "tissue"]

check <- matrix(nrow = 2, ncol = 12, data = unlist(strsplit(tse_roberts_636@colData@rownames, "636")))
check[2,]
check <- c("Ascending Colon", "Caecum", "Descending Colon", "Distal Sigmoid Colon", "Hepatic Flexure", "Mid Transverse Colon", "Proximal Descending Colon", "Proximal Transverse Colon", "Sigmoid Colon", "Splenic Flexure", "Terminal Ileum", "Tumour")
tse_roberts_636@colData@listData[["locations"]] <- check

boxplotdf <- data.frame("group" = colData(tse_roberts_636)$locations, relabund = getAbundanceFeature(tse_roberts_636, feature_id = "Fusobacterium", abund_values = "relabundance"))

boxplotdf <- boxplotdf[match(c("Terminal Ileum", "Caecum", "Ascending Colon", "Hepatic Flexure", "Proximal Transverse Colon", "Mid Transverse Colon", "Splenic Flexure", "Proximal Descending Colon", "Descending Colon", "Sigmoid Colon", "Distal Sigmoid Colon", "Tumour"), boxplotdf$group),]
boxplotdf2 <- data.frame(matrix(data = 0, nrow = 24, ncol = 1))
boxplotdf2$group <- rep(c("Terminal Ileum", "Caecum", "Ascending Colon", "Hepatic Flexure", "Proximal Transverse Colon", "Mid Transverse Colon", "Splenic Flexure", "Proximal Descending Colon", "Descending Colon", "Sigmoid Colon", "Distal Sigmoid Colon", "Tumour"), each=2)
boxplotdf2$Bacteria <- rep(c("Fusobacterium", "Other"), 12)
boxplotdf2$sample <- rep(1:12, each = 2)
colnames(boxplotdf2)[1] <- "value"
boxplotdf2$value <- rep(getAbundanceFeature(tse_roberts_636, feature_id = "Fusobacterium", abund_values = "relabundance"), each = 2)
boxplotdf2$value <- boxplotdf2$value*100

for(i in 1:(nrow(boxplotdf2)/2)){
  boxplotdf2[i*2,1] <- 100 - boxplotdf2[(i*2)-1,1]
}

ggplot(boxplotdf2, aes(x=as.numeric(factor(group)), y=value, fill = Bacteria)) + scale_x_continuous(breaks = unique(boxplotdf2$sample), labels = unique(boxplotdf2$group)) + xlab("Anatomical Region") + ylab("Relative Abundance (%)") +
  geom_area() + geom_vline(xintercept = 12, linetype="dotted", color = "red", size=1.5)

tse_roberts_ibd <- tse_roberts_genus[,tse_roberts_genus$Group %in% c("UC", "CD", "CD?")]
tse_roberts_cac_t <- tse_roberts_ibd[,tse_roberts_ibd$TvsN %in% "T"]

ibd_fuso <- data.frame(getAbundanceFeature(tse_roberts_ibd, feature_id = "Fusobacterium", abund_values = "relabundance"))

length(which(ibd_fuso != 0))

cac_fuso <- data.frame(getAbundanceFeature(tse_roberts_cac_t, feature_id = "Fusobacterium", abund_values = "relabundance"))

tse_roberts_scrc <- tse_roberts_genus[,tse_roberts_genus$Group %in% "sCRC"]
tse_roberts_scrc_t <- tse_roberts_scrc[,tse_roberts_scrc$TvsN %in% "T"]


crc_fuso <- data.frame(getAbundanceFeature(tse_roberts_scrc, feature_id = "Fusobacterium", abund_values = "relabundance"))

length(which(crc_fuso != 0))

crc_fuso_t <- data.frame(getAbundanceFeature(tse_roberts_scrc_t, feature_id = "Fusobacterium", abund_values = "relabundance"))

length(which(crc_fuso_t != 0))

set.seed(123)
x <- aldex.clr(
  reads = assay(tse_roberts_scrc),
  conds = colData(tse_roberts_scrc)$TvsN, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 1000, 
  denom = "all",
  verbose = FALSE
)
x_tt <- aldex.ttest(
  x, 
  paired.test = FALSE, 
  verbose = FALSE)

x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)

par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.05
)
aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.05
)

tse_roberts_onlyt <- tse_roberts_genus[,tse_roberts_genus$TvsN %in% "T"]
# remove 409, not a cac

tse_roberts_onlyt <- tse_roberts_onlyt[,-10]

set.seed(123)
x <- aldex.clr(
  reads = assay(tse_roberts_onlyt),
  conds = colData(tse_roberts_onlyt)$allCACvsCRC, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 1000, 
  denom = "all",
  verbose = FALSE
)
x_tt <- aldex.ttest(
  x, 
  paired.test = FALSE, 
  verbose = FALSE)

x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)

par(mfrow = c(1, 2))
aldex.plot(
  aldex_out, 
  type = "MA", 
  test = "welch", 
  xlab = "Log-ratio abundance",
  ylab = "Difference",
  cutoff = 0.05
)
aldex.plot(
  aldex_out, 
  type = "MW", 
  test = "welch",
  xlab = "Dispersion",
  ylab = "Difference",
  cutoff = 0.05
)



permanova_ibdvsnon <- vegan::adonis(t(tse_roberts_genus@assays@data$relabundance) ~ metadata_roberts$Group,
                                  data = colData(tse_roberts_genus),
                                  permutations = 9999)





turdpolisher <- function(taxatable, seqtab){
  exit <- 0
  addr <- NULL
  
  for (i in seq_along(taxatable[,1])){
    if (is.na(taxatable[i,1]) == TRUE){
      print(paste0("Fully NA row detected: #",i))
      addr <- c(addr, i)
      exit <- 1
    }
  }
  
  if (exit != 1){
    
    inc <- 1
    loopr <- 1
    
    for (i in seq_along(taxatable[,])){
      if (inc == (ncol(taxatable))){
        if (loopr != nrow(taxatable)){
          loopr <- loopr + 1
          #print(paste("loopr",loopr))
          tmpname <- NULL}
        #print("HELLO")
        inc <- 1 
      }
      else{
        inc <- inc + 1
        #print(paste("inc",inc))
      }
      if (is.na(taxatable[loopr,inc]) == TRUE){
        #print("did this")
        if (is.null(tmpname) == TRUE){tmpname <- paste0("Unclassified_", taxatable[loopr,(inc-1)])}
        
        taxatable[loopr,inc] <- tmpname
      }
    }
    return(taxatable)
  }
  else{
    print("Remove NA rows and try again")
    outputseqtab <<- seqtab[-c(addr),]
    outputtaxtable <<- taxtable[-c(addr),]
    print("Created new sequence/taxonomy objects without NA rows as outputseqtab/taxtable")
  }
}