turdpolisher <- function(taxatable, seqtab){
  exit <- 0
  addr <- NULL
  tmpname <- NULL
  
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

boxedIn <- function(input_tse, div_ind, fullmod, var, x_lab, y_lab) {
  boxplot <- data.frame(c(input_tse@colData@rownames, input_tse@colData@rownames))
  colnames(boxplot)[1] <- "names"
  boxplot$counts <- c(input_tse@colData@listData[[div_ind]], fitted(fullmod))
  boxplot$type <- c(rep("Raw data", ncol(input_tse)), rep("Fitted values", ncol(input_tse)))
  boxplot$grp <- unlist(rep(input_tse@colData[which(names(input_tse@colData@listData) == var)], 2))
  ggplot(boxplot,
         aes(x = grp,
             y = counts,
             fill = grp)) + labs(x = x_lab, y = y_lab, title = c(paste(div_ind, "fitted vs raw"))) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~type)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("mia")
BiocManager::install("miaViz")
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# load packages for processing
# we'll be using the miaverse suite, it's harder to get used to but very flexible
# all sorts of packages will be missing, add installs/loads as required

library(tidyverse)
library(dplyr)
library(mia)
library(miaViz)
library(scater)
library(indicspecies)
library(ggpubr)
library(sp)
library(pracma)
library(pairwiseAdonis)



# change path to wherever your files are, format might be different on mac
path_to_files <- "C:/users/eatmo/Desktop/CAC/Roberts/"

setwd(path_to_files)

# load in the components to construct our Tree Summarised Experiment object

taxtable <- read.csv("taxa_silva.csv")
seqtab <- read.csv("USV_counts_silva.csv")
metadata <- read.csv("tempmetadata.csv")

metadata$Serial <- 0
metadata[which(metadata$JLA == "578"),13] <- "Yes"
metadata[which(metadata$JLA == "582"),13] <- "Yes"
metadata[which(metadata$JLA == "636"),13] <- "Yes"
metadata[which(metadata$JLA == "646"),13] <- "Yes"

metadata[which(metadata$Serial == "0"),13] <- "No"


taxtable_abs <- read.csv("taxa_silva.csv")
seqtab_abs <- read.csv("absolute_usvs_silva.csv")
metadata_abs <- read.csv("tempmetadata_abs.csv")



metadata_abs[which(metadata_abs$Group == "CD?"),3] <- "CD"

#test <- t(as.data.frame(strsplit(metadata$SampleID, "_")))
#metadata$ptid <- test[,2]

#metadata[which(metadata$Diagnosis == "normal"),7] <- "Normal"

metadata[which(metadata$Group == "CD?"),3] <- "CD"

# fix imported non standard characters
taxtable <- column_to_rownames(taxtable, "X")
taxtable_abs <- column_to_rownames(taxtable_abs, "X")
seqtab <- column_to_rownames(seqtab, "X")
seqtab_abs <- column_to_rownames(seqtab_abs, "X")

# arrange seqtab columns in same order as metadata samples so it doesn't do wacky stuff



seqtab2 <- seqtab_abs
seqtab <- seqtab[c(metadata$Name),]
#seqtab_abs <- t(seqtab_abs)
seqtab_abs <- seqtab_abs[,c(metadata_abs$Name)]

# create a copy of the taxonomy table with NA values replaced with best available taxonomy
taxtable_nona <- as.matrix(taxtable)
taxtable_nona <- turdpolisher(taxtable_nona, seqtab)
# rerun with outputseqtab if generated
taxtable_nona <- as.matrix(taxtable)
taxtable_nona <- turdpolisher(as.matrix(outputtaxtable), as.matrix(outputseqtab))

# store sequences for later use
seqtab <- seqtab[rownames(seqtab) %in% metadata[,1],]

sequences <- as.data.frame(colnames(seqtab))

seqtab <- t(seqtab)

# give sequences more manageable names
for(i in 1:nrow(seqtab)){
  rownames(seqtab)[i] <- paste0("ASV_", i)
}

# give sequences more manageable names
for(i in 1:nrow(seqtab_abs)){
  rownames(seqtab_abs)[i] <- paste0("ASV_", i)
}

# store corresponding ASV number for ease of use in subsetting
rownames(taxtable_nona) <- rownames(seqtab)
rownames(taxtable_abs) <- rownames(seqtab_abs)
sequences[,2] <- rownames(seqtab)

# read in phylogenetic tree if applicable
#phylo <- ape::read.tree(file = "PhyloTree_Microbiomeanalyst.tre")

# replace tip label sequences with ASV number
#for(i in 1:length(phylo[["tip.label"]])){
#  phylo[["tip.label"]][i] <- sequences[which(sequences[,1] == phylo[["tip.label"]][i]),2]
#}

# tse structure:
# rowData - taxtable (samples as rows)
# colData - sample metadata (samples as rows)
# counts - sequence table (samples as columns)

# censor samples?
which(colnames(seqtab) == "Natasha75")

seqtabfilt <- seqtab[,-66]

which(colnames(seqtabfilt) == "Natasha96")

seqtabfilt <- seqtabfilt[,-83]

which(metadata$X.NAME == "Natasha75")

metadatafilt <- metadata[-27,]

which(metadatafilt$X.NAME == "Natasha96")



tse <- TreeSummarizedExperiment(assays = list(counts = seqtab),
                                colData = metadata,
                                rowData = taxtable_nona)

tse_abs <- TreeSummarizedExperiment(assays = list(counts = seqtab_abs),
                                colData = metadata_abs,
                                rowData = taxtable_abs)

#rowTree = phylo)

# remove Natasha75 and Natasha96

which(colnames(tse) == "Natasha75")

tse <- tse[,-27]

which(colnames(tse) == "Natasha96")

tse <- tse[,-83]

# ensure mia doesn't fill its nappy, convert data frame to matrix
tse@assays@data@listData[["counts"]] <- as.matrix(tse@assays@data@listData[["counts"]])

# estimate shannon (measure of alpha diversity using number and distribution of taxa) and faith 
# (which uses phylogenetic relationships and will fail without a tree object)
indices <- c("shannon", "inverse_simpson")
names <- c("Shannon_diversity", "Inverse_Simpson_Index")
tse <- estimateDiversity(tse, index = indices, name = names)


# assess distribution, right skew is normal

shannon_hist <- ggplot(as.data.frame(colData(tse)), 
                       aes(x = Shannon_diversity)) + 
  geom_histogram(bins = 25, fill = "gray", color = "black") +
  labs(x = "Shannon index", y = "Sample frequency")

shannon_hist

#analysis subsets:
  
#  CAC vs CRC, - stool, - NAs

tse_cac_vs_crc <- tse[,-which(is.na(tse$allCACvsCRC))]
tse_cac_vs_crc <- tse_cac_vs_crc[,-which(tse_cac_vs_crc$SampleType == "stool"))]
tse_cac_vs_crc <- estimateDiversity(tse_cac_vs_crc, index = indices, name = names)

shannon_box_loc <- ggplot(as.data.frame(colData(tse_cac_vs_crc)),
                          aes(x = Group, 
                              y = Shannon_diversity,
                              fill = Group)) + xlab("Disease Status") +
  ylab("Shannon Diversity") +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  guides(fill = "none") +
  theme(title = element_text(size = 12)) # makes titles smaller

shannon_box_loc

FSA::dunnTest(tse_cac_vs_crc@colData@listData[["Shannon_diversity"]]~tse_cac_vs_crc@colData@listData[["Group"]])

# significant differences but many multiple samples, do mixed effects models

tempdf <- tse_cac_vs_crc@colData@listData
nullmod <- lmer(tempdf$Shannon_diversity ~ 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
fullmod <- lmer(tempdf$Shannon_diversity ~ Group + 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
test <- anova(nullmod,fullmod)

# visualise fitted vs raw

boxplot <- data.frame(c(tse_cac_vs_crc@colData@rownames, tse_cac_vs_crc@colData@rownames))
colnames(boxplot)[1] <- "names"
boxplot$counts <- c(tse_cac_vs_crc$Shannon_diversity, fitted(fullmod))
boxplot$type <- c(rep("Raw data", ncol(tse_cac_vs_crc)), rep("Fitted values", ncol(tse_cac_vs_crc)))
boxplot$grp <- unlist(rep(tse_cac_vs_crc@colData[which(names(tse_cac_vs_crc@colData@listData) == "Group")], 2))
ggplot(boxplot,
       aes(x = grp,
           y = counts,
           fill = grp)) + labs(title = "Shannon diversity fitted vs raw") +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~type)

boxedIn(tse_cac_vs_crc, "Shannon_diversity", fullmod, "Group")

# only one cd patient, 2 very low diversity outliers in the scrc group

kruskal.test(Shannon_diversity ~ allCACvsCRC, data = colData(tse_cac_vs_crc))
# very significant, check lme

tempdf <- tse_cac_vs_crc@colData@listData
tempdf$JLA <- factor(tempdf$JLA)
nullmod <- lmer(tempdf$Shannon_diversity ~ 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
fullmod <- lmer(tempdf$Shannon_diversity ~ allCACvsCRC + 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
test <- anova(nullmod,fullmod)

# nearly significant but model has issues, possibly due to large cis in pts with 2 samples?

boxedIn(tse_cac_vs_crc, "Shannon_diversity", fullmod, "allCACvsCRC", x_lab = "Group", y_lab = "Shannon diversity value")

kruskal.test(chao1 ~ allCACvsCRC, data = colData(tse_cac_vs_crc))
# not at all significant, check lme

tempdf <- tse_cac_vs_crc@colData@listData
tempdf$JLA <- factor(tempdf$JLA)
nullmod <- lmer(tempdf$chao1 ~ 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
fullmod <- lmer(tempdf$chao1 ~ Group + 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
test <- anova(nullmod,fullmod)

boxedIn(tse_cac_vs_crc, "chao1", fullmod, "allCACvsCRC", x_lab = "Group", y_lab = "Chao1 diversity value")

# pval narrowly non-significant

#UC vs CD vs sCRC, - stool

tse_uc_vs_cd_vs_crc <- tse[,-which(tse$SampleType == "stool")]
tse_uc_vs_cd_vs_crc <- estimateDiversity(tse_uc_vs_cd_vs_crc, index = indices, name = names)

shannon_box_group_2 <- ggplot(as.data.frame(colData(tse_uc_vs_cd_vs_crc)),
                          aes(x = Group, 
                              y = Shannon_diversity,
                              fill = Group)) + xlab("Disease Status") +
  ylab("Shannon Diversity") +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  guides(fill = "none") +
  theme(title = element_text(size = 12)) # makes titles smaller

shannon_box_group_2

FSA::dunnTest(tse_uc_vs_cd_vs_crc@colData@listData[["Shannon_diversity"]]~tse_uc_vs_cd_vs_crc@colData@listData[["Group"]])
# not significant, check lme, which will also be not significant

tempdf <- tse_uc_vs_cd_vs_crc@colData@listData
tempdf$JLA <- factor(tempdf$JLA)
nullmod <- lmer(tempdf$Shannon_diversity ~ 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
fullmod <- lmer(tempdf$Shannon_diversity ~ Group + 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
test <- anova(nullmod,fullmod)
# closer to but not significant

boxedIn(tse_uc_vs_cd_vs_crc, "Shannon_diversity", fullmod, "Group", x_lab = "Group", y_lab = "Shannon diversity value")

# to do, chao, inv simpson, etc

#AgeBMI vs Histo, -stool

tse_age_vs_histo <- tse[,-which(tse$Type == "")]
tse_age_vs_histo <- tse_age_vs_histo[,-which(tse_age_vs_histo$SampleType == "stool")]
tse_age_vs_histo <- estimateDiversity(tse_age_vs_histo, index = indices, name = names)

kruskal.test(Shannon_diversity ~ Type, data = colData(tse_age_vs_histo))

tempdf <- tse_age_vs_histo@colData@listData
tempdf$JLA <- factor(tempdf$JLA)
nullmod <- lmer(tempdf$Shannon_diversity ~ 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
fullmod <- lmer(tempdf$Shannon_diversity ~ Type + 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
test <- anova(nullmod,fullmod)

# nearly significant but messed up models

shannon_box_abmi <- ggplot(as.data.frame(colData(tse_age_vs_histo)),
aes(x = Type,
y = Shannon_diversity,
fill = Type)) + xlab("Disease Status") +
ylab("Shannon Diversity") +
geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
guides(fill = "none") +
theme(title = element_text(size = 12)) # makes titles smaller
shannon_box_abmi

tse_age_vs_histo_minus <- tse_age_vs_histo[,-which(tse_age_vs_histo$JLA %in% c(340, 263))]
kruskal.test(Shannon_diversity ~ Type, data = colData(tse_age_vs_histo_minus))

shannon_box_abmiminus <- ggplot(as.data.frame(colData(tse_age_vs_histo_minus)),
                           aes(x = Type,
                               y = Shannon_diversity,
                               fill = Type)) + xlab("Disease Status") +
  ylab("Shannon Diversity") +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  guides(fill = "none") +
  theme(title = element_text(size = 12)) # makes titles smaller
shannon_box_abmiminus

tempdf <- tse_age_vs_histo_minus@colData@listData
tempdf$JLA <- factor(tempdf$JLA)
nullmod <- lmer(tempdf$Shannon_diversity ~ 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
fullmod <- lmer(tempdf$Shannon_diversity ~ Type + 1|JLA, data=tempdf, REML = F, control = lmerControl(optimizer ="Nelder_Mead"))
test <- anova(nullmod,fullmod)

# ns plus singular model



#AgeBMI vs CAC, -stool
#Histo vs CAC, -stool
#T vs N, -stool
#T vs N, -stool, only sCRC
#T vs N, -stool, only CAC



# assess alpha diversity by categorical metadata variables

shannon_box_loc <- ggplot(as.data.frame(colData(tse)),
                          aes(x = Group, 
                              y = Shannon_diversity,
                              fill = Group)) + xlab("Disease Status") +
  ylab("Shannon Diversity") +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  guides(fill = "none") +
  theme(title = element_text(size = 12)) # makes titles smaller

shannon_box_loc

ggsave(filename = "shannonboxdisease.png", plot = last_plot(), dpi = 300)

tempdf <- as.data.frame(colData(tse)[-which(is.na(colData(tse)$allCACvsCRC) == TRUE),])
tempdf <- as.data.frame(colData(tse)[-which(is.na(colData(tse)$CDvsUC) == TRUE),])


shannon_box_sampletype <- shannon_box_loc <- ggplot(as.data.frame(tempdf),
                                                    aes(x = CDvsUC, 
                                                        y = Shannon_diversity,
                                                        fill = CDvsUC)) + xlab("Sample Type") +
  ylab("Shannon Diversity") +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  guides(fill = "none") +
  theme(title = element_text(size = 12)) # makes titles smaller


shannon_box_sampletype

ggsave(filename = "shannonboxsampletype.png", plot = last_plot(), dpi = 300)

shannon_box_pt <- ggplot(as.data.frame(colData(tse)),
                         aes(x = ptid, 
                             y = Shannon_diversity,
                             fill = ptid)) + scale_x_discrete(guide = guide_axis(angle = 90)) +
  geom_boxplot() +
  theme(title = element_text(size = 12)) # makes titles smaller

shannon_box_pt

shannon_box_diag <- ggplot(as.data.frame(colData(tse)), aes(x = Diagnosis, 
                                                            y = Shannon_diversity, 
                                                            fill = Diagnosis)) + 
  geom_boxplot() +
  theme(title = element_text(size = 12)) # makes titles smaller

shannon_box_diag

#faith_box_pt <- ggplot(as.data.frame(colData(tse)), aes(x = ptid, 
#                                                        y = Faith_diversity, 
#                                                        fill = ptid)) + 
#  geom_boxplot() +
#  theme(title = element_text(size = 12)) # makes titles smaller

#faith_box_pt

# test for significant differences with single or multiple group comparisons
#wilcox.test(Shannon_diversity ~ Diagnosis, data = colData(tse))
kruskal.test(Shannon_diversity ~ Disease, data = colData(tse))
kruskal.test(Shannon_diversity ~ Diagnosis, data = colData(tse))

tempdf <- data.frame(colData(tse)[-which(is.na(colData(tse)$allCACvsCRC) == TRUE),])


# FSA fails for two groups
FSA::dunnTest(tempdf$Shannon_diversity~tempdf$Group)
FSA::dunnTest(tempdf$Shannon_diversity~tempdf$Serial)

FSA::dunnTest(tse@colData@listData[["Shannon_diversity"]]~tse@colData@listData[["Group"]])

kruskal.test(Shannon_diversity ~ Sample_type, data = colData(tse))

tse <- mia::estimateRichness(tse, 
                             assay_name = "counts", 
                             index = "chao1", 
                             name="chao1")



head(colData(tse)$chao1)

# assess chao1 richness by categorical metadata variables

chao_box_dis <- ggplot(as.data.frame(colData(tse)),
                       aes(x = Group, 
                           y = chao1,
                           fill = Group)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  theme(title = element_text(size = 12)) # makes titles smaller

chao_box_dis

chao_box_sampletype <- ggplot(as.data.frame(colData(tse)),
                              aes(x = Sample_type, 
                                  y = chao1,
                                  fill = Sample_type)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + xlab("Sample Type") + ylab("Chao1 Index") +
  guides(fill = "none") +
  theme(title = element_text(size = 12)) # makes titles smaller

chao_box_sampletype

kruskal.test(chao1 ~ Disease, data = colData(tse))

FSA::dunnTest(tse@colData@listData[["chao1"]]~tse@colData@listData[["Disease"]])

FSA::dunnTest(tse@colData@listData[["chao1"]]~tse@colData@listData[["Sample_type"]], method = "bh")

# remove Natasha75 and Natasha96

# transform counts to relative abundance, original counts are kept
# tables can be accessed via tse@assays@data$counts/relabundance/clr etc
tse <- transformCounts(tse, abund_values = "counts", method = "relabundance")

# plot relative abundance bar charts

tse <- relAbundanceCounts(tse)

# Getting top taxa on a Phylum level
tse_phylum <- agglomerateByRank(tse, rank ="Phylum", onRankOnly=TRUE)
top_taxa <- getTopTaxa(tse_phylum,top = 5, assay_name = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse)$Phylum,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse)$Phylum <- as.character(phylum_renamed)

tse_genus <- agglomerateByRank(tse, rank ="Genus", onRankOnly=TRUE)
top_genera <- getTopTaxa(tse_genus,top = 20, assay_name = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
genus_renamed <- lapply(rowData(tse)$Genus,
                         function(x){if (x %in% top_genera) {x} else {"Other"}})
rowData(tse)$Genus <- as.character(genus_renamed)

# cribbed from history after a crash, cleanup
# Visualizing the composition barplot, with samples order by "Bacteroidetes"
testplot <- plotAbundance(tse, rank = "Phylum",
                          order_rank_by="abund", order_sample_by = "Fusobacteriota", add_x_text = TRUE)

df2 <- data.frame(x = c(levels(testplot[["plot_env"]][["object"]][["X"]])), colour = factor(tse@colData$Group[match(c(levels(testplot[["plot_env"]][["object"]][["X"]])), tse@colData$Name)]))
df2$pal <- 0
df2$pal[which(df2$colour == "CD")] <- "#00BA38"
df2$pal[which(df2$colour == "UC")] <- "#619CFF"
df2$pal[which(df2$colour == "sCRC")] <- "#F8766D"

df3 <- data.frame(x = c(levels(testplot[["plot_env"]][["object"]][["X"]])), colour = factor(tse@colData$TvsN[match(c(levels(testplot[["plot_env"]][["object"]][["X"]])), tse@colData$Name)]))
df3$pal <- 0
df3$pal[which(df3$colour == "N")] <- "#619CFF"
df3$pal[which(df3$colour == "T")] <- "#F8766D"


testplot2 <- testplot + geom_tile(data=df2, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df2$colour), values = df2$pal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.direction="horizontal") +labs(fill = "Diagnosis")
testplot2
testplot3 <- testplot + geom_tile(data=df3, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df3$colour), values = df3$pal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.direction="horizontal") +labs(fill = "Sample type")

testplot3
legend1 <- as_ggplot(get_legend(testplot))
legend2 <- as_ggplot(get_legend(testplot2))
legend3 <- as_ggplot(get_legend(testplot3))
#ggarrange(testplot, legend1, testplot2, legend2, ncol = 2, nrow = 2, heights = c(10, 1), widths = c(8,2), align = 'h')
testplot2 <- testplot + geom_tile(data=df2, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df2$colour), values = df2$pal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.position="none") +labs(fill = "Diagnosis")
testplot2

testplot3 <- testplot + geom_tile(data=df3, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df3$colour), values = df3$pal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.position="none") +labs(fill = "Sample type")
testplot3

testplot <- plotAbundance(tse, rank = "Phylum",
                          order_rank_by="abund", order_sample_by = "Fusobacteriota") + scale_x_discrete(labels = tse@colData$ptid[match(c(levels(testplot[["plot_env"]][["object"]][["X"]])), tse@colData$X.NAME)]) + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
testplot
ggarrange(testplot, legend1, testplot2, legend2, testplot3, legend3, ncol = 2, nrow = 3, heights = c(10, 1, 1), widths = c(8,2,2), align = 'h')
# ensure rendering window size is large enough to show the geom_tile layer before saving
ggsave(filename = "fusorelabund.png", plot = last_plot(), dpi = 300)

#############################################

# cribbed from history after a crash, cleanup
# Visualizing the composition barplot, with samples order by "Bacteroidetes"
testplot <- plotAbundance(tse, rank = "Genus",
                          order_rank_by="abund", order_sample_by = "Fusobacteria", add_x_text = TRUE)

df2 <- data.frame(x = c(levels(testplot[["plot_env"]][["object"]][["X"]])), colour = factor(tse@colData$Group[match(c(levels(testplot[["plot_env"]][["object"]][["X"]])), tse@colData$Name)]))
df2$pal <- 0
df2$pal[which(df2$colour == "CD")] <- "#00BA38"
df2$pal[which(df2$colour == "UC")] <- "#619CFF"
df2$pal[which(df2$colour == "sCRC")] <- "#F8766D"

df3 <- data.frame(x = c(levels(testplot[["plot_env"]][["object"]][["X"]])), colour = factor(tse@colData$TvsN[match(c(levels(testplot[["plot_env"]][["object"]][["X"]])), tse@colData$Name)]))
df3$pal <- 0
df3$pal[which(df3$colour == "N")] <- "#619CFF"
df3$pal[which(df3$colour == "T")] <- "#F8766D"


testplot2 <- testplot + geom_tile(data=df2, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df2$colour), values = df2$pal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.direction="horizontal") +labs(fill = "Diagnosis")
testplot2
testplot3 <- testplot + geom_tile(data=df3, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df3$colour), values = df3$pal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.direction="horizontal") +labs(fill = "Sample type")

testplot3
legend1 <- as_ggplot(get_legend(testplot))
legend2 <- as_ggplot(get_legend(testplot2))
legend3 <- as_ggplot(get_legend(testplot3))
#ggarrange(testplot, legend1, testplot2, legend2, ncol = 2, nrow = 2, heights = c(10, 1), widths = c(8,2), align = 'h')
testplot2 <- testplot + geom_tile(data=df2, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df2$colour), values = df2$pal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.position="none") +labs(fill = "Diagnosis")
testplot2

testplot3 <- testplot + geom_tile(data=df3, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df3$colour), values = df3$pal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.position="none") +labs(fill = "Sample type")
testplot3

testplot <- plotAbundance(tse, rank = "Phylum",
                          order_rank_by="abund", order_sample_by = "Fusobacteriota") + scale_x_discrete(labels = tse@colData$ptid[match(c(levels(testplot[["plot_env"]][["object"]][["X"]])), tse@colData$X.NAME)]) + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
testplot
ggarrange(testplot, legend1, testplot2, legend2, testplot3, legend3, ncol = 2, nrow = 3, heights = c(10, 1, 1), widths = c(8,2,2), align = 'h')
# ensure rendering window size is large enough to show the geom_tile layer before saving
ggsave(filename = "fusorelabund.png", plot = last_plot(), dpi = 300)


tse_naf <- tse[ , tse$Sample_type %in% "Aspirate_FluidSwab"]
tse_naf <- transformCounts(tse_naf, abund_values = "counts", method = "relabundance")
tse_naf_vs_skin <- tse[ , tse$Sample_type %in% c("Aspirate_FluidSwab", "Nipple_SkinSwab")]
tse_naf_vs_skin <- transformCounts(tse_naf_vs_skin, abund_values = "counts", method = "relabundance")


taxtable_abs <- read_csv("taxa_silva_Microbiomeanalyst.csv")
taxtable_abs <- data.frame(taxtable_abs)
seqtab_abs <- read_csv("abs_abundance_fixed.csv")
seqtab_abs <- data.frame(seqtab_abs)
metadata_abs <- read_csv("Meta_basic_MicrobiomeAnalyst_FULL.csv")
metadata_abs <- data.frame(metadata_abs)

test <- t(as.data.frame(strsplit(metadata_abs$SampleID, "_")))
metadata_abs$ptid <- test[,2]

metadata_abs[which(metadata_abs$Diagnosis == "normal"),7] <- "Normal"



# fix imported non standard characters
taxtable_abs <- column_to_rownames(taxtable_abs, "X.TAXONOMY")
seqtab_abs <- column_to_rownames(seqtab_abs, "X.NAME")
seqtab_abs <- t(seqtab_abs)

# remove samples with failed sequencing

`%notin%` <- Negate(`%in%`)
which(metadata_abs$X.NAME %notin% colnames(seqtab_abs))

metadata_abs <- metadata_abs[-which(metadata_abs$X.NAME %notin% colnames(seqtab_abs)),]

# arrange seqtab columns in same order as metadata samples so it doesn't do wacky stuff
seqtab_abs <- seqtab_abs[,c(metadata_abs$X.NAME)]


# create a copy of the taxonomy table with NA values replaced with best available taxonomy
taxtable_nona_abs <- as.matrix(taxtable_abs)
taxtable_nona_abs <- turdpolisher(taxtable_nona_abs, seqtab_abs)
# rerun with outputseqtab
taxtable_nona_abs <- as.matrix(taxtable_abs)
taxtable_nona_abs <- turdpolisher(as.matrix(outputtaxtable), as.matrix(outputseqtab))

#metadata_abs <- metadata_abs[which(metadata_abs$X.NAME %in% colnames(seqtab_abs)),]
rownames(metadata_abs) <- NULL
metadata_abs <- column_to_rownames(metadata_abs, var = "X.NAME")
metadata_abs <- metadata_abs[colnames(seqtab_abs),]

# store sequences for later use
seqtab_abs <- outputseqtab[,colnames(outputseqtab) %in% rownames(metadata_abs)]

sequences_abs <- as.data.frame(rownames(seqtab_abs))

# give sequences more manageable names
for(i in 1:nrow(seqtab_abs)){
  rownames(seqtab_abs)[i] <- paste0("ASV_", i)
}

# store corresponding ASV number for ease of use in subsetting
rownames(taxtable_nona) <- rownames(seqtab_abs)
sequences_abs[,2] <- rownames(seqtab_abs)

#seqtab <- read.csv("seqtab_fixed.csv")
#seqtab <- column_to_rownames(seqtab, "X.NAME")

tse_abs <- TreeSummarizedExperiment(assays = list(counts = seqtab_abs),
                                    colData = metadata_abs,
                                    rowData = taxtable_nona)
#rowTree = phylo)

#which(colnames(tse_abs) == "Natasha75")
#tse_abs <- tse_abs[,-73]
#
#which(colnames(tse) == "Natasha96")
#tse <- tse[,-95]

tse_abs <- tse_abs[,-which(tse_abs@colData@listData[["ptid"]] == "TB80")]

# Getting top taxa on a Phylum level
tse_abs_phylum <- agglomerateByRank(tse_abs, rank ="Phylum", onRankOnly=TRUE)
top_taxa <- getTopTaxa(tse_abs_phylum,top = 5, assay_name = "counts")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse_abs)$Phylum,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse_abs)$Phylum <- as.character(phylum_renamed)

# plot absolute abundance barplot for comparison

test <- plotAbundance(tse_abs, rank = "Phylum",
                      order_rank_by="abund", add_x_text = TRUE, use_relative = FALSE)

plotAbundance(tse_abs, rank = "Phylum",
              order_rank_by="abund", use_relative = FALSE) + scale_x_discrete(labels = tse_abs@colData$ptid[match(c(levels(test[["plot_env"]][["object"]][["X"]])), rownames(tse_abs@colData))]) + theme(axis.text.x = element_text(angle = 90), legend.position = "none")

# censor natasha69 and reproduce

plotAbundance(tse_abs[,-1], rank = "Phylum",
              order_rank_by="abund", use_relative = FALSE) + scale_x_discrete(labels = tse_abs@colData$ptid[match(c(levels(test[["plot_env"]][["object"]][["X"]])), rownames(tse_abs@colData))]) + theme(axis.text.x = element_text(angle = 90), legend.position = "none")


# add total sum to coldata for plotting absolute abs in order

tse_abs@colData$totalsum <- colSums(tse_abs@assays@data@listData[["counts"]])

test <- plotAbundance(tse_abs, rank = "Phylum",
                      order_rank_by="abund", order_sample_by = "totalsum", add_x_text = TRUE, use_relative = FALSE)


# add tile nonsense below

df2 <- data.frame(x = c(levels(test[["plot_env"]][["object"]][["X"]])), colour = factor(tse_abs@colData$Disease[match(c(levels(test[["plot_env"]][["object"]][["X"]])), colnames(tse_abs))]))
df2$pal <- 0
df2$pal[which(df2$colour == "Benign")] <- "#00BA38"
df2$pal[which(df2$colour == "Normal")] <- "#619CFF"
df2$pal[which(df2$colour == "Cancer")] <- "#F8766D"


bgplot <- plotAbundance(tse_abs, rank = "Phylum",
                        order_rank_by="abund", order_sample_by = "totalsum", add_x_text = TRUE) + theme(axis.title.x=element_blank(),
                                                                                                        axis.text.x=element_blank(),
                                                                                                        axis.ticks.x=element_blank(),
                                                                                                        axis.text.y = element_text(colour = "white"),
                                                                                                        axis.ticks.y=element_line(colour = "white"),
                                                                                                        axis.title.y=element_text(colour = "white"))

testplot2 <- bgplot + geom_tile(data=df2, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df2$colour), values = df2$pal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.direction="horizontal") +labs(fill = "Diagnosis")
testplot2

testplot <- plotAbundance(tse_abs, rank = "Phylum",
                          order_rank_by="abund",order_sample_by = "totalsum", add_x_text = TRUE, use_relative = FALSE)

legend1 <- as_ggplot(get_legend(testplot))
legend2 <- as_ggplot(get_legend(testplot2))
#ggarrange(testplot, legend1, testplot2, legend2, ncol = 2, nrow = 2, heights = c(10, 1), widths = c(8,2), align = 'h')
testplot2 <- bgplot + geom_tile(data=df2, aes(x = x, y = 0.5, fill = colour)) + scale_fill_manual(breaks=c(df2$colour), values = df2$pal) +
  theme(axis.title.x=element_blank(),
        line = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "white"),
        axis.ticks.y=element_line(colour = "white"),
        axis.title.y=element_text(colour = "white"), legend.position="none") +labs(fill = "Diagnosis")
testplot2
testplot <- plotAbundance(tse_abs, rank = "Phylum",
                          order_rank_by="abund",order_sample_by = "totalsum", use_relative = FALSE) + scale_x_discrete(labels = tse_abs@colData$ptid[match(c(levels(test[["plot_env"]][["object"]][["X"]])), rownames(tse_abs@colData))]) + theme(axis.text.x = element_text(angle = 90), legend.position = "none")
testplot <- testplot + ylab("Adjusted Counts")
ggarrange(testplot, legend1, testplot2, legend2, ncol = 2, nrow = 2, heights = c(10, 1), widths = c(11,3), align = 'hv')
# ensure rendering window size is large enough to show the geom_tile layer before saving
ggsave(filename = "absbarplotwdiags.png", plot = last_plot(), dpi = 300)

# prevalence filter

tse_genus <- agglomerateByRank(tse, "Genus")
tse_genus <- subsetByPrevalentTaxa(tse_genus, detection = 0, prevalence = 0.15)

tse_prev15 <- subsetByPrevalentTaxa(tse, detection = 0, prevalence = 0.15)

tse_prev15 <- transformCounts(tse, method = "clr", pseudocount = 1)
tse_genus <- transformCounts(tse_genus, method = "clr", pseudocount = 1)

doPCA <- function(input_tse, Group) {
  clr_pca <- prcomp(t(input_tse@assays@data@listData[["rclr"]]))
  df_out <- as.data.frame(clr_pca$x)
  
  # construct df of all pc combinations, warning, if you have a lot this might die
  #allpcs <- as.data.frame(t(combn(colnames(df_out), 2)))
  
  df_out$Group <- as.factor(input_tse@colData@listData$Group)
  
  pcx <- "PC1"
  pcy <- "PC2"
  
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
}

permanova_disease2 <- vegan::adonis(t(tse@assays@data$relabundance) ~ Disease,
                                    data = colData(tse),
                                    permutations = 9999)

permanova_disease2[["aov.tab"]]

# autopca

# perform pca of tse object to be used
clr_pca <- prcomp(t(tse_genus@assays@data@listData[["clr"]]))
pcasummary <- data.frame(summary(clr_pca)[["importance"]])
df_out <- as.data.frame(clr_pca$x)


# construct df of all pc combinations, warning, if you have a lot this might die
allpcs <- as.data.frame(t(combn(colnames(df_out), 2)))

# append metadata column as factor
df_out$disease <- as.factor(tse_genus@colData$Sample_type)

# remove benign for ease of viz

df_out <- df_out[-which(df_out$disease == "Benign"),]

pcx <- "PC1"
pcy <- "PC2"
pcresults <- allpcs
pcresults$areaoverlap <- 0
pcresults$areaoverlappct <- 0
pcresults$areanooverlap <- 0
pcresults$areanooverlappct <- 0
pcresults$variancenooverlap <- 0
pcresults$pctpointsoutside <- 0

condition1 <- "Cancer"
condition2 <- "Normal"

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=disease))
p<-p+geom_point()
p

for(i in 1:nrow(allpcs)){
  pcx <- allpcs[i,1]
  pcy <- allpcs[i,2]
  
  hull <- df_out %>%
    group_by(metadata$Disease[-which(metadata$Disease == "Benign")]) %>%
    dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))
  
  tempdf1 <- df_out[pcx]
  tempdf2 <- df_out[pcy]
  
  points1 <- point.in.polygon(data.frame(df_out[pcx])[which(df_out$disease == condition1),], data.frame(df_out[pcy])[which(df_out$disease == condition1),], data.frame(hull[pcx])[which(hull$disease == condition2),], data.frame(hull[pcy])[which(hull$disease == condition2),])
  points1df <- df_out[which(df_out$disease == condition1),]
  points1df <- points1df[which(points1 != 0),]
  
  points2 <- point.in.polygon(data.frame(df_out[pcx])[which(df_out$disease == condition2),], data.frame(df_out[pcy])[which(df_out$disease == condition2),], data.frame(hull[pcx])[which(hull$disease == condition1),], data.frame(hull[pcy])[which(hull$disease == condition1),])
  points2df <- df_out[which(df_out$disease == condition2),]
  points2df <- points2df[which(points2 != 0),]
  
  # poly_crossings needs filled polygons, duplicate first row as last row
  poly1 <- as.matrix(cbind(data.frame(hull[pcx])[which(hull$disease == condition1),], data.frame(hull[pcy])[which(hull$disease == condition1),]))
  poly2 <- as.matrix(cbind(data.frame(hull[pcx])[which(hull$disease == condition2),], data.frame(hull[pcy])[which(hull$disease == condition2),]))
  poly1 <- as.data.frame(poly1)
  poly2 <- as.data.frame(poly2)
  poly1[(nrow(poly1) +1),] <- poly1[1,]
  poly2[(nrow(poly2) +1),] <- poly2[1,]
  
  
  
  crossings <- poly_crossings(t(poly1), t(poly2))
  crossings <- as.data.frame(crossings)
  if(length(crossings) != 0){  
    colnames(crossings) <- c(pcx, pcy)
    crossings$disease <- "Overlap"
    
    overlappoints <- rbind(points1df[,which(colnames(points1df) %in% colnames(crossings))], points2df[,which(colnames(points2df) %in% colnames(crossings))], crossings)
    overlappoints$disease <- "Overlap"
    
    pcresults$pctpointsoutside[i] <- ((nrow(tempdf1) - nrow(overlappoints))/nrow(tempdf1))*100
    
    overlaphull <- overlappoints %>%
      dplyr::slice(chull(overlappoints))
    
    
    
    vec1 <- as.vector(overlaphull[,1])
    vec2 <- as.vector(overlaphull[,2])
    
    (diff(range(tempdf1)) * diff(range(tempdf2)))
    
    pcresults$areaoverlap[i] <- abs(polyarea(vec1, vec2))
    pcresults$areaoverlappct[i] <- (abs(polyarea(vec1, vec2)) / (diff(range(tempdf1)) * diff(range(tempdf2)))) *100
    
    poly1areaminus <- abs(polyarea(poly1[,1], poly1[,2])) - abs(polyarea(overlaphull[,1], overlaphull[,2]))
    poly2areaminus <- abs(polyarea(poly2[,1], poly2[,2])) - abs(polyarea(overlaphull[,1], overlaphull[,2]))
    pcresults$areanooverlap[i] <- (poly1areaminus + poly2areaminus)
    
    pcresults$areanooverlappct[i] <- (poly1areaminus + poly2areaminus) / (diff(range(tempdf1)) * diff(range(tempdf2))) * 100
    pcresults$variancenooverlap[i] <- ((pcasummary[2,pcx] + pcasummary[2,pcy])/100)*pcresults$areanooverlappct[i]
  }
  
  else
    
  {
    if(isTRUE(0 %in% points1 == FALSE)){
      pcresults$areaoverlap[i] <- paste(condition1, "100% within", condition2)
    }
    
    if(isTRUE(0 %in% points2 == FALSE)){
      pcresults$areaoverlap[i] <- paste(condition2, "100% within", condition1)
    }
    
  }
  
}

pcresults2 <- as.data.frame(matrix(ncol = 4, data = as.numeric(unlist(pcresults))))
pcresults2[,1:2] <- pcresults[,1:2]
pcresults2 <- pcresults2[-which(is.na(pcresults2[,3]) == TRUE),]

# get convex hulls

p<-ggplot(df_out,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=disease))
p<-p+geom_point()
p

hull <- df_out %>%
  group_by(metadata$Disease[-which(metadata$Disease == "Benign")]) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p + aes(fill = factor(disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") + geom_polygon(data = overlaphull)

condition1 <- "Cancer"
condition2 <- "Normal"

points1 <- point.in.polygon(df_out$PC1[which(df_out$disease == condition1)], df_out$PC2[which(df_out$disease == condition1)], hull$PC1[which(hull$disease == condition2)], hull$PC2[which(hull$disease == condition2)])
points1df <- df_out[which(df_out$disease == condition1),]
points1df <- points1df[which(points1 != 0),]

points2 <- point.in.polygon(df_out$PC1[which(df_out$disease == condition2)], df_out$PC2[which(df_out$disease == condition2)], hull$PC1[which(hull$disease == condition1)], hull$PC2[which(hull$disease == condition1)])
points2df <- df_out[which(df_out$disease == condition2),]
points2df <- points2df[which(points2 != 0),]

# poly_crossings needs filled polygons, duplicate first row as last row
poly1 <- poly1 <- as.matrix(cbind(hull$PC1[which(hull$disease == condition1)], hull$PC2[which(hull$disease == condition1)]))
poly2 <- as.matrix(cbind(hull$PC1[which(hull$disease == condition2)], hull$PC2[which(hull$disease == condition2)]))
poly1 <- as.data.frame(poly1)
poly2 <- as.data.frame(poly2)
poly1[(nrow(poly1) +1),] <- poly1[1,]
poly2[(nrow(poly2) +1),] <- poly2[1,]

crossings <- poly_crossings(t(poly1), t(poly2))
crossings <- as.data.frame(crossings)

colnames(crossings) <- c("PC1", "PC2")
crossings$disease <- "Overlap"

overlappoints <- rbind(points1df[,which(colnames(points1df) %in% colnames(crossings))], points2df[,which(colnames(points2df) %in% colnames(crossings))], crossings)
overlappoints$disease <- "Overlap"

overlaphull <- overlappoints %>%
  dplyr::slice(chull(overlappoints))

p + aes(fill = factor(disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") + geom_polygon(data = overlaphull)

perfect_separation <- 
  
  
  
  
  
  
  
  tse_naf_abs <- tse_abs[ , tse_abs$Sample_type %in% "Aspirate_FluidSwab"]

tse_naf_abs_prev10 <- subsetByPrevalentTaxa(tse_naf_abs, detection = 0, prevalence = 0.1)
tse_naf_abs_prev10 <- tse_naf_abs_prev10[ , tse_naf_abs_prev10$Disease %in% c("Cancer", "Normal")]
tse_naf_abs_prev30 <- subsetByPrevalentTaxa(tse_naf_abs, detection = 0, prevalence = 0.3)
tse_naf_abs_prev30 <- tse_naf_abs_prev30[ , tse_naf_abs_prev30$Disease %in% c("Cancer", "Normal")]
tse_naf_abs_prev50 <- subsetByPrevalentTaxa(tse_naf_abs, detection = 0, prevalence = 0.5)
tse_naf_abs_prev50 <- tse_naf_abs_prev50[ , tse_naf_abs_prev50$Disease %in% c("Cancer", "Normal")]

tse_abs_prev10 <- tse_abs[ , tse_abs$Disease != "Benign"]
tse_abs_prev10 <- subsetByPrevalentTaxa(tse_abs_prev10, detection = 0, prevalence = 0.1)
tse_abs_prev15 <- tse_abs[ , tse_abs$Disease != "Benign"]
tse_abs_prev15 <- subsetByPrevalentTaxa(tse_abs_prev15, detection = 0, prevalence = 0.15)
tse_abs_prev30 <- tse_abs[ , tse_abs$Disease != "Benign"]
tse_abs_prev30 <- subsetByPrevalentTaxa(tse_abs_prev30, detection = 0, prevalence = 0.3)
tse_abs_prev50 <- tse_abs[ , tse_abs$Disease != "Benign"]
tse_abs_prev50 <- subsetByPrevalentTaxa(tse_abs_prev50, detection = 0, prevalence = 0.5)


# you can store individual tables as objects if needed via 2 methods
rel_abund_assay <- assays(tse_abs)$relabundance
rel_abund_assay2 <- tse@assays@data$relabundance
# you can see both objects are identical
which((rel_abund_assay == rel_abund_assay2) == "FALSE")

# create a bray curtis distance matrix, this is used to plot PCoA
# but distance matrices do not retain the information about the features that
# contribute to the distances, so traditional loadings plots cannot be generated

bray_curtis_dist <- vegan::vegdist(t(rel_abund_assay), method = "bray")
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)

# you can access different prinicpal coordinate vectors via bray_curtis_pcoa$vectors
# however, PCoA maximises DISTANCE not variance, the first two therefore should
# give you maximal separation, but not necessarily for your variable of interest

bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])

# Create a plot, change colour value to analyse different groups
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point(aes(color = tse@colData@listData$Disease)) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA",
       color = "Disease") +
  theme(title = element_text(size = 10)) # makes titles smaller

bray_curtis_plot

# gethulls code is slightly messy, col1/2 are the column names of your columns of interest

bray_curtis_pcoa_df$disease <- tse@colData@listData$Disease
bray_curtis_pcoa_df$disease <- as.factor(bray_curtis_pcoa_df$disease)

pcx <- "pcoa1"
pcy <- "pcoa2"

hull <- bray_curtis_pcoa_df %>%
  group_by(disease) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

#hull <- gethulls(bray_curtis_pcoa_df, tse@colData@listData$Disease, col1 = "pcoa1", col2 = "pcoa2")

bray_curtis_plot + geom_polygon(data = hull, aes(fill = factor(disease)), alpha = 0.1) + guides(fill = "none")

# clr transformed data has statistical advantages and is the preferred way of
# assessing compositional data these days
# WARNING - clr uses the geometric mean, this will be affected by prevalence filtering
# consider re-running transformations on your count table after filtering
tse <- transformCounts(tse, method = "clr", pseudocount = 1)

# Gets clr table
clr_assay <- assays(tse_prev15)$clr

# Transposes it to get taxa to columns
clr_assay <- t(clr_assay)

# Calculates Euclidean distances between samples. Because taxa is in columns,
# it is used to compare different samples.
euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")

# Does principal coordinate analysis
euclidean_pcoa <- ecodist::pco(euclidean_dist)

# Creates a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1], 
                                pcoa2 = euclidean_pcoa$vectors[,2])

# Creates the plot
euclidean_plot <- ggplot(data = euclidean_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point(aes(color = factor(colData(tse)$Disease))) +
  labs(x = "PC1",
       y = "PC2",
       title = "Euclidean PCoA with CLR transformation") +
  labs(color = "Disease")
theme(title = element_text(size = 12)) # makes titles smaller

euclidean_plot

euclidean_pcoa_df$disease <- tse_prev15@colData@listData$Disease
euclidean_pcoa_df$disease <- as.factor(euclidean_pcoa_df$disease)

pcx <- "pcoa1"
pcy <- "pcoa2"

hull <- euclidean_pcoa_df %>%
  group_by(disease) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

euclidean_plot + geom_polygon(data = hull, aes(fill = factor(disease)), alpha = 0.1) + guides(fill = "none", color = guide_legend(title = "T stage"))



# PERMANOVA to test association with conditions/taxa
# do not use CLR transformed data for this

# ensure factors are factors
#tse@colData@listData$Disease <- as.factor(tse@colData@listData$Disease)

set.seed(123)
permanova_disease2 <- vegan::adonis(t(tse_genus@assays@data$relabundance) ~ Disease,
                                    data = colData(tse),
                                    permutations = 99999)

permanova_disease2[["aov.tab"]]


# you can see that though there are significant differences between T stages
# the R2 value is relatively low, this is common with microbiome data but what
# other factors may explain variance in the data?

# fix ptid to factors

#tse@colData@listData[["SampleID"]] <- as.factor(tse@colData@listData[["SampleID"]])

permanova_sample <- vegan::adonis(t(tse@assays@data$relabundance) ~ ptid,
                                  data = colData(tse),
                                  permutations = 9999)

View(permanova_sample[["aov.tab"]])

# from this model we see that while lots of variables have significantly different
# microbiota between them, that a huge amount of the difference between them is
# simply the between patient differences.

# be warned that adonis will produce different results based on the order of the model call
# use adonis2 for consistency, HOWEVER you will not be able to extract coefficients

# extract the model coefficients from the first model, filter to top 10

coef <- coefficients(permanova_disease)["Disease1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

# make the names look pretty
names[which(startsWith(rowData(tse_genus_abs_nostool)[names(top.coef), ][,"Genus"], "Unclassified_") == FALSE)] <- paste(rowData(tse_genus_abs_nostool)[names(top.coef), ][,"Family"][which(startsWith(rowData(tse_genus_abs_nostool)[names(top.coef), ][,"Genus"], "Unclassified_") == FALSE)], rowData(tse_genus_abs_nostool)[names(top.coef), ][,"Genus"][which(startsWith(rowData(tse_genus_abs_nostool)[names(top.coef), ][,"Genus"], "Unclassified_") == FALSE)])

# plot coefficients with family+genus name
top_taxa_coeffient_plot <- ggplot(data.frame(x = top.coef,
                                             y = factor(names(top.coef),
                                                        unique(names(top.coef)))),
                                  aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="", y="", title="") +
  scale_y_discrete(labels = names) +
  theme_bw()

top_taxa_coeffient_plot

# look at pairwise differences

test <- pairwise.adonis(x=t(tse_genus_abs_nostool@assays@data$counts),factors=tse_genus_abs_nostool@colData@listData$allCACvsCRC,sim.function='vegdist',
                        sim.method='euclidean',p.adjust.m='BH')

# nd

test <- pairwise.adonis(x=t(tse@assays@data$relabundance),factors=tse@colData@listData$Diagnosis,sim.function='vegdist',
                        sim.method='euclidean',p.adjust.m='BH')

# no fdr corrected differences, only non-corrected differences are between a group with one patient

test <- pairwise.adonis(x=t(tse@assays@data$relabundance),factors=tse@colData@listData$Smoker,sim.function='vegdist',
                        sim.method='euclidean',p.adjust.m='BH')

# nd

test <- pairwise.adonis(x=t(tse@assays@data$relabundance),factors=tse@colData@listData$Grade,sim.function='vegdist',
                        sim.method='euclidean',p.adjust.m='BH')

# need to discuss

test <- pairwise.adonis(x=t(tse@assays@data$relabundance),factors=tse@colData@listData$Ethnicity,sim.function='vegdist',
                        sim.method='euclidean',p.adjust.m='BH')

# nd

test <- pairwise.adonis(x=t(tse@assays@data$relabundance),factors=tse@colData@listData$ptid,sim.function='vegdist',
                        sim.method='euclidean',p.adjust.m='BH')

# again explains the most but dies on fdr


tse_genus <- agglomerateByRank(tse, "Genus")
tse_genus <- subsetByPrevalentTaxa(tse_genus, detection = 0, prevalence = 0.15)

associations <- multipatt(t(tse@assays@data$clr), droplevels(as.factor(tse@colData@listData$Disease)), func = "r.g", control = how(nperm=9999))
sigassoc <- associations[["sign"]][which(associations[["sign"]][,"p.value"] < 0.05),]

associations_genus <- multipatt(t(tse_genus@assays@data$clr), droplevels(as.factor(tse_genus@colData@listData$Disease)), func = "r.g", control = how(nperm=9999))
sigassoc_genus <- associations_genus[["sign"]][which(associations_genus[["sign"]][,"p.value"] < 0.05),]

# make pretty names as before
prettynames <- rowData(tse)[,"Genus"]
#prettynames[which(startsWith(rowData(tse)[,"Genus"], "Unclassified_") == TRUE)] <- paste(rowData(tse)[which(startsWith(rowData(tse)[,"Genus"], "Unclassified_") == TRUE),"Family"], rowData(tse)[which(startsWith(rowData(tse)[,"Family"], "Unclassified_") == TRUE),"Genus"])
prettynames <- gsub("_", " ", prettynames)

index <- gsub("ASV_", "", rownames(sigassoc))
for(i in 1:length(index)){
  index[i] <- prettynames[as.numeric(index[i])]
}

coefs <- associations[["str"]][which(rownames(associations[["str"]]) %in% rownames(sigassoc) == TRUE),]
coefs1 <- coefs[,1]
coefs2 <- coefs[,2]
coefs3 <- coefs[,3]
coefs1 <- coefs1[order(coefs[,1])]
coefs2 <- coefs2[order(coefs[,2])]
coefs3 <- coefs3[order(coefs[,3])]

names <- NULL
for(i in 1:length(coefs1)){
  names[i] <- tse@rowRanges@elementMetadata@listData[["Genus"]][which(names(coefs1)[i] == rownames(tse@assays@data$relabundance))]
}

ggplot(data.frame(x = coefs1,
                  y = factor(names(coefs1), unique(names(coefs1)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="Point biserial correlation coefficient", y="", title="Significant taxa associations with benign") +
  scale_y_discrete(labels = names) +
  theme_bw()

names <- NULL
for(i in 1:length(coefs2)){
  names[i] <- tse@rowRanges@elementMetadata@listData[["Genus"]][which(names(coefs2)[i] == rownames(tse@assays@data$relabundance))]
}

ggplot(data.frame(x = coefs2,
                  y = factor(names(coefs2), unique(names(coefs2)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="Point biserial correlation coefficient", y="", title="Significant taxa associations with cancer") +
  scale_y_discrete(labels = names) +
  theme_bw()

names <- NULL
for(i in 1:length(coefs3)){
  names[i] <- tse@rowRanges@elementMetadata@listData[["Genus"]][which(names(coefs3)[i] == rownames(tse@assays@data$relabundance))]
}

ggplot(data.frame(x = coefs3,
                  y = factor(names(coefs3), unique(names(coefs3)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="Point biserial correlation coefficient", y="", title="Significant taxa associations with normal") +
  scale_y_discrete(labels = names) +
  theme_bw()

coef_pca <- prcomp(coefs[,1:2])
df_out <- as.data.frame(coef_pca$x)
df_out$disease <- sigassoc$index

hull <- df_out %>%
  group_by(disease) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=as.factor(disease)))
p<-p+geom_point() + geom_polygon(data = hull, aes(fill = factor(disease)), alpha = 0.1) + guides(fill = "none")

p

# explore whether associations are plausible with bar plots, check for other asvs with same taxonomy

which(tse@rowRanges@elementMetadata@listData[["Genus"]] == "Gemella")

taxa <- rep("ASV_210 (Gemella)", length(tse@assays@data@listData[["clr"]][210,]))
disease <- tse@colData@listData[["Disease"]]
value <- tse@assays@data@listData[["relabundance"]][210,]
data <- data.frame(taxa,disease,value)

data_summary <- data %>% # the names of the new data frame and the data frame to be summarised
  group_by(disease) %>%   # the grouping variable
  dplyr::summarise(mean = mean(value),  # calculates the mean of each group
                   sd = sd(value), # calculates the standard deviation of each group
                   n = n(),  # calculates the sample size per group
                   SE = sd(value)/sqrt(n())) # calculates the standard error of each group

ggplot(data, aes(fill=disease, y=value, x=disease)) + 
  geom_col(position="dodge") #+ geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2)

which(tse_genus@rowRanges@elementMetadata@listData[["Genus"]] == "Aggregatibacter")

taxa <- rep("Aggregatibacter", length(tse_genus@assays@data@listData[["clr"]][62,]))
disease <- tse@colData@listData[["Disease"]]
value <- tse@assays@data@listData[["relabundance"]][62,]
data <- data.frame(taxa,disease,value)

data_summary <- data %>% # the names of the new data frame and the data frame to be summarised
  group_by(disease) %>%   # the grouping variable
  dplyr::summarise(mean = mean(value),  # calculates the mean of each group
                   median = median(value),
                   sd = sd(value), # calculates the standard deviation of each group
                   n = n(),  # calculates the sample size per group
                   SE = sd(value)/sqrt(n())) # calculates the standard error of each group

ggplot(data, aes(fill=disease, y=value, x=disease)) + 
  geom_bar(position="dodge", stat="identity")# + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2)

# to do, paste asv number, assign to some kind of barplot + grouping

associations <- multipatt(t(tse_genus@assays@data$relabundance[,-c(4,22)]), droplevels(as.factor(tse_genus@colData@listData$Disease[-c(4,22)])), func = "r.g", control = how(nperm=9999))
sigassoc <- associations[["sign"]][which(associations[["sign"]][,"p.value"] < 0.15),]
coefs <- associations[["str"]][which(rownames(associations[["str"]]) %in% rownames(sigassoc) == TRUE),]
coefs1 <- coefs[,1]
coefs2 <- coefs[,2]
coefs3 <- coefs[,3]
coefs1 <- coefs1[order(coefs[,1])]
coefs2 <- coefs2[order(coefs[,2])]
coefs3 <- coefs3[order(coefs[,3])]

names <- NULL
for(i in 1:length(coefs1)){
  names[i] <- tse_genus@rowRanges@elementMetadata@listData[["Genus"]][which(names(coefs1)[i] == rownames(tse_genus@assays@data$relabundance))]
}

ggplot(data.frame(x = coefs1,
                  y = factor(names(coefs1), unique(names(coefs1)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="Point biserial correlation coefficient", y="", title="p < 0.15 taxa associations with cancer") +
  scale_y_discrete(labels = names) +
  theme_bw()

names <- NULL
for(i in 1:length(coefs2)){
  names[i] <- tse_genus@rowRanges@elementMetadata@listData[["Genus"]][which(names(coefs2)[i] == rownames(tse_genus@assays@data$relabundance))]
}

ggplot(data.frame(x = coefs2,
                  y = factor(names(coefs2), unique(names(coefs2)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="Point biserial correlation coefficient", y="", title="p < 0.15 taxa associations with normal") +
  scale_y_discrete(labels = names) +
  theme_bw()

names <- NULL
for(i in 1:length(coefs3)){
  names[i] <- tse_genus@rowRanges@elementMetadata@listData[["Genus"]][which(names(coefs3)[i] == rownames(tse_genus@assays@data$relabundance))]
}

ggplot(data.frame(x = coefs3,
                  y = factor(names(coefs3), unique(names(coefs3)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="Point biserial correlation coefficient", y="", title="Significant taxa associations with normal") +
  scale_y_discrete(labels = names) +
  theme_bw()

# change as applicable
ggplot(as.data.frame(tse_genus@assays@data$clr[which(tse_genus@rowRanges@elementMetadata@listData$Genus == "Achromobacter"),]), aes(x=as.numeric(as.factor(tse_genus@colData@listData$Disease)), y=tse_genus@assays@data$clr[which(tse_genus@rowRanges@elementMetadata@listData$Genus == "Achromobacter"),])) +
  xlab("Disease") +
  ylab("CLR transformed abundance") +
  labs(title = "Achromobacter vs disease") +
  scale_x_discrete(limits = c(1, 2, 3), labels = c("Benign", "Cancer", "Normal")) +
  geom_point()+
  geom_smooth(method=lm)

coef_pca <- prcomp(coefs[,1:3])
df_out <- as.data.frame(coef_pca$x)
df_out$group <- groupcol

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=as.factor(group)))
p<-p+geom_point()
p



set.seed(123)
testuni <- runMDS(tse, FUN = calculateUniFrac, name = "UniFrac",
                  tree = rowTree(tse),
                  ntop = nrow(tse),
                  exprs_values = "counts")

# test to see if the scores from the unifrac calculated here are same as phyloseq
# IGNORE THIS WHOLE SECTION UNTIL I REMEMBER WHAT I DID

set.seed(123)
OTU.UF = phyloseq::otu_table(as.matrix(tse@assays@data@listData[["counts"]]), taxa_are_rows=TRUE)
tax.UF = phyloseq::tax_table(as.matrix(taxtable_nona))
meta.UF = phyloseq::sample_data(metadata)

phyloseq::sample_names(meta.UF) <- metadata[,1]
physeq = phyloseq::phyloseq(OTU.UF, tax.UF, meta.UF)
physeq.tree = phyloseq::merge_phyloseq(physeq, tse@rowTree[["phylo"]])
physeq.tree@sam_data[["tumourvsnon"]] <- as.factor(physeq.tree@sam_data[["tumourvsnon"]])


uwUF.ordu = phyloseq::ordinate(physeq.tree, method="NMDS", distance="unifrac", weighted=FALSE)
par(mfrow=c(1,1))
plot(uwUF.ordu, type="n", main="Unweighted UniFrac")
points(uwUF.ordu, pch=20, display="sites", col=c("blue", "red")[physeq.tree@sam_data[["tumourvsnon"]]])


fit <- vegan::envfit(uwUF.ordu, meta.UF)
fit.tumourstatus <- vegan::envfit(uwUF.ordu, meta.UF[,"tumourvsnon"])
fitasv <- vegan::envfit(uwUF.ordu, t(OTU.UF))

fit.scrs <- as.data.frame(scores(fit, display = "vectors"))
asv.scrs <- as.data.frame(scores(fitasv, display = "vectors"))

plot(fit.tumourstatus, pch=5, col="black")

# make a nice ggplot of an ordination, use stat_ellipse for pretty but not very meaningful ellipses
# or add hulls, which isn't working right now due to plot_ordination putting geom_point in the wrong place
# hull_df <- as.data.frame(cbind(plotord[["data"]][["NMDS1"]], plotord[["data"]][["NMDS2"]]))
# colnames(hull_df) <- c("NMDS1", "NMDS2")

hull <- gethulls(hull_df, tse@colData@listData$location_simple, "NMDS1", "NMDS2")

plotord <- plot_ordination(physeq.tree, uwUF.ordu, type="sites", color = "location_simple") + 
  theme_bw() + 
  stat_ellipse() + 
  geom_point(aes(color = physeq.tree@sam_data$location_simple)) +
  ggtitle("Unweighted UniFrac")

plotord + #geom_segment(data = allasv.scrs[which(rownames(allasv.scrs) == "ASV_3619"),],
  #aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
  #arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  #geom_polygon(data = hull, aes(fill = groups), alpha = 0.1) +
  geom_text(inherit.aes = FALSE, data=as.data.frame(uwUF.ordu$points),aes(x=MDS1,y=MDS2,label=rownames(uwUF.ordu$points)))

actinobacillus <- allasv.scrs[c(2528, 2765, 3619),]

plotord + ggtitle("Unweighted UniFrac - Actinobacillus") +
  geom_segment(data = actinobacillus,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(inherit.aes = FALSE, data=as.data.frame(actinobacillus),aes(x=NMDS1,y=NMDS2,label=rownames(actinobacillus)))

taxtable_nona2 <- as.data.frame(taxtable_nona)

roseburia <- which(taxtable_nona2$Species == "Unclassified_Roseburia")
roseburia <- allasv.scrs[roseburia,]

plotord + ggtitle("Unweighted UniFrac - Roseburia") +
  geom_segment(data = roseburia,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(inherit.aes = FALSE, data=as.data.frame(roseburia),aes(x=NMDS1,y=NMDS2,label=rownames(roseburia)))

lachnospiraceae <- which(taxtable_nona2$Species == "Unclassified_Lachnospiraceae")
lachnospiraceae <- allasv.scrs[lachnospiraceae,]

plotord + ggtitle("Unweighted UniFrac - Lachnospiraceae") +
  geom_segment(data = lachnospiraceae,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(inherit.aes = FALSE, data=as.data.frame(lachnospiraceae),aes(x=NMDS1,y=NMDS2,label=rownames(lachnospiraceae)))

prevotella <- which(taxtable_nona2$Species == "Unclassified_Prevotella")
prevotella <- allasv.scrs[prevotella,]

plotord + ggtitle("Unweighted UniFrac - Prevotella") +
  geom_segment(data = prevotella,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(inherit.aes = FALSE, data=as.data.frame(prevotella),aes(x=NMDS1,y=NMDS2,label=rownames(prevotella)))

howardella <- which(taxtable_nona2$Species == "Unclassified_Howardella")
howardella <- allasv.scrs[howardella,]

plotord + ggtitle("Unweighted UniFrac - Howardella") +
  geom_segment(data = howardella,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(inherit.aes = FALSE, data=as.data.frame(howardella),aes(x=NMDS1,y=NMDS2,label=rownames(howardella)))

metabolites <- 
  
  
  
  
  clr_pca <- prcomp(t(tse_genus@assays@data@listData[["clr"]]))
df_out <- as.data.frame(clr_pca$x)

# construct df of all pc combinations, warning, if you have a lot this might die
allpcs <- as.data.frame(t(combn(colnames(df_out), 2)))

df_out$Disease <- as.factor(tse_genus@colData@listData$Disease)

pcx <- "PC1"
pcy <- "PC66"

p<-ggplot(df_out,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=Disease))
p<-p+geom_point()
p

# get convex hulls

hull <- df_out %>%
  group_by(Disease) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p + aes(fill = factor(Disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(x = paste(pcx, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcx, "PC")[[1]][2])][2])*100, "% variance explained"),
       y = paste(pcy, as.numeric(summary(clr_pca)[["importance"]][,as.numeric(strsplit(pcy, "PC")[[1]][2])][2])*100, "% variance explained"), fill = "Smoker")

################ TEST OVERLAPPING POLYGON BS

clr_pca <- prcomp(t(tse_genus@assays@data@listData[["clr"]]))
pcasummary <- data.frame(summary(clr_pca)[["importance"]])
df_out <- as.data.frame(clr_pca$x)

# construct df of all pc combinations, warning, if you have a lot this might die
allpcs <- as.data.frame(t(combn(colnames(df_out), 2)))

df_out$Disease <- as.factor(tse_genus@colData@listData$Disease)

# remove benign for ease of viz

df_out <- df_out[-which(df_out$Disease == "Benign"),]

p<-ggplot(df_out,aes(x=PC1,y=PC2,color=Disease))
p<-p+geom_point()
p

library(sp)

pcx <- "PC1"
pcy <- "PC2"
pcresults <- allpcs
pcresults$areaoverlap <- 0
pcresults$areaoverlappct <- 0
pcresults$areanooverlap <- 0
pcresults$areanooverlappct <- 0
pcresults$variancenooverlap <- 0
pcresults$pctpointsoutside <- 0

condition1 <- "Cancer"
condition2 <- "Normal"

for(i in 1:nrow(allpcs)){
  pcx <- allpcs[i,1]
  pcy <- allpcs[i,2]
  
  hull <- df_out %>%
    group_by(Disease) %>%
    dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))
  
  tempdf1 <- df_out[pcx]
  tempdf2 <- df_out[pcy]
  
  points1 <- point.in.polygon(data.frame(df_out[pcx])[which(df_out$disease == condition1),], data.frame(df_out[pcy])[which(df_out$disease == condition1),], data.frame(hull[pcx])[which(hull$disease == condition2),], data.frame(hull[pcy])[which(hull$disease == condition2),])
  points1df <- df_out[which(df_out$disease == condition1),]
  points1df <- points1df[which(points1 != 0),]
  
  points2 <- point.in.polygon(data.frame(df_out[pcx])[which(df_out$disease == condition2),], data.frame(df_out[pcy])[which(df_out$disease == condition2),], data.frame(hull[pcx])[which(hull$disease == condition1),], data.frame(hull[pcy])[which(hull$disease == condition1),])
  points2df <- df_out[which(df_out$disease == condition2),]
  points2df <- points2df[which(points2 != 0),]
  
  # poly_crossings needs filled polygons, duplicate first row as last row
  poly1 <- as.matrix(cbind(data.frame(hull[pcx])[which(hull$disease == condition1),], data.frame(hull[pcy])[which(hull$disease == condition1),]))
  poly2 <- as.matrix(cbind(data.frame(hull[pcx])[which(hull$disease == condition2),], data.frame(hull[pcy])[which(hull$disease == condition2),]))
  poly1 <- as.data.frame(poly1)
  poly2 <- as.data.frame(poly2)
  poly1[(nrow(poly1) +1),] <- poly1[1,]
  poly2[(nrow(poly2) +1),] <- poly2[1,]
  
  
  
  crossings <- poly_crossings(t(poly1), t(poly2))
  crossings <- as.data.frame(crossings)
  if(length(crossings) != 0){  
    colnames(crossings) <- c(pcx, pcy)
    crossings$disease <- "Overlap"
    
    overlappoints <- rbind(points1df[,which(colnames(points1df) %in% colnames(crossings))], points2df[,which(colnames(points2df) %in% colnames(crossings))], crossings)
    overlappoints$disease <- "Overlap"
    
    pcresults$pctpointsoutside[i] <- ((nrow(tempdf1) - nrow(overlappoints))/nrow(tempdf1))*100
    
    overlaphull <- overlappoints %>%
      dplyr::slice(chull(overlappoints))
    
    
    
    vec1 <- as.vector(overlaphull[,1])
    vec2 <- as.vector(overlaphull[,2])
    
    (diff(range(tempdf1)) * diff(range(tempdf2)))
    
    pcresults$areaoverlap[i] <- abs(polyarea(vec1, vec2))
    pcresults$areaoverlappct[i] <- (abs(polyarea(vec1, vec2)) / (diff(range(tempdf1)) * diff(range(tempdf2)))) *100
    
    poly1areaminus <- abs(polyarea(poly1[,1], poly1[,2])) - abs(polyarea(overlaphull[,1], overlaphull[,2]))
    poly2areaminus <- abs(polyarea(poly2[,1], poly2[,2])) - abs(polyarea(overlaphull[,1], overlaphull[,2]))
    pcresults$areanooverlap[i] <- (poly1areaminus + poly2areaminus)
    
    pcresults$areanooverlappct[i] <- (poly1areaminus + poly2areaminus) / (diff(range(tempdf1)) * diff(range(tempdf2))) * 100
    pcresults$variancenooverlap[i] <- ((pcasummary[2,pcx] + pcasummary[2,pcy])/100)*pcresults$areanooverlappct[i]
  }
  
  else
    
  {
    if(isTRUE(0 %in% points1 == FALSE)){
      pcresults$areaoverlap[i] <- paste(condition1, "100% within", condition2)
    }
    
    if(isTRUE(0 %in% points2 == FALSE)){
      pcresults$areaoverlap[i] <- paste(condition2, "100% within", condition1)
    }
    
  }
  
}

pcresults2 <- as.data.frame(matrix(ncol = 4, data = as.numeric(unlist(pcresults))))
pcresults2[,1:2] <- pcresults[,1:2]
pcresults2 <- pcresults2[-which(is.na(pcresults2[,3]) == TRUE),]

# get convex hulls

p<-ggplot(df_out,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=Disease))
p<-p+geom_point()
p

hull <- df_out %>%
  group_by(metadata$Disease[-which(metadata$Disease == "Benign")]) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p + aes(fill = factor(disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") + geom_polygon(data = overlaphull)

condition1 <- "Cancer"
condition2 <- "Normal"

points1 <- point.in.polygon(df_out$PC1[which(df_out$disease == condition1)], df_out$PC2[which(df_out$disease == condition1)], hull$PC1[which(hull$disease == condition2)], hull$PC2[which(hull$disease == condition2)])
points1df <- df_out[which(df_out$disease == condition1),]
points1df <- points1df[which(points1 != 0),]

points2 <- point.in.polygon(df_out$PC1[which(df_out$disease == condition2)], df_out$PC2[which(df_out$disease == condition2)], hull$PC1[which(hull$disease == condition1)], hull$PC2[which(hull$disease == condition1)])
points2df <- df_out[which(df_out$disease == condition2),]
points2df <- points2df[which(points2 != 0),]

# poly_crossings needs filled polygons, duplicate first row as last row
poly1 <- poly1 <- as.matrix(cbind(hull$PC1[which(hull$disease == condition1)], hull$PC2[which(hull$disease == condition1)]))
poly2 <- as.matrix(cbind(hull$PC1[which(hull$disease == condition2)], hull$PC2[which(hull$disease == condition2)]))
poly1 <- as.data.frame(poly1)
poly2 <- as.data.frame(poly2)
poly1[(nrow(poly1) +1),] <- poly1[1,]
poly2[(nrow(poly2) +1),] <- poly2[1,]

crossings <- poly_crossings(t(poly1), t(poly2))
crossings <- as.data.frame(crossings)

colnames(crossings) <- c("PC1", "PC2")
crossings$disease <- "Overlap"

overlappoints <- rbind(points1df[,which(colnames(points1df) %in% colnames(crossings))], points2df[,which(colnames(points2df) %in% colnames(crossings))], crossings)
overlappoints$disease <- "Overlap"

overlaphull <- overlappoints %>%
  dplyr::slice(chull(overlappoints))

p + aes(fill = factor(disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") + geom_polygon(data = overlaphull)

perfect_separation <- 
  
  # make into ellipses
  
  get_ellipse <- function(data, fill){
    edata <- as.matrix(data)
    ehull <- ellipsoidhull(edata)
    phull <- as.data.frame(predict(ehull))
    data.frame(
      x=phull$V1, 
      y=phull$y, 
      fill=rep(fill, nrow(phull))
    )
  }

ellipse1 <- get_ellipse(hull[1:length(which(hull$tstage == "1")),c(1,8)], fill = "1")
ellipse3 <- get_ellipse(hull[min(range(which(hull$tstage == 3))):max(range(which(hull$tstage == 3))),c(1,8)], fill = "3")
ellipse4 <- get_ellipse(hull[min(range(which(hull$tstage == 4))):max(range(which(hull$tstage == 4))),c(1,8)], fill = "4")

ellipse <- rbind(ellipse1, ellipse3, ellipse4)
p + geom_polygon(data=ellipse, aes(x=x, y=y, colour = fill), fill = NA)

df_out2 <- df_out
df_out2$tstage <- as.character(df_out2$tstage)
df_out2$tstage[which(df_out2$tstage != "1")] <- ">1"
df_out2$tstage <- as.factor(df_out2$tstage)

ellipse1 <- get_ellipse(hull[1:length(which(hull$tstage == "1")),c(1,8)], fill = "1")
ellipse2 <- get_ellipse(hull[min(range(which(hull$tstage !=1))):max(range(which(hull$tstage != 1))),c(1,8)], fill = ">1")
ellipse <- rbind(ellipse1, ellipse2)

p + geom_polygon(data=ellipse, aes(x=x, y=y, colour = fill), fill = NA)

# loadings in rotation, eyeball some pc1 associated taxa, check at genus level

View(clr_pca[["rotation"]])

abundance_analysis_data <- data.frame(t(assay(tse_genus, "clr")))

abundance_analysis_data <- cbind(
  abundance_analysis_data, 
  tstage = df_out2$tstage
)

conds

permanova_tstage <- vegan::adonis(t(tse_genus_prev30@assays@data$relabundance) ~ metadata$T,
                                  data = colData(tse_genus),
                                  permutations = 9999)

print(paste0("Different different cohorts and variance of abundance ",
             "between samples, p-value: ", 
             as.data.frame(permanova_tstage$aov.tab)["T", "Pr(>F)"]))

coefs <- coefficients(permanova_tstage)["metadata$T",]

top.coefs <- sort(head(coefs[rev(order(abs(coefs)))],20))

ggplot(data.frame(x = top.coefs,
                  y = factor(names(top.coefs),
                             unique(names(top.coefs)))),
       aes(x = x, y = y)) +
  geom_bar(stat="identity") +
  labs(x="",y="",title="Top Taxa") +
  theme_bw()

tse_prevotella <- tse_genus["Prevotella"]

genera <- names(abundance_analysis_data[, !names(abundance_analysis_data) %in% "tstage"])

wilcoxon_p <- c() # Initialize empty vector for p-values

# Do "for loop" over selected column names
for (i in genera) {
  
  result <- wilcox.test(abundance_analysis_data[, i] ~ tstage,
                        data = abundance_analysis_data)
  
  # Stores p-value to the vector with this column name
  wilcoxon_p[[i]]  <- result$p.value
  
}

wilcoxon_p <- data.frame(taxa =  names(wilcoxon_p),
                         p_raw = unlist(wilcoxon_p))

test <- seqtab %>% 
  filter_all(any_vars(. != 0))

set.seed(1)
testord <- metaMDS(test[1:10,], distance = 'bray', k = 2)

vegantestord <- vegan::vegdist(test[1:10,], method = 'bray')






sequences <- as.data.frame(rownames(seqtab))

seqtab <- seqtab[,colnames(seqtab) %in% metadata[,1]]

for(i in 1:4398){
  rownames(seqtab)[i] <- paste0("ASV_", i)
}

rownames(taxtable_nona) <- rownames(seqtab)
sequences[,2] <- rownames(seqtab)

se <- SummarizedExperiment(assays = list(counts = seqtab),
                           colData = metadata,
                           rowData = taxtable_nona)
se <- mia::estimateRichness(se, 
                            abund_values = "counts", 
                            index = "observed", 
                            name="observed")

# metadata uses numbers currently for location_2, convert to chr to make categorical for plots
se@colData@listData[["location_2"]] <- as.character(se@colData@listData[["location_2"]])


library(scater)
# scater demands object of class "SingleCellExperiment", assign class(se) <- "SingleCellExperiment" then change back if necessary
class(se) <- "SingleCellExperiment"

plotColData(se, 
            "observed", 
            "location_2", 
            colour_by = "location_2") +
  theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  ylab(expression(Richness[Observed]))














turdpolisher <- function(taxatable){
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
    outputseqtab <<- seqtab.nochim[,-c(addr)]
    print("Created new seqtab object outputseqtab without NA rows")
  }
}

plotColData <- function(object, y, x = NULL, 
                        colour_by = NULL, shape_by = NULL, size_by = NULL, 
                        by_exprs_values = "logcounts", other_fields=list(),
                        swap_rownames = NULL, ...)
{
  #if (!is(object, "SingleCellExperiment")) {
  #  stop("object must be an SingleCellExperiment object.")
  #}
  
  ## Define dataframe to pass to plotMetadata
  y_by_out <- retrieveCellInfo(object, y, search = "colData")
  y_lab <- y_by_out$name
  if (is.null(y_lab)) {
    stop(sprintf("could not find '%s' in 'colData(object)'", y))
  }
  df_to_plot <- data.frame(Y=y_by_out$val)
  
  if (!is.null(x)) {
    x_by_out <- retrieveCellInfo(object, x, search = "colData")
    x_lab <- x_by_out$name
    if (is.null(x_lab)) {
      stop(sprintf("could not find '%s' in 'rowData(object)'", x))
    }
    df_to_plot$X <- x_by_out$val
  } else {
    x_lab <- NULL
    df_to_plot$X <- factor(character(ncol(object)))
  }
  
  ## checking visualization arguments
  vis_out <- scater::.incorporate_common_vis_col(df_to_plot, se = object, 
                                                 colour_by = colour_by, shape_by = shape_by, size_by = size_by, 
                                                 by_exprs_values = by_exprs_values, other_fields = other_fields,
                                                 swap_rownames = swap_rownames)
  
  df_to_plot <- vis_out$df
  colour_by <- vis_out$colour_by
  shape_by <- vis_out$shape_by
  size_by <- vis_out$size_by
  
  # Creating the plot object:
  scater::.central_plotter(df_to_plot, xlab = x_lab, ylab = y_lab,
                           colour_by = colour_by, size_by = size_by, shape_by = shape_by, 
                           ..., point_FUN=NULL)
}

turdpolisher <- function(taxatable, seqtab){
  exit <- 0
  addr <- NULL
  tmpname <- NULL
  
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

















for(i in 1:nrow(allpcs)){
  pcx <- allpcs[i,1]
  pcy <- allpcs[i,2]
  
  hull <- df_out %>%
    group_by(metadata$Disease[-which(metadata$Disease == "Benign")]) %>%
    dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))
  
  tempdf1 <- df_out[pcx]
  tempdf2 <- df_out[pcy]
  
  points1 <- point.in.polygon(data.frame(df_out[pcx])[which(df_out$disease == condition1),], data.frame(df_out[pcy])[which(df_out$disease == condition1),], data.frame(hull[pcx])[which(hull$disease == condition2),], data.frame(hull[pcy])[which(hull$disease == condition2),])
  points1df <- df_out[which(df_out$disease == condition1),]
  points1df <- points1df[which(points1 != 0),]
  
  points2 <- point.in.polygon(data.frame(df_out[pcx])[which(df_out$disease == condition2),], data.frame(df_out[pcy])[which(df_out$disease == condition2),], data.frame(hull[pcx])[which(hull$disease == condition1),], data.frame(hull[pcy])[which(hull$disease == condition1),])
  points2df <- df_out[which(df_out$disease == condition2),]
  points2df <- points2df[which(points2 != 0),]
  
  # poly_crossings needs filled polygons, duplicate first row as last row
  poly1 <- as.matrix(cbind(data.frame(hull[pcx])[which(hull$disease == condition1),], data.frame(hull[pcy])[which(hull$disease == condition1),]))
  poly2 <- as.matrix(cbind(data.frame(hull[pcx])[which(hull$disease == condition2),], data.frame(hull[pcy])[which(hull$disease == condition2),]))
  poly1 <- as.data.frame(poly1)
  poly2 <- as.data.frame(poly2)
  poly1[(nrow(poly1) +1),] <- poly1[1,]
  poly2[(nrow(poly2) +1),] <- poly2[1,]
  
  
  
  crossings <- poly_crossings(t(poly1), t(poly2))
  crossings <- as.data.frame(crossings)
  if(length(crossings) != 0){  
    colnames(crossings) <- c(pcx, pcy)
    crossings$disease <- "Overlap"
    
    overlappoints <- rbind(points1df[,which(colnames(points1df) %in% colnames(crossings))], points2df[,which(colnames(points2df) %in% colnames(crossings))], crossings)
    overlappoints$disease <- "Overlap"
    
    pcresults$pctpointsoutside[i] <- ((nrow(tempdf1) - nrow(overlappoints))/nrow(tempdf1))*100
    
    overlaphull <- overlappoints %>%
      dplyr::slice(chull(overlappoints))
    
    
    
    vec1 <- as.vector(overlaphull[,1])
    vec2 <- as.vector(overlaphull[,2])
    
    (diff(range(tempdf1)) * diff(range(tempdf2)))
    
    pcresults$areaoverlap[i] <- abs(polyarea(vec1, vec2))
    pcresults$areaoverlappct[i] <- (abs(polyarea(vec1, vec2)) / (diff(range(tempdf1)) * diff(range(tempdf2)))) *100
    
    poly1areaminus <- abs(polyarea(poly1[,1], poly1[,2])) - abs(polyarea(overlaphull[,1], overlaphull[,2]))
    poly2areaminus <- abs(polyarea(poly2[,1], poly2[,2])) - abs(polyarea(overlaphull[,1], overlaphull[,2]))
    pcresults$areanooverlap[i] <- (poly1areaminus + poly2areaminus)
    
    pcresults$areanooverlappct[i] <- (poly1areaminus + poly2areaminus) / (diff(range(tempdf1)) * diff(range(tempdf2))) * 100
    pcresults$variancenooverlap[i] <- ((pcasummary[2,pcx] + pcasummary[2,pcy])/100)*pcresults$areanooverlappct[i]
  }
  
  else
    
  {
    if(isTRUE(0 %in% points1 == FALSE)){
      pcresults$areaoverlap[i] <- paste(condition1, "100% within", condition2)
    }
    
    if(isTRUE(0 %in% points2 == FALSE)){
      pcresults$areaoverlap[i] <- paste(condition2, "100% within", condition1)
    }
    
  }
  
}

pcresults2 <- as.data.frame(matrix(ncol = 4, data = as.numeric(unlist(pcresults))))
pcresults2[,1:2] <- pcresults[,1:2]
pcresults2 <- pcresults2[-which(is.na(pcresults2[,3]) == TRUE),]

# get convex hulls

p<-ggplot(df_out,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=disease))
p<-p+geom_point()
p

hull <- df_out %>%
  group_by(metadata$Disease[-which(metadata$Disease == "Benign")]) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p + aes(fill = factor(disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") + geom_polygon(data = overlaphull)

condition1 <- "Cancer"
condition2 <- "Normal"

points1 <- point.in.polygon(df_out$PC1[which(df_out$disease == condition1)], df_out$PC2[which(df_out$disease == condition1)], hull$PC1[which(hull$disease == condition2)], hull$PC2[which(hull$disease == condition2)])
points1df <- df_out[which(df_out$disease == condition1),]
points1df <- points1df[which(points1 != 0),]

points2 <- point.in.polygon(df_out$PC1[which(df_out$disease == condition2)], df_out$PC2[which(df_out$disease == condition2)], hull$PC1[which(hull$disease == condition1)], hull$PC2[which(hull$disease == condition1)])
points2df <- df_out[which(df_out$disease == condition2),]
points2df <- points2df[which(points2 != 0),]

# poly_crossings needs filled polygons, duplicate first row as last row
poly1 <- poly1 <- as.matrix(cbind(hull$PC1[which(hull$disease == condition1)], hull$PC2[which(hull$disease == condition1)]))
poly2 <- as.matrix(cbind(hull$PC1[which(hull$disease == condition2)], hull$PC2[which(hull$disease == condition2)]))
poly1 <- as.data.frame(poly1)
poly2 <- as.data.frame(poly2)
poly1[(nrow(poly1) +1),] <- poly1[1,]
poly2[(nrow(poly2) +1),] <- poly2[1,]

crossings <- poly_crossings(t(poly1), t(poly2))
crossings <- as.data.frame(crossings)

colnames(crossings) <- c("PC1", "PC2")
crossings$disease <- "Overlap"

overlappoints <- rbind(points1df[,which(colnames(points1df) %in% colnames(crossings))], points2df[,which(colnames(points2df) %in% colnames(crossings))], crossings)
overlappoints$disease <- "Overlap"

overlaphull <- overlappoints %>%
  dplyr::slice(chull(overlappoints))

p + aes(fill = factor(disease)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") + geom_polygon(data = overlaphull)

perfect_separation <- 