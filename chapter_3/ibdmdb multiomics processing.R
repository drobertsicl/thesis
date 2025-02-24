library(tidyverse)
library(mia)

# process biopsy 16s taxonomic_profile.tsv file to regular taxonomy table

taxonomic_profiles <- read_tsv(file = "C:/Users/eatmo/Desktop/disease duration paper/taxonomic_profiles.txt")

taxtable2 <- data.frame(matrix(data = NA, ncol = 6, nrow = nrow(taxonomic_profiles)))
for(i in 1:nrow(taxonomic_profiles)){
  taxtable2[i,1:6] <- t(data.frame(c(strsplit(taxonomic_profiles$taxonomy[i], "\\;"))))
}

for (i in 1:nrow(taxtable2)){
  taxtable2[i,2] <- data.frame(strsplit(taxtable2[i,2], " __"))[2,]
  taxtable2[i,3] <- data.frame(strsplit(taxtable2[i,3], " __"))[2,]
  taxtable2[i,4] <- data.frame(strsplit(taxtable2[i,4], " __"))[2,]
  taxtable2[i,5] <- data.frame(strsplit(taxtable2[i,5], " __"))[2,]
  taxtable2[i,6] <- data.frame(strsplit(taxtable2[i,6], " __"))[2,]
}


for(i in 1:length(which(taxtable2[,3] == "c"))){
  taxtable2[which(taxtable2[,3] == "c"),][1,] <- c(taxtable2[which(taxtable2[,3] == "c"),][1,1:2], rep(paste0("Unclassified_", taxtable2[which(taxtable2[,3] == "c"),2][1]), (6-2)))
}


for(i in 1:length(which(taxtable2[,4] == "o"))){
  taxtable2[which(taxtable2[,4] == "o"),][1,] <- c(taxtable2[which(taxtable2[,4] == "o"),][1,1:3], rep(paste0("Unclassified_", taxtable2[which(taxtable2[,4] == "o"),3][1]), (6-3)))
}

for(i in 1:length(which(taxtable2[,5] == "f"))){
  taxtable2[which(taxtable2[,5] == "f"),][1,] <- c(taxtable2[which(taxtable2[,5] == "f"),][1,1:4], rep(paste0("Unclassified_", taxtable2[which(taxtable2[,5] == "f"),4][1]), (6-4)))
}

for(i in 1:length(which(taxtable2[,6] == "g"))){
  taxtable2[which(taxtable2[,6] == "g"),][1,] <- c(taxtable2[which(taxtable2[,6] == "g"),][1,1:5], rep(paste0("Unclassified_", taxtable2[which(taxtable2[,6] == "g"),5][1]), (6-5)))
}

colnames(taxtable2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

write_csv(taxtable2, file = "C:/Users/eatmo/Desktop/disease duration paper/biopsy_taxtable.csv")

taxonomic_profiles <- taxonomic_profiles[,-1]

for(i in 1:nrow(taxonomic_profiles)){
  rownames(taxonomic_profiles)[i] <- paste0("ASV_", i)
}

write_csv(taxonomic_profiles, file = "C:/Users/eatmo/Desktop/disease duration paper/biopsy_seqtable.csv")

# process shotgun samples

# load as shotgun_x, remove ncbi column, set column to rownames for later merging
temprownames <- NULL

# THIS IS INSANELY STUPID

path <- "C:/Users/eatmo/Desktop/disease duration paper/stool_shotgun/"

setwd(path)

files <- list.files(path)

for(i in 1:length(files)){
  assign(paste0("shotgun_", i), read_tsv(files[i], skip = 3, col_select = c("#clade_name", "relative_abundance")))
  assign(paste0("shotgun_", i), column_to_rownames(get(paste0("shotgun_", i)), var = "#clade_name"))
  temprownames <<- rownames(get(paste0("shotgun_", i)))
  assign(paste0("shotgun_", i), data.frame(t(data.frame(strsplit(get(paste0("shotgun_", i))$relative_abundance, "\t#")))[,1]))
  assign(paste0("shotgun_", i), remove_rownames(cbind(temprownames, get(paste0("shotgun_", i)))))
  #rm(temprownames)
  #assign(paste0("shotgun_", i), column_to_rownames(get(paste0("shotgun_", i)), var = "temprownames"))
  assign(paste0("shotgun_", i), `colnames<-`(get(paste0("shotgun_", i)), c("temprownames","relative_abundance")))
}

for(i in 1:(length(files) -1)){
  if (i == 1) {
    test2 <- merge(get(paste0("shotgun_", i)), get(paste0("shotgun_", i+1)), by = "temprownames", all = TRUE)
    rownames(test2) <- test2[,1]
    #test2 <- test2[,-1]
  } else {
    test2 <- merge(test2, get(paste0("shotgun_", i+1)), by = "temprownames", all = TRUE)
    rownames(test2) <- test2[,1]
    #test2 <- test2[,-1]
  }
}

test2[is.na(test2)] <- 0
test2 <- test2[,-1]
colnames(test2) <- t(data.frame(strsplit(files, "_taxonomic_profile")))[,1]

taxtable2 <- data.frame(matrix(data = NA, ncol = 8, nrow = nrow(test2)))
for(i in 1:nrow(test2)){
  taxtable2[i,1:8] <- c(t(as.data.frame(strsplit(rownames(test2)[i], "\\|"))), rep("NA", 8-length(t(as.data.frame(strsplit(rownames(test2)[i], "\\|"))))))
}

fixTable <- function(input_table) {
  tempRownames <- rownames(input_table)
  tempColnames <- colnames(input_table)
  input_table <- data.frame(matrix(data = as.numeric(unlist(input_table)), nrow = nrow(input_table), ncol = ncol(input_table)))
  rownames(input_table) <- tempRownames
  colnames(input_table) <- tempColnames
  return(input_table)
}

genustab <- test2[which(taxtable2[,6] != "NA" & taxtable2[,7] == "NA"),]
genustab <- fixTable(genustab)

speciestab <- test2[which(taxtable2[,7] != "NA" & taxtable2[,8] == "NA"),]
speciestab <- fixTable(speciestab)

# check that column sums ~100%

sumsgenus <- NULL 
for(i in 1:ncol(genustab)){sumsgenus[i] <- sum(genustab[,i])}
sumsgenus

sumsspecies <- NULL 
for(i in 1:ncol(speciestab)){sumsspecies[i] <- sum(speciestab[,i])}
sumsspecies

genustax <- taxtable2[which(taxtable2[,6] != "NA" & taxtable2[,7] == "NA"),-c(7,8)]
speciestax <- taxtable2[which(taxtable2[,7] != "NA" & taxtable2[,8] == "NA"),-8]

colnames(genustax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
colnames(speciestax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

for (i in 1:nrow(genustax)){
  genustax[i,1] <- data.frame(strsplit(genustax[i,1], "__"))[2,]
  genustax[i,2] <- data.frame(strsplit(genustax[i,2], "__"))[2,]
  genustax[i,3] <- data.frame(strsplit(genustax[i,3], "__"))[2,]
  genustax[i,4] <- data.frame(strsplit(genustax[i,4], "__"))[2,]
  genustax[i,5] <- data.frame(strsplit(genustax[i,5], "__"))[2,]
  genustax[i,6] <- data.frame(strsplit(genustax[i,6], "__"))[2,]
}

for (i in 1:nrow(speciestax)){
  speciestax[i,1] <- data.frame(strsplit(speciestax[i,1], "__"))[2,]
  speciestax[i,2] <- data.frame(strsplit(speciestax[i,2], "__"))[2,]
  speciestax[i,3] <- data.frame(strsplit(speciestax[i,3], "__"))[2,]
  speciestax[i,4] <- data.frame(strsplit(speciestax[i,4], "__"))[2,]
  speciestax[i,5] <- data.frame(strsplit(speciestax[i,5], "__"))[2,]
  speciestax[i,6] <- data.frame(strsplit(speciestax[i,6], "__"))[2,]
  speciestax[i,7] <- data.frame(strsplit(speciestax[i,7], "__"))[2,]
}

rm(list=paste0("shotgun_", 1:length(files)))

##########################
# load metabolomics data #
##########################

# get raw filenames from metadata

mbx_metadata <- integrated_metadata[-which(is.na(integrated_metadata$raw_filename)), ]

c8pos <- read_csv("C:/Users/eatmo/Desktop/disease duration paper/raw ms/C8p_NN.csv")
c8pos <- cbind(c8pos[,1:6], c8pos[,as.vector(unlist(mbx_metadata[,6]))])
write_csv(c8pos, file = "C:/Users/eatmo/Desktop/disease duration paper/raw ms/c8p_onlysamples.csv")


c18neg <- read_csv("C:/Users/eatmo/Desktop/disease duration paper/raw ms/C18n_NN.csv")
c18neg <- cbind(c18neg[,1:6], c18neg[,as.vector(unlist(mbx_metadata[,6]))])
write_csv(c18neg, file = "C:/Users/eatmo/Desktop/disease duration paper/raw ms/c18n_onlysamples.csv")

hilicpos <- read_csv("C:/Users/eatmo/Desktop/disease duration paper/raw ms/HILp_NN.csv")
hilicpos <- cbind(hilicpos[,1:6], hilicpos[,as.vector(unlist(mbx_metadata[,6]))])
write_csv(hilicpos, file = "C:/Users/eatmo/Desktop/disease duration paper/raw ms/HILp_onlysamples.csv")

hilicneg <- read_csv("C:/Users/eatmo/Desktop/disease duration paper/raw ms/HILn_NN.csv")
hilicneg <- cbind(hilicneg[,1:6], hilicneg[,as.vector(unlist(mbx_metadata[,6]))])
write_csv(hilicneg, file = "C:/Users/eatmo/Desktop/disease duration paper/raw ms/HILn_onlysamples.csv")

# make subsets of only the numeric data, overwrite into original table later

c8pos_mat <- c8pos[,7:29]
c18neg_mat <- c18neg[,7:29]
hilicpos_mat <- hilicpos[,7:29]
hilicneg_mat <- hilicneg[,7:29]

# remove features with >50% missing values

over50 <- NULL
for(i in 1:nrow(c8pos_mat)){
  over50 <- c(over50, length(which(is.na(c8pos_mat[i,]))))
}

c8pos_mat <- c8pos_mat[-which(over50 > 11),]
c8pos <- c8pos[-which(over50 > 11),]

over50 <- NULL
for(i in 1:nrow(c18neg_mat)){
  over50 <- c(over50, length(which(is.na(c18neg_mat[i,]))))
}

c18neg_mat <- c18neg_mat[-which(over50 > 11),]
c18neg <- c18neg[-which(over50 > 11),]

over50 <- NULL
for(i in 1:nrow(hilicpos_mat)){
  over50 <- c(over50, length(which(is.na(hilicpos_mat[i,]))))
}

hilicpos_mat <- hilicpos_mat[-which(over50 > 11),]
hilicpos <- hilicpos[-which(over50 > 11),]

over50 <- NULL
for(i in 1:nrow(hilicneg_mat)){
  over50 <- c(over50, length(which(is.na(hilicneg_mat[i,]))))
}

hilicneg_mat <- hilicneg_mat[-which(over50 > 11),]
hilicneg <- hilicneg[-which(over50 > 11),]

# replace remaining NA values with 1/5th of lowest feature value
lowval <- NULL
for(i in 1:nrow(c8pos_mat)){
  if(length(which(is.na(c8pos_mat[i,]))) > 0){
    lowval <- c(lowval, (min(as.numeric(c8pos_mat[i,])[-which(is.na(c8pos_mat[i,]))])/5))
  }
  else
  {
    lowval <- c(lowval, NA)
  }
}

# find minimum non NA value
min(lowval[-which(is.na(lowval))])

# add offset of 1 to allow log normalisation if this is 0

c8pos_mat[!is.na(c8pos_mat)] <- c8pos_mat[!is.na(c8pos_mat)] +1

# scale data by 5 so that adding 1/5th of lowest value will yield 1, which logs to 0

c8pos_mat[!is.na(c8pos_mat)] <- c8pos_mat[!is.na(c8pos_mat)] *5

# re run lowval code

lowval <- NULL
for(i in 1:nrow(c8pos_mat)){
  if(length(which(is.na(c8pos_mat[i,]))) > 0){
    lowval <- c(lowval, (min(as.numeric(c8pos_mat[i,])[-which(is.na(c8pos_mat[i,]))])/5))
  }
  else
  {
    lowval <- c(lowval, NA)
  }
}

# find minimum non NA value
min(lowval[-which(is.na(lowval))])

# pass to metaboanalyst to use built in filtering etc

c8pos_mat <- c8pos_mat[,mbx_metadata$raw_filename]
c8pos_mat <- cbind(c8pos$Compound, c8pos_mat)
colnames(c8pos_mat[1]) <- "Compound"
c8pos_mat <- rbind(c("", mbx_metadata$over10), c8pos_mat)

write_csv(c8pos_mat, file = "C:/Users/eatmo/Desktop/disease duration paper/preprocessed_c8.csv")

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "C:/Users/eatmo/Desktop/disease duration paper/preprocessed_c8.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)

# remove features with >50% missing, add offset, log transform, inspect distribution, export processed data via qs

mSet<-ReplaceMin(mSet)

# don't use variable filter, we don't need to remove thousands of features

mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

# play with normalisation/scaling/etc and check working directory for exported plots
mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet <- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

# extract processed metaboanalyst table
int.mat <- qs::qread(file = "complete_norm.qs")

# check for differences in order of features
"FALSE" %in% (colnames(int.mat) == c18neg_mat$Compound[2:nrow(c18neg_mat)])

# check fold change etc, likely nothing significant in univariate

mSet <- FC.Anal(mSet, 2.0, 0, paired = FALSE)
mSet <- PlotFC(mSet, "fc_0_", "png", 72, width=NA)

mSet <- Ttests.Anal(mSet, F, 0.05, FALSE, TRUE)
mSet <- PlotTT(mSet, "tt_0_", "png", 72, width=NA)

# check metaboanalyst PCA

mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
mSet <- PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
mSet <- PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
mSet <- PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)

# check pca pairs plot to see if anything of interest pops out

c8_pca <- prcomp(int.mat)
c8_scores <- data.frame(c8_pca$x)


# rearrange columns to match metadata order
c8_scores <- c8_scores[mbx_metadata$raw_filename,]
c8_scores$varofinterest <- mbx_metadata$over10

# get convex hulls

pcx <- "PC2"
pcy <- "PC3"

hull <- c8_scores %>%
  group_by(varofinterest) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p<-ggplot(c8_scores,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=varofinterest))
p<-p+geom_point()
#p<-p+geom_line(mapping = test)
p

p + aes(fill = factor(varofinterest)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none")

# transpose
int.mat <- t(data.frame(int.mat))
# change order
int.mat <- int.mat[,mbx_metadata$raw_filename]
# overwrite matrix area in duplicate c8 dataframe

# check sample order is the same
colnames(c8pos[7:ncol(c8pos)]) == colnames(int.mat)
# check feature order is the same
"FALSE" %in% (rownames(int.mat) == c8pos$Compound)

c8pos_processed <- cbind(c8pos[,1:6], int.mat)
write_csv(c8pos_processed, file = "C:/Users/eatmo/Desktop/disease duration paper/c8pos_processed_scaledx5.csv")

# negative mode

lowval <- NULL
for(i in 1:nrow(c18neg_mat)){
  if(length(which(is.na(c18neg_mat[i,]))) > 0){
    lowval <- c(lowval, (min(as.numeric(c18neg_mat[i,])[-which(is.na(c18neg_mat[i,]))])/5))
  }
  else
  {
    lowval <- c(lowval, NA)
  }
}

# find minimum non NA value
min(lowval[-which(is.na(lowval))])

# add offset of 1 in this case to maintain the same preprocessing?
# or should data be scaled such that 1/5th of smallest value is 1?

c18neg_mat[!is.na(c18neg_mat)] <- c18neg_mat[!is.na(c18neg_mat)] +1

# scale data by 5 so that adding 1/5th of values close to 1 won't go negative

c18neg_mat[!is.na(c18neg_mat)] <- c18neg_mat[!is.na(c18neg_mat)] *5

lowval <- NULL
for(i in 1:nrow(c18neg_mat)){
  if(length(which(is.na(c18neg_mat[i,]))) > 0){
    lowval <- c(lowval, (min(as.numeric(c18neg_mat[i,])[-which(is.na(c18neg_mat[i,]))])/5))
  }
  else
  {
    lowval <- c(lowval, NA)
  }
}

# find minimum non NA value
min(lowval[-which(is.na(lowval))])

#construct pre-processed data
c18neg_mat <- c18neg_mat[,mbx_metadata$raw_filename]
c18neg_mat <- cbind(c18neg$Compound, c18neg_mat)
colnames(c18neg_mat)[1] <- "Compound"
c18neg_mat <- rbind(c("", mbx_metadata$over10), c18neg_mat)

write_csv(c18neg_mat, file = "C:/Users/eatmo/Desktop/disease duration paper/preprocessed_c18.csv")

# handoff to metaboanalyst, remove mSet first so it doesn't go insane

rm(mSet)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "C:/Users/eatmo/Desktop/disease duration paper/preprocessed_c18.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)

# remove features with >50% missing, add offset, log transform, inspect distribution, export processed data via qs

mSet<-ReplaceMin(mSet)

# don't use variable filter, we don't need to remove thousands of features

mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

# extract processed metaboanalyst table
int.mat <- qs::qread(file = "complete_norm.qs")

# remove features taken out by metaboanalyst
#c18neg_mat <- c18neg_mat[which(c18neg_mat$Compound %in% colnames(int.mat)),]
# re add feature label row
#c18neg_mat <- rbind(c("", mbx_metadata$over10), c18neg_mat)

# check for differences in order of features
"FALSE" %in% (colnames(int.mat) == c18neg_mat$Compound[2:nrow(c18neg_mat)])

# check fold change etc, likely nothing significant in univariate

mSet <- FC.Anal(mSet, 2.0, 0, paired = FALSE)
mSet <- PlotFC(mSet, "fc_0_", "png", 72, width=NA)

mSet <- Ttests.Anal(mSet, F, 0.05, FALSE, TRUE)
mSet <- PlotTT(mSet, "tt_0_", "png", 72, width=NA)

# check metaboanalyst PCA

mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
mSet <- PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
mSet <- PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
mSet <- PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)

# check pca pairs plot to see if anything of interest pops out

c18_pca <- prcomp(int.mat)
c18_scores <- data.frame(c18_pca$x)


# rearrange columns to match metadata order
c18_scores <- c18_scores[mbx_metadata$raw_filename,]
c18_scores$varofinterest <- mbx_metadata$over10

# get convex hulls

pcx <- "PC3"
pcy <- "PC4"

hull <- c18_scores %>%
  group_by(varofinterest) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p<-ggplot(c18_scores,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=varofinterest))
p<-p+geom_point()
#p<-p+geom_line(mapping = test)
p

p + aes(fill = factor(varofinterest)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none")

# transpose
int.mat <- t(data.frame(int.mat))
# change order
int.mat <- int.mat[,mbx_metadata$raw_filename]
# overwrite matrix area in duplicate c8 dataframe

# check sample order is the same
colnames(c18neg[7:ncol(c18neg)]) == colnames(int.mat)
# check feature order is the same
"FALSE" %in% (rownames(int.mat) == c18neg$Compound)

c18neg_processed <- cbind(c18neg[,1:6], int.mat)
write_csv(c18neg_processed, file = "C:/Users/eatmo/Desktop/disease duration paper/c18neg_processed_scaledx5.csv")

# hilic pos

lowval <- NULL
for(i in 1:nrow(hilicpos_mat)){
  if(length(which(is.na(hilicpos_mat[i,]))) > 0){
    lowval <- c(lowval, (min(as.numeric(hilicpos_mat[i,])[-which(is.na(hilicpos_mat[i,]))])/5))
  }
  else
  {
    lowval <- c(lowval, NA)
  }
}

# find minimum non NA value
min(lowval[-which(is.na(lowval))])

# add offset of 1 to allow log normalisation if this is 0

hilicpos_mat[!is.na(hilicpos_mat)] <- hilicpos_mat[!is.na(hilicpos_mat)] +1

# scale data by 5 so that adding 1/5th of lowest value will yield 1, which logs to 0

hilicpos_mat[!is.na(hilicpos_mat)] <- hilicpos_mat[!is.na(hilicpos_mat)] *5

lowval <- NULL
for(i in 1:nrow(hilicpos_mat)){
  if(length(which(is.na(hilicpos_mat[i,]))) > 0){
    lowval <- c(lowval, (min(as.numeric(hilicpos_mat[i,])[-which(is.na(hilicpos_mat[i,]))])/5))
  }
  else
  {
    lowval <- c(lowval, NA)
  }
}

# find minimum non NA value
min(lowval[-which(is.na(lowval))])

#construct pre-processed data
hilicpos_mat <- hilicpos_mat[,mbx_metadata$raw_filename]
hilicpos_mat <- cbind(hilicpos$Compound, hilicpos_mat)
colnames(hilicpos_mat)[1] <- "Compound"
hilicpos_mat <- rbind(c("", mbx_metadata$over10), hilicpos_mat)

write_csv(hilicpos_mat, file = "C:/Users/eatmo/Desktop/disease duration paper/preprocessed_hilicpos.csv")

# handoff to metaboanalyst, remove mSet first so it doesn't go insane

rm(mSet)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "C:/Users/eatmo/Desktop/disease duration paper/preprocessed_hilicpos.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)

# remove features with >50% missing, add offset, log transform, inspect distribution, export processed data via qs

mSet<-ReplaceMin(mSet)

# don't use variable filter, we don't need to remove thousands of features

mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

# play with normalisation/scaling/etc and check working directory for exported plots
mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet <- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

# extract processed metaboanalyst table
int.mat <- qs::qread(file = "complete_norm.qs")

# remove features taken out by metaboanalyst
#c18neg_mat <- c18neg_mat[which(c18neg_mat$Compound %in% colnames(int.mat)),]
# re add feature label row
#c18neg_mat <- rbind(c("", mbx_metadata$over10), c18neg_mat)

# check for differences in order of features
"FALSE" %in% (colnames(int.mat) == hilicpos_mat$Compound[2:nrow(hilicpos_mat)])

# check fold change etc, likely nothing significant in univariate

mSet <- FC.Anal(mSet, 2.0, 0, paired = FALSE)
mSet <- PlotFC(mSet, "fc_0_", "png", 72, width=NA)

mSet <- Ttests.Anal(mSet, F, 0.05, FALSE, TRUE)
mSet <- PlotTT(mSet, "tt_0_", "png", 72, width=NA)

# check metaboanalyst PCA

mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
mSet <- PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
mSet <- PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
mSet <- PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)

# check pca pairs plot to see if anything of interest pops out

hilicpos_pca <- prcomp(int.mat)
hilicpos_scores <- data.frame(hilicpos_pca$x)


# rearrange columns to match metadata order
hilicpos_scores <- hilicpos_scores[mbx_metadata$raw_filename,]
hilicpos_scores$varofinterest <- mbx_metadata$over10

# get convex hulls

pcx <- "PC3"
pcy <- "PC5"

hull <- hilicpos_scores %>%
  group_by(varofinterest) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p<-ggplot(hilicpos_scores,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=varofinterest))
p<-p+geom_point()
#p<-p+geom_line(mapping = test)
p

p + aes(fill = factor(varofinterest)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none")

# transpose
int.mat <- t(data.frame(int.mat))
# change order
int.mat <- int.mat[,mbx_metadata$raw_filename]
# overwrite matrix area in duplicate c8 dataframe

# check sample order is the same
colnames(hilicpos[7:ncol(hilicpos)]) == colnames(int.mat)
# check feature order is the same
"FALSE" %in% (rownames(int.mat) == hilicpos$Compound)

hilicpos_processed <- cbind(hilicpos[,1:6], int.mat)
write_csv(hilicpos_processed, file = "C:/Users/eatmo/Desktop/disease duration paper/hilicpos_processed_scaledx5.csv")

# hilic neg

lowval <- NULL
for(i in 1:nrow(hilicneg_mat)){
  if(length(which(is.na(hilicneg_mat[i,]))) > 0){
    lowval <- c(lowval, (min(as.numeric(hilicneg_mat[i,])[-which(is.na(hilicneg_mat[i,]))])/5))
  }
  else
  {
    lowval <- c(lowval, NA)
  }
}

# find minimum non NA value
min(lowval[-which(is.na(lowval))])

# add offset of 1 to allow log normalisation if this is 0

hilicneg_mat[!is.na(hilicneg_mat)] <- hilicneg_mat[!is.na(hilicneg_mat)] +1

# scale data by 5 so that adding 1/5th of lowest value will yield 1, which logs to 0

hilicneg_mat[!is.na(hilicneg_mat)] <- hilicneg_mat[!is.na(hilicneg_mat)] *5

lowval <- NULL
for(i in 1:nrow(hilicneg_mat)){
  if(length(which(is.na(hilicneg_mat[i,]))) > 0){
    lowval <- c(lowval, (min(as.numeric(hilicneg_mat[i,])[-which(is.na(hilicneg_mat[i,]))])/5))
  }
  else
  {
    lowval <- c(lowval, NA)
  }
}

# find minimum non NA value
min(lowval[-which(is.na(lowval))])

#construct pre-processed data
hilicneg_mat <- hilicneg_mat[,mbx_metadata$raw_filename]
hilicneg_mat <- cbind(hilicneg$Compound, hilicneg_mat)
colnames(hilicneg_mat)[1] <- "Compound"
hilicneg_mat <- rbind(c("", mbx_metadata$over10), hilicneg_mat)

write_csv(hilicneg_mat, file = "C:/Users/eatmo/Desktop/disease duration paper/preprocessed_hilicneg.csv")

# handoff to metaboanalyst, remove mSet first so it doesn't go insane

rm(mSet)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "C:/Users/eatmo/Desktop/disease duration paper/preprocessed_hilicneg.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)

# remove features with >50% missing, add offset, log transform, inspect distribution, export processed data via qs

mSet<-ReplaceMin(mSet)

# don't use variable filter, we don't need to remove thousands of features

mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

# play with normalisation/scaling/etc and check working directory for exported plots
mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet <- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

# extract processed metaboanalyst table
int.mat <- qs::qread(file = "complete_norm.qs")

# remove features taken out by metaboanalyst
#c18neg_mat <- c18neg_mat[which(c18neg_mat$Compound %in% colnames(int.mat)),]
# re add feature label row
#c18neg_mat <- rbind(c("", mbx_metadata$over10), c18neg_mat)

# check for differences in order of features
"FALSE" %in% (colnames(int.mat) == hilicneg_mat$Compound[2:nrow(hilicneg_mat)])

# check fold change etc, likely nothing significant in univariate

mSet <- FC.Anal(mSet, 2.0, 0, paired = FALSE)
mSet <- PlotFC(mSet, "fc_0_", "png", 72, width=NA)

mSet <- Ttests.Anal(mSet, F, 0.05, FALSE, TRUE)
mSet <- PlotTT(mSet, "tt_0_", "png", 72, width=NA)

# check metaboanalyst PCA

mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
mSet <- PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
mSet <- PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
mSet <- PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)

# check pca pairs plot to see if anything of interest pops out

hilicneg_pca <- prcomp(int.mat)
hilicneg_scores <- data.frame(hilicneg_pca$x)


# rearrange columns to match metadata order
hilicneg_scores <- hilicneg_scores[mbx_metadata$raw_filename,]
hilicneg_scores$varofinterest <- mbx_metadata$over10

# get convex hulls

pcx <- "PC3"
pcy <- "PC5"

hull <- hilicneg_scores %>%
  group_by(varofinterest) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p<-ggplot(hilicneg_scores,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=varofinterest))
p<-p+geom_point()
#p<-p+geom_line(mapping = test)
p

p + aes(fill = factor(varofinterest)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none")

# transpose
int.mat <- t(data.frame(int.mat))
# change order
int.mat <- int.mat[,mbx_metadata$raw_filename]
# overwrite matrix area in duplicate c8 dataframe

# check sample order is the same
colnames(hilicneg[7:ncol(hilicneg)]) == colnames(int.mat)
# check feature order is the same
"FALSE" %in% (rownames(int.mat) == hilicneg$Compound)

hilicneg_processed <- cbind(hilicneg[,1:6], int.mat)
write_csv(hilicneg_processed, file = "C:/Users/eatmo/Desktop/disease duration paper/hilicneg_processed_scaledx5.csv")

# test multitable pca scaling

# get scaling factor to equalise PC1 between hilic datasets
svd(hilicpos_processed[,7:29])$d[1]/svd(hilicneg_processed[,7:29])$d[1]

# multiply lower value dataset by higher value
hilicneg_processed2 <- hilicneg_processed
hilicneg_processed2[,7:29] <- hilicneg_processed2[,7:29] * (svd(hilicpos_processed[,7:29])$d[1]/svd(hilicneg_processed[,7:29])$d[1])

# check svd
svd(hilicneg_processed2[,7:29])$d[1]

# rbind together, then check pca

test <- rbind(hilicpos_processed[,7:29], hilicneg_processed2[,7:29])

test_pca <- prcomp(t(test), scale = FALSE)

test_scores <- data.frame(test_pca$x)


# rearrange columns to match metadata order
test_scores <- test_scores[mbx_metadata$raw_filename,]
test_scores$varofinterest <- mbx_metadata$over10

# get convex hulls

pcx <- "PC1"
pcy <- "PC2"

hull <- test_scores %>%
  group_by(varofinterest) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p<-ggplot(test_scores,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=varofinterest))
p<-p+geom_point()
#p<-p+geom_line(mapping = test)
p

p + aes(fill = factor(varofinterest)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none")

#########################
# construct tse objects #
#########################

integrated_metadata <- read_csv(file = "C:/Users/eatmo/Desktop/disease duration paper/integrated_metadata.csv")

shotgun_metadata <- data.frame(integrated_metadata[which(integrated_metadata$data_type == "shotgun_stool"),])

biopsy_metadata <- data.frame(integrated_metadata[which(integrated_metadata$data_type == "biopsy_16s_rectum" | integrated_metadata$data_type == "biopsy_16s_ileum"),])

biopsy_metadata <- column_to_rownames(biopsy_metadata, "external_id")

biopsy_tab <- data.frame(taxonomic_profiles[,rownames(biopsy_metadata)])
rownames(biopsy_tab) <- rownames(taxonomic_profiles)
# remove taxa not found in subset data
biopsy_tax <- data.frame(taxtable2[-which(rowSums(biopsy_tab) == 0),])
biopsy_tab <- biopsy_tab[-which(rowSums(biopsy_tab) == 0),]
colnames(biopsy_tab) <- rownames(biopsy_metadata)
rownames(biopsy_tax) <- rownames(biopsy_tab)


biopsy_tse <- TreeSummarizedExperiment(assays = list(counts = biopsy_tab),
                                colData = biopsy_metadata,
                                rowData = biopsy_tax)

biopsy_tse@assays@data@listData[["counts"]] <- as.matrix(biopsy_tse@assays@data@listData[["counts"]])

biopsy_tse <- subsetByPrevalentTaxa(biopsy_tse, prevalence = 0.1)

biopsy_tse <- transformCounts(biopsy_tse, abund_values = "counts", method = "relabundance")

biopsy_tse <- transformCounts(biopsy_tse, method = "clr", pseudocount = 1)

rownames(speciestax) <- rownames(speciestab)

shotgun_tse <- TreeSummarizedExperiment(assays = list(counts = speciestab),
                                colData = shotgun_metadata,
                                rowData = speciestax)

shotgun_tse@assays@data@listData[["counts"]] <- as.matrix(shotgun_tse@assays@data@listData[["counts"]])

shotgun_tse <- subsetByPrevalentTaxa(shotgun_tse, prevalence = 0.1)

shotgun_tse <- transformCounts(shotgun_tse, method = "clr", pseudocount = 1)

mbx_match_ileum_biopsy_metadata <- mbx_metadata[which(mbx_metadata$Participant_ID %in% biopsy_metadata[which(biopsy_metadata$location == "ileum"),]$Participant_ID),]

ileum_biopsy_metadata <- biopsy_metadata[which(biopsy_metadata$location == "ileum"),]

ileum_16s_mat <- biopsy_tse@assays@data@listData[["clr"]][,rownames(ileum_biopsy_metadata)]

ileum_mbx_mat_unscaled <- rbind(c8pos_processed[,7:29][,mbx_match_ileum_biopsy_metadata$raw_filename], c18neg_processed[,7:29][,mbx_match_ileum_biopsy_metadata$raw_filename], hilicpos_processed[,7:29][,mbx_match_ileum_biopsy_metadata$raw_filename], hilicneg_processed[,7:29][,mbx_match_ileum_biopsy_metadata$raw_filename])

colnames(ileum_16s_mat) <- ileum_biopsy_metadata$Participant_ID
colnames(ileum_mbx_mat_unscaled) <- mbx_match_ileum_biopsy_metadata$Participant_ID

ileum_16s_mbx_unscaled <- rbind(ileum_16s_mat, ileum_mbx_mat_unscaled)

ileum_16s_mbx_pca <- prcomp(t(ileum_16s_mbx_unscaled), scale = FALSE)
ileum_16s_mbx_pca_unitvar <- prcomp(t(ileum_16s_mbx_unscaled), scale = TRUE)
ileum_16s_pca <- prcomp(t(ileum_16s_mat), scale = FALSE)
ileum_mbx_pca <- prcomp(t(ileum_mbx_mat_unscaled), scale = FALSE)


ileum_16s_mbx_scores <- data.frame(ileum_16s_mbx_pca$x)
ileum_16s_mbx_unitvar_scores <- data.frame(ileum_16s_mbx_pca_unitvar$x)

ileum_16s_scores <- data.frame(ileum_16s_pca$x)
ileum_mbx_scores <- data.frame(ileum_mbx_pca$x)


ileum_16s_mbx_scores$varofinterest <- ileum_biopsy_metadata$over10
ileum_16s_mbx_unitvar_scores$varofinterest <- ileum_biopsy_metadata$over10

ileum_16s_scores$varofinterest <- ileum_biopsy_metadata$over10
ileum_mbx_scores$varofinterest <- ileum_biopsy_metadata$over10

# get convex hulls

pcx <- "PC1"
pcy <- "PC2"

hull <- ileum_16s_mbx_scores %>%
  group_by(varofinterest) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p<-ggplot(ileum_16s_mbx_scores,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=varofinterest))
p<-p+geom_point()
#p<-p+geom_line(mapping = test)
p

p + aes(fill = factor(varofinterest)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none")

# add htx stuff

htx_ileum_metadata <- integrated_metadata[which(integrated_metadata$data_type == "htx_ileum"),]

htx_ileum_metadata <- htx_ileum_metadata[which(htx_ileum_metadata$Participant_ID %in% ileum_biopsy_metadata$Participant_ID),]

htx_counts <- read_tsv(file = "C:/Users/eatmo/Desktop/disease duration paper/host_tx_counts.tsv")
# remove duplicate MARCH1 and 2 genes, thanks excel

htx_counts <- htx_counts[-c(21820, 21823),]

htx_counts <- as.data.frame(htx_counts)
rownames(htx_counts) <- htx_counts[,1]

htx_counts <- htx_counts[,htx_ileum_metadata$external_id]

numzeros <- NULL
for (i in 1:nrow(htx_counts)){
  numzeros[i] <- length(which(htx_counts[i,] == 0))    
}
numzeros <- as.data.frame(numzeros)

# filter out genes not present in at least 50% of samples

htx_counts <- htx_counts[-which(numzeros > 8),]

dds <- DESeqDataSetFromMatrix(countData = htx_counts,
                              colData = htx_ileum_metadata,
                              design= ~ over10)

dds <- DESeq(dds)

rld <- rlog(dds, blind=FALSE)

# inspect pca

rld_pca <- prcomp(t(rld@assays@data@listData[[1]]))
rld_pca_scores <- rld_pca$x
rld_pca_scores <- data.frame(rld_pca_scores)
rownames(rld_pca_scores) == htx_ileum_metadata$external_id

rld_pca_scores$varofinterest <- htx_ileum_metadata$disease_duration

pcx <- "PC1"
pcy <- "PC2"

hull <- rld_pca_scores %>%
  group_by(varofinterest) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p<-ggplot(rld_pca_scores,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)), color = varofinterest)) +
  xlab(paste0(colnames(rld_pca_scores)[(as.numeric(strsplit(pcx, "PC")[[1]][2]))], " (", as.numeric(summary(rld_pca)[["importance"]][,1][2]*100), "% variance explained)")) +
  ylab(paste0(colnames(rld_pca_scores)[(as.numeric(strsplit(pcy, "PC")[[1]][2]))], " (", as.numeric(summary(rld_pca)[["importance"]][,2][2]*100), "% variance explained)"))
p<-p+geom_point()
p + aes(fill = factor(varofinterest)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(color="Over 10 years disease duration") 


#

rld_df <- data.frame(rld@assays@data@listData[[1]])

shotgun_df <- data.frame(shotgun_tse@assays@data@listData[["clr"]])

shotgun_ileal_metadata <- shotgun_metadata[which(shotgun_metadata$Participant_ID %in% ileum_biopsy_metadata$Participant_ID),]

shotgun_ileal_metadata <- shotgun_ileal_metadata[-c(13, 18),]

shotgun_df <- shotgun_df[,shotgun_ileal_metadata$external_id]

colnames(shotgun_df) <- shotgun_ileal_metadata$Participant_ID

rld_df <- rld_df[,htx_ileum_metadata$external_id]

colnames(rld_df) <- htx_ileum_metadata$Participant_ID

shotgun_pca <- prcomp(t(shotgun_df))
shotgun_pca_scores <- data.frame(shotgun_pca$x)
shotgun_pca_scores$varofinterest <- shotgun_ileal_metadata$over10

boxplotdf <- plotCounts(dds, gene="UBD", intgroup="over10",
                        returnData=TRUE)

ggplot(boxplotdf,
       aes(x = over10, 
           y = count,
           fill = over10, label = htx_ileum_metadata$Participant_ID)) + 
  geom_boxplot() + geom_point() + geom_text()
  theme(title = element_text(size = 12)) # makes titles smaller

all_dfs <- rbind(ileum_16s_mat, ileum_mbx_mat_unscaled, shotgun_df, rld_df)

all_pca <- prcomp(t(all_dfs), scale = TRUE)
all_scores <- data.frame(all_pca$x)

all_scores$varofinterest <- htx_ileum_metadata$over10

# get convex hulls

pcx <- "PC1"
pcy <- "PC2"

hull <- all_scores %>%
  group_by(varofinterest) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))

p<-ggplot(all_scores,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)), color = varofinterest)) +
  xlab(paste0(colnames(all_scores)[(as.numeric(strsplit(pcx, "PC")[[1]][2]))], " (", as.numeric(summary(all_pca)[["importance"]][,1][2]*100), "% variance explained)")) +
  ylab(paste0(colnames(all_scores)[(as.numeric(strsplit(pcy, "PC")[[1]][2]))], " (", as.numeric(summary(all_pca)[["importance"]][,2][2]*100), "% variance explained)"))
p<-p+geom_point()
p

p + aes(fill = factor(varofinterest)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(color="Over 10 years disease duration") 



ileum_16s_mbx_htx_unscaled <- rbind(ileum_16s_mat, ileum_mbx_mat_unscaled)

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
