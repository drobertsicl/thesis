uniqueIds <- c("CSM5MCTZ",  "CSM5MCUM",  "CSM5MCX3",  "CSM67UFZ",  "CSM79HMN",  "CSM79HQR",  "CSM7KORU",  "CSM9X1XW",  "CSM9X21J",  "CSM9X21P",  "MSM79HD6",  "MSM79H7M",  "MSM79H7G",  "MSMA26AZ",  "MSMB4LZ4",  "MSM9VZNX",  "MSM9VZIM",  "MSM9VZIQ",  "MSMAPC5D",  "MSMB4LZR",  "PSMA265X",  "PSMA2675",  "PSMA269W",  "PSMA26A3",  "PSMB4MC5",  "CSM5FZ48",  "CSM5FZ4K",  "CSM5MCYK",  "CSM67UBZ",  "CSM79HLM",  "CSM5FZ4A",  "CSM5MCUC",  "CSM5MCXL",  "CSM67UDN",  "CSM79HLA",  "CSM5FZ4O",  "CSM5MCUU",  "CSM5MCY8",  "CSM67UEA",  "CSM67UEI",  "CSM5MCXB",  "CSM5MCYO",  "CSM67UEW",  "CSM67UF5",  "CSM7KOMT",  "CSM5MCYU",  "CSM67U9P",  "CSM79HGH",  "CSM79HPS",  "CSM7KONS",  "CSM67U9T",  "CSM67UAQ",  "CSM79HIJ",  "CSM7KOPI",  "CSM7KOPQ",  "CSM79HJG",  "CSM7KOLY",  "CSM9X1Y5",  "CSM79HQL",  "CSM7KOP8",  "CSM7KOPC",  "CSM7KOPE",  "CSM7KOU9",  "CSM7KOUF",  "CSM9X22G",  "CSM9X22Q",  "CSM7KOQ7",  "CSM7KOUL",  "CSM7KOUR",  "CSM9X213",  "CSM9X211",  "CSM9X22U",  "CSM9X23L",  "CSM9X23P")

# PATIENT ID M2008 IS A CROHN'S CASE. DO NOT USE

metadata <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/ibdmdb metadata.csv", col_names = TRUE)
metadatauc <- metadata[which(metadata$diagnosis == "UC"),]

metadatatest <- metadatauc[which(metadatauc$`External ID` %in% uniqueIds),]

which(colnames(metadata) == "Location")
metadatabxlocs <- metadata[,299:311]

metadatabxs <- metadatauc[,433:456]
metadatabxs <- cbind

sampleids <- c("CSM5FZ48","CSM5FZ4K","CSM5MCYK","CSM67UBZ","CSM79HLM","CSM5FZ4A","CSM5MCUC","CSM5MCXL","CSM67UDN","CSM79HLA","CSM5FZ4O","CSM5MCUU","CSM5MCY8","CSM67UEA","CSM67UEI","CSM5MCTZ","CSM5MCUM","CSM5MCX3","CSM67UFZ","CSM79HMN","CSM5MCXB","CSM5MCYO","CSM67UEW","CSM67UF5","CSM7KOMT","CSM5MCYU","CSM67U9P","CSM79HGH","CSM79HPS","CSM7KONS","CSM67U9T","CSM67UAQ","CSM79HIJ","CSM7KOPI","CSM7KOPQ","CSM79HJG","CSM7KOLY","CSM9X1Y5","CSM79HQL","CSM7KOP8","CSM7KOPC","CSM7KOPE","CSM7KOU9","CSM7KOUF","CSM9X22G","CSM9X22Q","CSM79HQR","CSM7KORU","CSM9X1XW","CSM9X21J","CSM9X21P","CSM7KOQ7","CSM7KOUL","CSM7KOUR","CSM9X213","CSM9X211","CSM9X22U","CSM9X23L","CSM9X23P","ESM5MEEG","ESM5GEZ6","ESM5MEBU","ESM5MEC5","ESM718SY","HSM5MD87","HSM5MD43","HSM5MD3Y","HSM5MD41","HSM6XRTM","HSM6XRTQ","HSM67VH1","HSM6XRQE","HSM6XRRD","HSM6XRV2","HSM7J4PE","HSM7J4M4","HSM7J4ME","HSM7CZ1T","HSM7J4GR","HSM7J4NE","HSMA33O1","HSMA33O5","HSM7J4HO","HSM7J4HS","HSM7J4JZ","HSM7J4K8","HSMA33OR","HSMA33OZ","HSMA33MI","HSM7J4JT","HSM7J4NU","HSMA33OJ","HSMA33M8","HSM7J4JV","HSMA33MZ","HSMA33IY","HSMA33RL","HSMA33NW","HSMA33RX","MSM5LLFK","MSM5LLFO","MSM6J2HF","MSM6J2OH","MSM6J2OP","MSM79HAJ","MSM79HAN","MSM9VZLP","MSM9VZLV","MSMA26ET","MSM79HDA","MSM79HDE","MSM9VZEW","MSM9VZOU","MSM9VZP3","MSMA26EJ","PSM7J12V","PSM7J134")

metabolomicsdata <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/iHMP_metabolomics.csv", col_names = TRUE)
selector <- c(colnames(metabolomicsdata[,1:7]), sampleids)

#library(data.table)

metsubset <- subset(metabolomicsdata, select = selector)
metsubset <- setcolorder(metsubset, as.character(selector))

write_csv(metsubset, file = "D:/Users/Duncan Roberts/Desktop/disease duration paper/mbx_established_plusnew_nonaextent.csv")

###### go time

mSet <- InitDataObjects("pktable", "utils", FALSE)
mSet <- Read.TextData(mSet, "D:/Users/Duncan Roberts/Desktop/disease duration paper/mbx_established_plusnew_nonaextent_formetabo.csv", "col", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-FilterVariable(mSet, "iqr", "F", 25)
mSet<-PreparePrenormData(mSet)
mSet <- Normalization(mSet, "MedianNorm", "LogNorm", "ParetoNorm")

mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet <- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

normalized_set <- mSet[["dataSet"]][["norm"]]
ordered_normalized_set <- normalized_set[order(row.names(normalized_set)), ]

metadata <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/metabolomics_metadata_establisheduc_plusnew_nonaextent.csv", col_names = TRUE)

newmetid <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/C18n_Metabolites_ID_AfterPublication_2021-05-17.csv", col_names = TRUE)

for(i in 1:nrow(newmetid)){
  newmetid[i,2] <- paste0(newmetid[i,1], "_", newmetid[i,2])
}

for(i in 1:nrow(newmetid)){
  metabolomicsdata[which(metabolomicsdata[,7] == as.character(newmetid[i,2])),6] <- newmetid[i,5]
}

write_csv(metabolomicsdata, file = "iHMP2 new c18 annotations.csv")

onlyident <- metabolomicsdata[which(nchar(metabolomicsdata[,6]) > 0),]

onlyidentuc <- subset(onlyident, select = selector)

write_csv(onlyidentuc, file = "ucpts only identified metabolites.csv")

mSet<-InitDataObjects("pktable", "ts", FALSE)
mSet<-SetDesignType(mSet, "multi")
mSet<-Read.TextDataTs(mSet, "D:/Users/Duncan Roberts/Desktop/disease duration paper/data_normalized_quintile_log_pareto.csv", "colu");
mSet<-ReadMetaData(mSet, "D:/Users/Duncan Roberts/Desktop/disease duration paper/simplemetadata.csv");
mSet<-SanityCheckData(mSet)
mSet<-SanityCheckData(mSet)
mSet<-RemoveMissingPercent(mSet, percent=0.5)
mSet<-ImputeMissingVar(mSet, method="min")
mSet<-SanityCheckMeta(mSet, 1)
mSet<-SanityCheckData(mSet)

corrs <- read_csv("D:/Users/Duncan Roberts/Downloads/correlation_feature.csv", col_names = TRUE)
corrs <- as.data.frame(corrs)
corrs <- corrs[which(corrs[,5] < 0.05),]

metabolomicsdata <- as.data.frame(metabolomicsdata)

corrs[,6] <- "NA"

for(i in 1:nrow(corrs)){
  corrs[i,6] <- metabolomicsdata[which(metabolomicsdata[,7] == as.character(corrs[i,1])),6]
}

splot <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/oplsda_splot_data.csv", col_names = TRUE)
splotnamed <- as.data.frame(matrix(data = NA, nrow = 267, ncol = 4))
j <- 0
for(i in 1:nrow(splot)){
  if(!is.na(metabolomicsdata[which(metabolomicsdata[,7] == as.character(splot[i,1])),6])){
    #splot[i,1] <- metabolomicsdata[which(metabolomicsdata[,7] == as.character(splot[i,1])),6]
    j <- j + 1
    splotnamed[j,] <- splot[i,]
    splotnamed[j,4] <- metabolomicsdata[which(metabolomicsdata[,7] == as.character(splot[i,1])),6]
  }
}

metsubset <- subset(metabolomicsdata, select = selector)

normeddata <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/data_normalized_quintile_log_pareto.csv", col_names = TRUE)

duration <- metadata[,1]
yoink <- which(metabolomicsdata[,6] == as.character(splotnamed[18,1]))
metab <- normeddata[which(normeddata[,1] == as.character(metabolomicsdata[yoink,7])),2:118]
test <- cbind(duration, t(metab))
test[,2] <- as.numeric(as.character(test[,2]))
colnames(test)[2] <- "metabolite"

which(duration != 0)
test2 <- test[which(duration != 0),]

ggplot(test2, aes(x=Duration, y=metabolite)) + 
  geom_point()+
  geom_smooth(method=lm)

establisheduc <- rownames(test[which(test[,1] != 0),])
establisheduc <- c(colnames(metabolomicsdata[,1:7]), establisheduc)

metsubset2 <- subset(metabolomicsdata, select = establisheduc)

write_csv(metsubset2, file = "D:/Users/Duncan Roberts/Desktop/disease duration paper/mbx_established_nonaextent.csv")
test3 <- rbind(establisheduc, metsubset2)
over10 <- metadata[which(test[,1] != 0),2]
over10 <- as.data.frame(over10)
over10 <- t(over10)
test3[1,8:70] <- over10
write_csv(test3, file = "D:/Users/Duncan Roberts/Desktop/disease duration paper/mbx_established_nonaextent.csv")

splot2 <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/splot_data_nonewdiags.csv", col_names = TRUE)
splotnamed2 <- as.data.frame(matrix(data = NA, nrow = 267, ncol = 3))
j <- 0
for(i in 1:nrow(splot2)){
  if(!is.na(metabolomicsdata[which(metabolomicsdata[,7] == as.character(splot2[i,1])),6])){
    splot2[i,1] <- metabolomicsdata[which(metabolomicsdata[,7] == as.character(splot2[i,1])),6]
    j <- j + 1
    splotnamed2[j,] <- splot2[i,]
  }
}

duration <- metadata[,1]
yoink <- which(metabolomicsdata[,6] == as.character(splotnamed2[1,1]))
metab <- normeddata[which(normeddata[,1] == as.character(metabolomicsdata[yoink,7])),2:118]
test <- cbind(duration, t(metab))
test[,2] <- as.numeric(as.character(test[,2]))
colnames(test)[2] <- "metabolite"

which(duration != 0)
test2 <- test[which(duration != 0),]

ggplot(test2, aes(x=Duration, y=metabolite)) + 
  geom_point()+
  geom_smooth(method=lm)

test2[which(test2[,1] < 10),1] <- "Under 10"
test2[which(test2[,1] != "Under 10"),1] <- "Over 10"

vecselect <- as.character(unlist(metadatauc[,9]))

which(vecselect %in% sampleids)

metadatauc <- metadata[which(metadata$`External ID` %in% sampleids),]

metadatauc <- metadatauc[which(metadatauc$data_type == "metabolomics")]

#???????? why are cols wrong
dairy <- metadata[which(metadata[,9] %in% sampleids), c(2,77:78)]
dairy <- metadata[, c(2,84:85)]
dairy <- dairy[!duplicated(dairy),]

dairy <- metadatauc[which(vecselect %in% sampleids), c(1:2,10,84:85)]
#dairy <- dairy[!duplicated(dairy$`Participant ID`),]

coffee <- metadatauc[which(vecselect %in% sampleids), c(1:2,10,108)]

#dairy <- cbind(dairy, duration)

dairy[which(dairy[,1] == "YES"),1] <- "Over 10"
dairy[which(dairy[,1] != "Over 10"),1] <- "Under 10"

dairy <- as.data.frame(dairy)

dairyover10s <- dairy[which(dairy[,1] == "Over 10"),]
dairyunder10s <- dairy[which(dairy[,1] != "Over 10"),]
dairytableover10s_yog <- table(dairyover10s[,2])
dairytableunder10s_yog <- table(dairyunder10s[,2])
dairytableover10s_milk <- table(dairyover10s[,3])
dairytableunder10s_milk <- table(dairyunder10s[,3])

dairyyogbars <- as.data.frame(matrix(data = NA, nrow = 5, ncol = 2))
colnames(dairyyogbars) <- c("Under 10", "Over 10")
rownames(dairyyogbars) <- names(dairytableunder10s_yog)

dairymilkbars <- as.data.frame(matrix(data = NA, nrow = 5, ncol = 2))
colnames(dairymilkbars) <- c("Under 10", "Over 10")
rownames(dairymilkbars) <- names(dairytableunder10s_milk)

for(i in 1:5){
  dairyyogbars[i,1] <- (dairytableunder10s_yog[i]/sum(dairytableunder10s_yog))*100
  dairyyogbars[i,2] <- (dairytableover10s_yog[i]/sum(dairytableover10s_yog))*100
}
dairyyogbars[5,2] <- 0

for(i in 1:5){
  dairymilkbars[i,1] <- (dairytableunder10s_milk[i]/sum(dairytableunder10s_milk))*100
  dairymilkbars[i,2] <- (dairytableover10s_milk[i]/sum(dairytableover10s_milk))*100
}

dairy_chi <- dairymilkbars

for(i in 1:5){
  dairy_chi[i,1] <- (dairytableunder10s_milk[i])
  dairy_chi[i,2] <- (dairytableover10s_milk[i])
}

dairy_chi <- as.matrix(dairy_chi)
rownames(dairy_chi) <- rownames(dairymilkbars)
dairy_chi <- as.table(dairy_chi)

dairyyogbars2 <- c(dairyyogbars[,1], dairyyogbars[,2])
dairyyogbars2 <- cbind(dairyyogbars2, dairyyogbars2, dairyyogbars2)
colnames(dairyyogbars2) <- c("percent", "answer", "group")
dairyyogbars2[,2] <- c(rownames(dairyyogbars, rownames(dairyyogbars)))
dairyyogbars2[1:5,3] <- "Under 10"
dairyyogbars2[6:10,3] <- "Over 10"

dairyyogbars2 <- as.data.frame(dairyyogbars2)
dairyyogbars2[,1] <- as.numeric(as.character(dairyyogbars2[,1]))

dairymilkbars2 <- c(dairymilkbars[,1], dairymilkbars[,2])
dairymilkbars2 <- cbind(dairymilkbars2, dairymilkbars2, dairymilkbars2)
colnames(dairymilkbars2) <- c("percent", "answer", "group")
dairymilkbars2[,2] <- c(rownames(dairymilkbars, rownames(dairymilkbars)))
dairymilkbars2[1:5,3] <- "Under 10"
dairymilkbars2[6:10,3] <- "Over 10"

dairymilkbars2 <- as.data.frame(dairymilkbars2)
dairymilkbars2[,1] <- as.numeric(as.character(dairymilkbars2[,1]))

dairyyogbars3 <- dairyyogbars2
dairyyogbars3$group <- with(dairyyogbars2, reorder(group, percent))

ggplot(dairyyogbars3, aes(fill=answer, y=percent, x=group)) + 
  geom_bar(position="stack", stat="identity") + xlab("Disease duration") + ylab("Proportion") + labs(fill = NULL, title = "(B) Yogurts or fermented dairy") 

dairymilkbars3 <- dairymilkbars2
dairymilkbars3$group <- with(dairymilkbars2, reorder(group, percent))

ggplot(dairymilkbars3, aes(fill=answer, y=percent, x=group)) + 
  geom_bar(position="stack", stat="identity") + xlab("Disease duration") + ylab("Proportion") + labs(fill = NULL, title = "(A) Milk or milk products") 

lct <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/lactase_counts.csv", col_names = TRUE)
fxrpxr <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/fxrpxr.csv", col_names = TRUE)

ptids <- as.vector(unlist(distinct(metadata[,10])))

fullmeta <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/hmp2_metadata.csv", col_names = TRUE)
transcripts <- subset(fullmeta, biopsy_location == "Ileum" & data_type == "host_transcriptomics")
transcriptsrectum <- subset(fullmeta, biopsy_location == "Rectum" & data_type == "host_transcriptomics")
alltranscripts <- subset(fullmeta, data_type == "host_transcriptomics")

samples <- as.character(unlist(transcripts[,3]))
samplesrectum <- as.character(unlist(transcriptsrectum[,3]))
transcripts <- transcripts[which(samples %in% ptids),]
transcriptsrectum <- transcriptsrectum[which(samples %in% ptids),]

lctselect <- as.character(unlist(transcripts[,2]))
rectselect <- as.character(unlist(transcriptsrectum[,2]))


lct <- subset(lct, select = lctselect)
fxr <- subset(fxrpxr[1,], select = lctselect)
pxr <- subset(fxrpxr[2,], select = lctselect)

fxrrectum <- subset(fxrpxr[1,], select = rectselect)
pxrrectum <- subset(fxrpxr[2,], select = rectselect)

fxrptids <- as.character(unlist(transcriptsrectum[,3]))

# remove additional biopsy
lct <- lct[,-15]
fxr <- fxr[,-15]
pxr <- pxr[,-15]
transcripts <- transcripts[-15,]
# 4042 and 2064 missing
lctptids <- ptids[-21]
lctptids <- lctptids[-18]

tempsampleids <- as.character(unlist(metadata[,10]))
tempsampleids %in% lctptids
tempsampleids %in% fxrptids


over10short <- metadata[which(tempsampleids %in% lctptids),c(2,10)]
over10short <- over10short[!duplicated(over10short),]
over10short[,2] <- as.numeric(lct[1,])
rownames(over10short) <- lctptids
colnames(over10short)[2] <- "Lactase gene count"
subsetmeta <- metadatauc[which(metadata)]

over10short$Over10 <- with(over10short, reorder(Over10, `Lactase gene count`)
                           
over10short %>%
 ggplot(aes(x=Over10, y=`Lactase gene count`, fill=Over10)) +
 geom_boxplot() + stat_boxplot(geom = "errorbar", width = 0.15) + guides(fill = "none") +
 labs(x="Years of disease activity", y = "LCT gene transcription (copy number)")

over10fxr <- over10short
over10fxr[,2] <- as.numeric(fxr[1,])
colnames(over10fxr)[2] <- "FXR gene count"

over10fxr %>%
 ggplot( aes(x=Over10, y=`FXR gene count`, fill=Over10)) +
 geom_boxplot()

over10pxr <- over10short
over10pxr[,2] <- as.numeric(pxr[1,])
colnames(over10pxr)[2] <- "PXR gene count"

over10pxr %>%
 ggplot( aes(x=Over10, y=`PXR gene count`, fill=Over10)) +
 geom_boxplot()

over10fxrrectum <- metadata[which(tempsampleids %in% fxrptids),c(2,10)]
over10fxrrectum <- over10fxrrectum[!duplicated(over10fxrrectum),]
rownames(over10fxrrectum) <- c("C3003", "C3004", "C3005", "C3006", "C3011", "C3013", "C3015", "C3029", "C3032", "C3034", "C3037", "E5004", "H4040")
fxrrectselect <- as.character(unlist(transcriptsrectum[which(as.character(unlist(transcriptsrectum[,3])) %in% c("C3003", "C3004", "C3005", "C3006", "C3011", "C3013", "C3015", "C3029", "C3032", "C3034", "C3037", "E5004", "H4040")),2]))

fxrrectum <- subset(fxrpxr[1,], select = fxrrectselect)
over10fxrrectum[,2] <- as.numeric(fxrrectum)
colnames(over10fxrrectum)[2] <- "FXR gene count"

over10fxrrectum %>%
 ggplot( aes(x=Over10, y=`FXR gene count`, fill=Over10)) +
 geom_boxplot()

allbxfxr <- alltranscripts[which(as.character(unlist(alltranscripts[,3])) %in% ptids),2]
allbxfxrlocs <- alltranscripts[which(as.character(unlist(alltranscripts[,3])) %in% ptids),41]
allbxfxrptids <- alltranscripts[which(as.character(unlist(alltranscripts[,3])) %in% ptids),3]
allbxfxr <- as.character(unlist(allbxfxr))
allbxfxrover10 <- allbxfxrptids
for(i in 1:nrow(allbxfxrptids)){
 allbxfxrover10[i,1] <- ptover10[which(ptids == as.character(allbxfxrptids[i,1]))]
}

allfxr <- subset(fxrpxr[1,], select = allbxfxr)
allfxr2 <- cbind(allbxfxrover10, t(allfxr))
colnames(allfxr2) <- c("Over10", "FXR gene count")

allfxr2 %>%
 ggplot( aes(x=Over10, y=`FXR gene count`, fill=Over10)) +
 geom_boxplot()

ptover10 <- as.character(unlist(metadata[which(!duplicated(metadata[,10])),2]))

rectal16s <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/rectal16s.csv", col_names = TRUE)

rowstodelete <- NULL
for(i in 1:nrow(rectal16s)){
 if(isTRUE(rowSums(rectal16s[i,7:26]) == 0) == TRUE){
   rowstodelete <- c(rowstodelete, i)
 }
 else{
   
 }
}

rectal16s <- rectal16s[-rowstodelete,]

which(tempsampleids %in% colnames(rectal16s[7:26]))

tocorrelate <- tempsampleids[which(tempsampleids %in% colnames(rectal16s[7:26]))]

mbx117 <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/data_processed_117.csv", col_names = TRUE)

tempsubset <- mbx117[-1,2:118]
tempsubset <- tempsubset[, which(tempsampleids %in% colnames(rectal16s[7:26]))]

tempsubset2 <- matrix(data = as.matrix(tempsubset), nrow = 2499, ncol = 102)
tempsubset2 <- matrix(as.numeric(tempsubset2), nrow = 2499, ncol = 102)
rownames(tempsubset2) <- as.character(unlist(mbx117[-1,1]))
colnames(tempsubset2) <- colnames(tempsubset)

cor.test(x = tempsubset2[26,], )

corrs <- read_csv("D:/Users/Duncan Roberts/Desktop/disease duration paper/correlation_feature(1).csv", col_names = TRUE)
corrsnamed <- as.data.frame(matrix(data = NA, nrow = 2500, ncol = 6))
j <- 0
for(i in 1:nrow(corrs)){
 if(!is.na(metabolomicsdata[which(metabolomicsdata[,7] == as.character(corrs[i,1])),6])){
   #splot[i,1] <- metabolomicsdata[which(metabolomicsdata[,7] == as.character(splot[i,1])),6]
   j <- j + 1
   corrsnamed[j,] <- corrs[i,]
   corrsnamed[j,6] <- metabolomicsdata[which(metabolomicsdata[,7] == as.character(corrs[i,1])),6]
 }
}

colnames(corrsnamed) <- colnames(corrs)
colnames(corrsnamed)[6] <- "metabolite"

corrsnamed <- corrsnamed[which(corrsnamed$FDR < 0.05),]

write.csv(corrsnamed, file = "D:/Users/Duncan Roberts/Desktop/disease duration paper/correlationsid.csv")

test <- corrsnamed[!duplicated(corrsnamed[,6]),]

write.csv(test, file = "D:/Users/Duncan Roberts/Desktop/disease duration paper/correlationsidnodupes.csv")

simplemetadata <- read_csv(file = "D:/Users/Duncan Roberts/Desktop/disease duration paper/simplemetadata.csv")

##### diet

procmeat <- diet[,c(1, 17)]
redmeat <- diet[,c(1, 18)]
whitemeat<- diet[,c(1, 19)]
shellfish<- diet[,c(1, 20)]
fish<- diet[,c(1, 21)]
wholegrain<- diet[,c(1, 14)]
beans<- diet[,c(1, 13)]
starch<- diet[,c(1, 15)]
eggs<- diet[,c(1, 16)]

# ctrl h replace with var above

whitemeat[which(whitemeat[,1] == "YES"),1] <- "Over 10"
whitemeat[which(whitemeat[,1] != "Over 10"),1] <- "Under 10"

whitemeatover10s <- whitemeat[which(whitemeat[,1] == "Over 10"),]
whitemeatunder10s <- whitemeat[which(whitemeat[,1] != "Over 10"),]

whitemeattableover10s <- table(whitemeatover10s[,2])
whitemeattableunder10s <- table(whitemeatunder10s[,2])

whitemeatbars <- as.data.frame(matrix(data = NA, nrow = 5, ncol = 2))
colnames(whitemeatbars) <- c("Under 10", "Over 10")
rownames(whitemeatbars) <- names(whitemeattableover10s)

# total num for chi sq
for(i in 1:5){
 whitemeatbars[i,1] <- (whitemeattableunder10s[i])
 whitemeatbars[i,2] <- (whitemeattableover10s[i])
}


whitemeatbars[5,2] <- 0

# percent for plotting
for(i in 1:5){
 whitemeatbars[i,1] <- (whitemeattableunder10s[i]/sum(whitemeattableunder10s))*100
 whitemeatbars[i,2] <- (whitemeattableover10s[i]/sum(whitemeattableover10s))*100
}

# plug na holes with 0

whitemeatbars[5,1] <- 0

whitemeatbars2 <- c(whitemeatbars[,1], whitemeatbars[,2])
whitemeatbars2 <- cbind(whitemeatbars2, whitemeatbars2, whitemeatbars2)
colnames(whitemeatbars2) <- c("percent", "answer", "group")
whitemeatbars2[,2] <- c(rownames(whitemeatbars, rownames(whitemeatbars)))
whitemeatbars2[1:5,3] <- "Under 10"
whitemeatbars2[6:10,3] <- "Over 10"

whitemeatbars2 <- as.data.frame(whitemeatbars2)

whitemeatbars2[,1] <- as.numeric(as.character(whitemeatbars2[,1]))

whitemeatbars2$group <- with(whitemeatbars2, reorder(group, percent))

whitemeatbars2[,1] <- as.numeric(as.character(whitemeatbars2[,1]))

ggplot(whitemeatbars2, aes(fill=answer, y=percent, x=group)) + 
 geom_bar(position="stack", stat="identity") + xlab("Disease duration") + ylab("Proportion") + labs(fill = NULL, title = "Red meat intake") 



whitemeat2[which(whitemeat2$`Tea or whitemeat no sugar and no sugar replacement` == "No, I did not consume these products in the last 7 days"),4] <- "None"













for(i in 1:nrow(rectal16s))
 
 library(lme4)
