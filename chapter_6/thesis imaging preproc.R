# CAC imaging data processing

# step 1: get data for tissue regions in tumours

library(hdf5r)
library(R.matlab)
library(imager)
library(viridis)
library(ggplot2)
#library(mia)
library(MatrixGenerics)
library(gatepoints)
library(irlba)
library(dplyr)

metadata <- read.csv("C:/Users/eatmo/Desktop/CAC/imaging/files/metadata.csv")

# our CAC tumours are:

metadata$Filename[which(metadata$Label_annotation == "Tumour" & metadata$Disease != "sCRC")]

# our sporadic tumours are:

metadata$Filename[which(metadata$Label_annotation == "Tumour" & metadata$Disease == "sCRC")]

#path_to_functions <- "C:/Users/eatmo/Desktop/lauren/imaging_funcs.R"
#source(path_to_functions)

path_to_folder <- "C:/Users/eatmo/Desktop/CAC/imaging/files/"

peaklist <- readMat(paste0(path_to_folder, "cwtpeaks.mat"))

# let's load only our annotated CAC tumours with histological matches (578, 599, 646)

peaklist_578 <- readMat(paste0(path_to_folder, "578Caec.mat"))
peaklist_599 <- readMat(paste0(path_to_folder, "599T.mat"))
peaklist_646 <- readMat(paste0(path_to_folder, "646DescT.mat"))
peaklist_496 <- readMat(paste0(path_to_folder, "496T.mat"))

names(peaklist_578) <- "list.of.peaks"
names(peaklist_599) <- "list.of.peaks"
names(peaklist_646) <- "list.of.peaks"
names(peaklist_496) <- "list.of.peaks"

getMax <- function(peaklist) {
  maxScans <- which(lengths(lapply(peaklist[["list.of.peaks"]], '[[', 1)) == max(lengths(lapply(peaklist[["list.of.peaks"]], '[[', 1))))
}

getMaxmz <- function(peaklist) {
  maxScans <- which(lengths(lapply(peaklist[["list.of.peaks"]], '[[', 1)) == max(lengths(lapply(peaklist[["list.of.peaks"]], '[[', 1))))
  mz_data <- peaklist[["list.of.peaks"]][[maxScans[1]]][[1]][1,]
  return(mz_data)
}

max_578 <- getMax(peaklist_578)
max_599 <- getMax(peaklist_599)
max_646 <- getMax(peaklist_646)
max_496 <- getMax(peaklist_496)

mz_578 <- getMaxmz(peaklist_578)
mz_599 <- getMaxmz(peaklist_599)
mz_646 <- getMaxmz(peaklist_646)
mz_496 <- getMaxmz(peaklist_496)

# fill missing pixels

fillMissing <- function(peaklist) {
  missing <- which(lengths(sapply(peaklist$list.of.peaks, "[[", 1)) < 5)
  notmissing <- which(1:length(peaklist$list.of.peaks) %notin% missing)
  
  # to do: retain pixels with exactly 2 ions by matching to the nearest peak?
  # unlikely to do much, realistically speaking
  
  for(i in 1:length(missing)){
    peaklist$list.of.peaks[[missing[i]]] <- peaklist$list.of.peaks[[as.numeric(notmissing[which.min(abs(missing[i] - notmissing))[1]])]]
  }
  
  return(peaklist)
}

peaklist_578 <- fillMissing(peaklist_578)

dims_578 <- c(200,276)
dims_599 = c(242,182)
dims_646 = c(308,284)
dims_496 <- c(132,124)

img_578 <- filterImg(peaklist_578, numScans = dims_578, maxScans = max_578[1])
img_599 <- filterImg(peaklist_599, numScans = c(242,182), maxScans = max_599)
img_646 <- filterImg(peaklist_646, numScans = c(308,284), maxScans = max_646)
img_496 <- filterImg(peaklist_496, numScans = c(132,124), maxScans = max_496)

# no need to check check if max scan is actually within the tissue region as we
# will recalibrate with a joint m/z axis

# load sporadics

peaklist_211 <- readMat(paste0(path_to_folder, "211T.mat"))
peaklist_550 <- readMat(paste0(path_to_folder, "550T.mat"))
peaklist_339 <- readMat(paste0(path_to_folder, "339T.mat"))
peaklist_263 <- readMat(paste0(path_to_folder, "263T.mat"))

max_211 <- getMax(peaklist_211)
max_550 <- getMax(peaklist_550)
max_339 <- getMax(peaklist_339)
max_263 <- getMax(peaklist_263)

mz_211 <- getMaxmz(peaklist_211)
mz_550 <- getMaxmz(peaklist_550)
mz_339 <- getMaxmz(peaklist_339)
mz_263 <- getMaxmz(peaklist_263)

img_211 <- filterImg(peaklist_211, numScans = c(184,152), maxScans = max_211)
img_550 <- filterImg(peaklist_550, numScans = c(168,172), maxScans = max_550)

peaklist_339 <- fillMissing(peaklist_339)

img_339 <- filterImg(peaklist_339, numScans = c(162,184), maxScans = max_339)
img_263 <- filterImg(peaklist_263, numScans = c(162,146), maxScans = max_263)

dims_211 = c(184,152)
dims_550 = c(168,172)
dims_339 = c(162,184)
dims_263 <- c(162,146)

ticPlot <- function(dims, img) {
  cimgplot <- matrix(data = NA, nrow = dims[1], ncol = dims[2])
  
  for(i in 1:prod(dims)){
    cimgplot[i] <- sum(img[[i]][["intensity"]])
  }
  
  cimgplot <- imager::as.cimg(cimgplot)
  cimgplot <- imager::mirror(cimgplot, "y")
  p <- ggplot(data = as.data.frame(cimgplot)) +
    geom_raster(aes(x = x, y = y, fill = log(value))) + 
    scale_fill_viridis_c(option = "magma", name = "Intensity (A.U)") + coord_fixed()
  p+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
    xlab("Bottom label, if you want one") + ylab(NULL) + ggtitle("Untransformed TIC image") +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

ticPlot(dims_578, img_578)
ticPlot(dims_599, img_599)
ticPlot(dims_646, img_646)
ticPlot(c(182,242), img_646)
ticPlot(dims_646, img_646)

ticPlot(dims_646, img_646)
ticPlot(dims_496, img_496)
ticPlot(dims_211, img_211)
ticPlot(dims_550, img_550)
ticPlot(dims_339, img_339)
ticPlot(dims_263, img_263)

#realign mz scales?

df_211 <- sapply(img_211, '[[', 1)
df_211 <- t(df_211)

mz_211 <- mz_211[-which(colSums(df_211) == 0)]
df_211 <- df_211[,-which(colSums(df_211) == 0)]

cimgplot <- matrix(data = NA, nrow = dims_211[1], ncol = dims_211[2])
for(i in 1:prod(dims_211)){
  cimgplot[i] <- log(((df_211[i,which.min(abs(mz_211 - 885.5499))] / df_211[i,which.min(abs(mz_211 - 255.2330))]) + (df_211[i,which.min(abs(mz_211 - 255.2330))] / df_211[i,which.min(abs(mz_211 - 885.5499))])))
}
cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = log(value+1))) +
  scale_fill_viridis_c(option = "magma") + coord_fixed()
p

tissuepix_211 <- which(is.finite(log(((df_211[,which.min(abs(mz_211 - 885.5499))] / df_211[,which.min(abs(mz_211 - 255.2330))]) + (df_211[,which.min(abs(mz_211 - 255.2330))] / df_211[,which.min(abs(mz_211 - 885.5499))])))) == TRUE)


df_263 <- sapply(img_263, '[[', 1)
df_263 <- t(df_263)

cimgplot <- matrix(data = NA, nrow = dims_263[1], ncol = dims_263[2])
for(i in 1:prod(dims_263)){
  cimgplot[i] <- log(((df_263[i,which.min(abs(mz_263 - 885.5499))] / df_263[i,which.min(abs(mz_263 - 255.2330))]) + (df_263[i,which.min(abs(mz_263 - 255.2330))] / df_263[i,which.min(abs(mz_263 - 885.5499))])))
}
cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = log(value+1))) +
  scale_fill_viridis_c(option = "magma") + coord_fixed()
p

tissuepix_263 <- which(is.finite(log(((df_263[,which.min(abs(mz_263 - 885.5499))] / df_263[,which.min(abs(mz_263 - 255.2330))]) + (df_263[,which.min(abs(mz_263 - 255.2330))] / df_263[,which.min(abs(mz_263 - 885.5499))])))) == TRUE)


df_339 <- sapply(img_339, '[[', 1)
df_339 <- t(df_339)

mz_339 <- mz_339[-which(colSums(df_339) == 0)]
df_339 <- df_339[,-which(colSums(df_339) == 0)]

pca_339 <- prcomp_irlba(log(df_339+1), n = 3)

values <- (as.data.frame(pca_339$x)[,1])
#distances <- dist(values, method = "euclidean")
tree <- genieclust::gclust(values, k = 5)
#tree <- hclust(distances, method = "ward.D2")
cut <- cutree(tree, k = 5)



cimgplot <- matrix(data = NA, nrow = dims_339[1], ncol = dims_339[2])

for(i in 1:prod(dims_339)){
  cimgplot[i] <- cut[i]
}

cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
# tinker with rotations

p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = value)) + 
  scale_fill_viridis_c(option = "magma") + coord_fixed()
p

tissuepix_339 <- which(cut == 5 | cut == 4)


df_496 <- sapply(img_496, '[[', 1)
df_496 <- t(df_496)

pca_496 <- prcomp_irlba(log(df_496+1), n = 5)

values <- (as.data.frame(pca_496$x)[,1:2])
#distances <- dist(values, method = "euclidean")
tree <- genieclust::gclust(values, k = 3)
#tree <- hclust(distances, method = "ward.D2")
cut <- cutree(tree, k = 3)



cimgplot <- matrix(data = NA, nrow = dims_496[1], ncol = dims_496[2])

for(i in 1:prod(dims_496)){
  cimgplot[i] <- cut[i]
}

cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
# tinker with rotations

p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = value)) + 
  scale_fill_viridis_c(option = "magma") + coord_fixed()
p
values[which(cut == 2),] <- 0
tree <- genieclust::gclust(values, k = 5)
#tree <- hclust(distances, method = "ward.D2")
cut <- cutree(tree, k = 5)

cimgplot <- matrix(data = NA, nrow = dims_496[1], ncol = dims_496[2])

for(i in 1:prod(dims_496)){
  cimgplot[i] <- cut[i]
}

cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
# tinker with rotations

p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = value)) + 
  scale_fill_viridis_c(option = "magma") + coord_fixed()
p

values[which(cut == 1),] <- 0
values[which(cut == 2),] <- 0
values[which(cut == 3),] <- 0

tree <- genieclust::gclust(values, k = 5)
#tree <- hclust(distances, method = "ward.D2")
cut <- cutree(tree, k = 5)

cimgplot <- matrix(data = NA, nrow = dims_496[1], ncol = dims_496[2])

for(i in 1:prod(dims_496)){
  cimgplot[i] <- cut[i]
}

cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
# tinker with rotations

p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = value)) + 
  scale_fill_viridis_c(option = "magma") + coord_fixed()
p

values[which(cut == 4),] <- 0
values[which(cut == 5),] <- 0

tree <- genieclust::gclust(values, k = 5)
#tree <- hclust(distances, method = "ward.D2")
cut <- cutree(tree, k = 5)

cimgplot <- matrix(data = NA, nrow = dims_496[1], ncol = dims_496[2])

for(i in 1:prod(dims_496)){
  cimgplot[i] <- cut[i]
}

cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
# tinker with rotations

p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = value)) + 
  scale_fill_viridis_c(option = "magma") + coord_fixed()
p

values[which(cut == 5),] <- 0

tree <- genieclust::gclust(values, k = 5)
#tree <- hclust(distances, method = "ward.D2")
cut <- cutree(tree, k = 5)

cimgplot <- matrix(data = NA, nrow = dims_496[1], ncol = dims_496[2])

for(i in 1:prod(dims_496)){
  cimgplot[i] <- cut[i]
}

cimgplot <- imager::as.cimg(cimgplot)
cimgplot <- imager::mirror(cimgplot, "y")
# tinker with rotations

p <- ggplot(data = as.data.frame(cimgplot)) +
  geom_raster(aes(x = x, y = y, fill = value)) + 
  scale_fill_viridis_c(option = "magma") + coord_fixed()
p

values[which(cut == 5),] <- 0

tissuepix_496 <- which(cut > 1 & cut < 5)

df_550 <- sapply(img_550, '[[', 1)
df_550 <- t(df_550)

tissuepix_550 <- which(is.finite(log(((df_550[,which.min(abs(mz_550 - 885.5499))] / df_550[,which.min(abs(mz_550 - 255.2330))]) + (df_550[,which.min(abs(mz_550 - 255.2330))] / df_550[,which.min(abs(mz_550 - 885.5499))])))) == TRUE)

df_578 <- sapply(img_578, '[[', 1)
df_578 <- t(df_578)

mz_578 <- mz_578[-which(colSums(df_578) == 0)]
df_578 <- df_578[,-which(colSums(df_578) == 0)]

cimgplot <- matrix(data = NA, nrow = dims_578[1], ncol = dims_578[2])
for(i in 1:prod(dims_578)){
  cimgplot[i] <- sum(df_578[i,98])
}
cimgplot <- imager::as.cimg(cimgplot)
px <- as.pixset(cimgplot)
sp <- split_connected(fill(px,3))
tissuepix_578 <- which(sp[[4]] == TRUE)

df_599 <- sapply(img_599, '[[', 1)
df_599 <- t(df_599)

tissuepix_599 <- which(is.finite(log(((df_599[,which.min(abs(mz_599 - 885.5499))] / df_599[,which.min(abs(mz_599 - 255.2330))]) + (df_599[,which.min(abs(mz_599 - 255.2330))] / df_599[,which.min(abs(mz_599 - 885.5499))])))) == TRUE)

df_646 <- sapply(img_646, '[[', 1)
df_646 <- t(df_646)

mz_646 <- mz_646[-which(colSums(df_646) == 0)]
df_646 <- df_646[,-which(colSums(df_646) == 0)]

tissuepix_646 <- which(is.finite(log(((df_646[,which.min(abs(mz_646 - 885.5499))] / df_646[,which.min(abs(mz_646 - 255.2330))]) + (df_646[,which.min(abs(mz_646 - 255.2330))] / df_646[,which.min(abs(mz_646 - 885.5499))])))) == TRUE)

# bullshit clustering for new mz axis

all_mz <- c(mz_211, mz_263, mz_339, mz_496, mz_550, mz_578, mz_599, mz_646)

values <- all_mz
tree <- genieclust::gclust(values, k = 3000)
cut <- cutree(tree, k = 1500)
new_mz <- NULL
for(i in 1:1500){
  new_mz[i] <- median(all_mz[which(cut == i)])
}

img_211 <- filterImg(peaklist_211, numScans = c(184,152), maxScans = max_211, recal = TRUE, new_mz = new_mz)
img_550 <- filterImg(peaklist_550, numScans = c(168,172), maxScans = max_550, recal = TRUE, new_mz = new_mz)
img_339 <- filterImg(peaklist_339, numScans = c(162,184), maxScans = max_339, recal = TRUE, new_mz = new_mz)
img_263 <- filterImg(peaklist_263, numScans = c(162,146), maxScans = max_263, recal = TRUE, new_mz = new_mz)
img_578 <- filterImg(peaklist_578, numScans = dims_578, maxScans = max_578, recal = TRUE, new_mz = new_mz)
img_599 <- filterImg(peaklist_599, numScans = c(242,182), maxScans = max_599, recal = TRUE, new_mz = new_mz)
img_646 <- filterImg(peaklist_646, numScans = c(308,284), maxScans = max_646, recal = TRUE, new_mz = new_mz)
img_496 <- filterImg(peaklist_496, numScans = c(132,124), maxScans = max_496, recal = TRUE, new_mz = new_mz)

df_211 <- sapply(img_211, '[[', 1)
df_211 <- t(df_211)
df_211 <- as.data.frame(df_211)
df_211$Class <- "sCRC"
df_211$Patient <- 211
sub_211 <- df_211[tissuepix_211,]

df_263 <- sapply(img_263, '[[', 1)
df_263 <- t(df_263)
df_263 <- as.data.frame(df_263)
df_263$Class <- "sCRC"
df_263$Patient <- 263
sub_263 <- df_263[tissuepix_263,]

df_339 <- sapply(img_339, '[[', 1)
df_339 <- t(df_339)
df_339 <- as.data.frame(df_339)
df_339$Class <- "sCRC"
df_339$Patient <- 339
sub_339 <- df_339[tissuepix_339,]

df_550 <- sapply(img_550, '[[', 1)
df_550 <- t(df_550)
df_550 <- as.data.frame(df_550)
df_550$Class <- "sCRC"
df_550$Patient <- 550
sub_550 <- df_550[tissuepix_550,]

df_496 <- sapply(img_496, '[[', 1)
df_496 <- t(df_496)
df_496 <- as.data.frame(df_496)
df_496$Class <- "CAC"
df_496$Patient <- 496
sub_496 <- df_496[tissuepix_496,]

df_578 <- sapply(img_578, '[[', 1)
df_578 <- t(df_578)
df_578 <- as.data.frame(df_578)
df_578$Class <- "CAC"
df_578$Patient <- 578
sub_578 <- df_578[tissuepix_578,]

df_599 <- sapply(img_599, '[[', 1)
df_599 <- t(df_599)
df_599 <- as.data.frame(df_599)
df_599$Class <- "CAC"
df_599$Patient <- 599
sub_599 <- df_599[tissuepix_599,]

df_646 <- sapply(img_646, '[[', 1)
df_646 <- t(df_646)
df_646 <- as.data.frame(df_646)
df_646$Class <- "CAC"
df_646$Patient <- 646
sub_646 <- df_646[tissuepix_646,]

all_tissue_df <- rbind(sub_211, sub_263, sub_339, sub_496, sub_550, sub_578, sub_599, sub_646)

all_tissue_df$Side <- "Left"
all_tissue_df$Side[which(all_tissue_df$Patient == 211)] <- "Right"
all_tissue_df$Side[which(all_tissue_df$Patient == 550)] <- "Right"
all_tissue_df$Side[which(all_tissue_df$Patient == 578)] <- "Right"
all_tissue_df$Side[which(all_tissue_df$Patient == 646)] <- "Right"

all_cac <- all_tissue_df[which(all_tissue_df$Class == "CAC"),]
all_crc <- all_tissue_df[which(all_tissue_df$Class == "sCRC"),]

#load("C:/Users/eatmo/Desktop/CAC/mergeddata.RData")

which(rowSums(all_tissue_df[,1:1500]) == 0)

all_tissue_df <- all_tissue_df[-which(rowSums(all_tissue_df[,1:1500]) == 0),]

all_tissue_df_uv <- all_tissue_df[,1:1500]

all_tissue_df_uv <- t(as.data.frame(apply(all_tissue_df_uv, 1, function(row) row / sd(row))))

all_tissue_df_tic <- all_tissue_df[,1:1500]

all_tissue_df_tic <- apply(all_tissue_df_tic, 1, function(row) row / sum(row))
all_tissue_df_tic <- t(all_tissue_df_tic)

all_tissue_df_clr <- calc_clr(t(all_tissue_df_tic), pseudocount = 1)
all_tissue_df_rclr <- calc_rclr(t(all_tissue_df_tic))

# which col has the most var

mz_vars <- NULL

for(i in 1:1500){
  mz_vars[i] <- var(all_tissue_df[,i])
}

# obj_in <- all_tissue_df[,1:1500]
# obj_in <- log(all_tissue_df_uv[,1:1500]+1)
obj_in <- all_tissue_df_uv
set.seed(12345)
pca_obj <- prcomp_irlba(log(obj_in+1), n = 3)
df_out <- as.data.frame(pca_obj$x)
df_out$Class <- all_tissue_df$Class
pcx <- "PC1"
pcy <- "PC2"
p<-ggplot(df_out,aes(x=eval(parse(text=pcx)),y=eval(parse(text=pcy)),color=Class))
p<-p+geom_point()
#p
# get convex hulls
hull <- df_out %>%
  group_by(Class) %>%
  dplyr::slice(chull(eval(parse(text=pcx)), eval(parse(text=pcy))))
p + aes(fill = factor(Class)) + geom_polygon(data = hull, alpha = 0.1) + guides(fill = "none") +
  labs(x = paste(pcx, as.numeric(summary(pca_obj)[["importance"]][,as.numeric(strsplit(pcx, "PC")[[1]][2])][2])*100, "% variance explained"),
       y = paste(pcy, as.numeric(summary(pca_obj)[["importance"]][,as.numeric(strsplit(pcy, "PC")[[1]][2])][2])*100, "% variance explained"))
loadings_df <- as.data.frame(pca_obj$rotation)
loadings_df$mzs <- new_mz
ggplot(loadings_df, aes(x = mzs, y = PC1)) +
  labs(x = "m/z", y = "PC1 Score") +
  geom_segment(aes(x = mzs, xend = mzs, yend = 0, y = PC1), linewidth = 0.1)
ggplot(loadings_df, aes(x = mzs, y = PC2)) +
  labs(x = "m/z", y = "PC2 Score") +
  geom_segment(aes(x = mzs, xend = mzs, yend = 0, y = PC2), linewidth = 0.1)
ggplot(loadings_df, aes(x = mzs, y = PC3)) +
  labs(x = "m/z", y = "PC3 Score") +
  geom_segment(aes(x = mzs, xend = mzs, yend = 0, y = PC3), linewidth = 0.1)
ggplot(loadings_df, aes(x = PC1, y = PC2)) +
  labs(x = "PC1 Score", y = "PC2 Score") +
  geom_point()

avg_211 <- colMeans(sub_211)
avg_263 <- colMeans(sub_263)
avg_339 <- colMeans(sub_339)
avg_496 <- colMeans(sub_496)
avg_550 <- colMeans(sub_550)
avg_578 <- colMeans(sub_578)
avg_599 <- colMeans(sub_599)
avg_646 <- colMeans(sub_646)









roi_df_211 <- df_211[roi_211,]
roi_df_263 <- df_263[roi_263,]
roi_df_339 <- df_339[roi_339,]
roi_df_496 <- df_496[roi_496,]
roi_df_550 <- df_550[roi_550,]
roi_df_578 <- df_578[roi_578,]
roi_df_599 <- df_599[roi_599,]
roi_df_646 <- df_646[roi_646,]

roi_df_all <- rbind(roi_df_211,roi_df_263,roi_df_339,roi_df_496,roi_df_550,roi_df_578,roi_df_599,roi_df_646)

which(rowSums(roi_df_all[,1:1500]) == 0)

roi_df_all <- roi_df_all[-which(rowSums(roi_df_all[,1:1500]) == 0),]

roi_df_uv <- roi_df_all[,1:1497]

roi_df_uv <- t(as.data.frame(apply(roi_df_uv, 1, function(row) row / sd(row))))

roi_df_tic <- roi_df_all[,1:1497]

roi_df_tic <- apply(roi_df_tic, 1, function(row) row / sum(row))
roi_df_tic <- t(roi_df_tic)

roi_df_clr <- calc_clr(t(roi_df_tic), pseudocount = 1)
roi_df_rclr <- calc_rclr(t(roi_df_tic))

# mean spectra

roi_df_all <- roi_df_all[,keepcol]
roi_df_clr <- roi_df_clr[keepcol,]
roi_df_rclr <- roi_df_rclr[keepcol,]
roi_df_tic <- roi_df_tic[,keepcol]
roi_df_uv <- roi_df_uv[,keepcol]
