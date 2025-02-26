library("tidyverse")
library("tidyr")

refseqlist <- readRDS("C:/16s/blasttest/refseqlist.rds")

for(i in 1:length(refseqlist)){
  names(refseqlist[[i]])[1] <- "Gene"
  names(refseqlist[[i]])[2] <- "Name"
  names(refseqlist[[i]])[3] <- "Definition"
  names(refseqlist[[i]])[4] <- "copy number"
}

sckos <- c("K00057","K00059","K00075","K00088","K00133","K00134","K00384","K00525","K00554","K00566","K00600","K00604","K00609","K00615","K00655","K00762","K00789","K00790","K00791","K00820","K00858","K00859","K00927","K00939","K00942","K00948","K00954","K00962","K00981","K01000","K01056","K01265","K01358","K01409","K01462","K01465","K01491","K01591","K01689","K01714","K01756","K01775","K01783","K01803","K01810","K01866","K01867","K01868","K01869","K01870","K01872","K01873","K01874","K01875","K01876","K01881","K01883","K01887","K01889","K01890","K01892","K01915","K01921","K01924","K01925","K01929","K01937","K01939","K01945","K01951","K01955","K01956","K01972","K01990","K01992","K02003","K02004","K02078",
           "K02108","K02109","K02110","K02111","K02112","K02114","K02115","K02313","K02314","K02316","K02335","K02337","K02338","K02340","K02341","K02343","K02355","K02356","K02357","K02358","K02469","K02470","K02493","K02495","K02518","K02519","K02520","K02528","K02563","K02600","K02601","K02834","K02835","K02836","K02838","K02860","K02863","K02864","K02867","K02871","K02874","K02876","K02878","K02879","K02881","K02884","K02886","K02887","K02888","K02890","K02892","K02895","K02899","K02902","K02904","K02906","K02909","K02911","K02913","K02916","K02926","K02931","K02933","K02935","K02939","K02945","K02946","K02948","K02950","K02952","K02954","K02956","K02959","K02961","K02963","K02965","K02967","K02968","K02982","K02986","K02988","K02990","K02992","K02994","K02996","K03040","K03043","K03046","K03070","K03075","K03076","K03086","K03100","K03101","K03106","K03110","K03111","K03168","K03177","K03217","K03218","K03424","K03438","K03466","K03470","K03501","K03530","K03531","K03544","K03545","K03550","K03551","K03553","K03588","K03595","K03596","K03624","K03625","K03631","K03655","K03657","K03664","K03671","K03685","K03686","K03687","K03701","K03702","K03703","K03723","K03798","K03977",
           "K03979","K04043","K04066","K04075","K04077","K04078","K04096","K04485","K04487","K06173","K06178","K06180","K06187","K06207","K06942","K07056","K07447","K07478","K09710","K09761","K09903","K10773","K11753","K11754")

medians <- NULL
for(i in 1:length(refseqlist)){
  medians[i] <- median(unlist(refseqlist[[1]][which(unlist(refseqlist[[i]][,1]) %in% sckos),4]))
}

# replace nas with 0?

medians <- NULL
for(i in 1:length(refseqlist)){
  tempval <- unlist(refseqlist[[1]][which(unlist(refseqlist[[i]][,1]) %in% sckos),4])
  tempval[which(is.na(tempval))] <- 0
  medians[i] <- median(tempval)
}

nummissing <- NULL
for(i in 1:length(refseqlist)){
  nummissing[i] <- length(which(is.na(unlist(refseqlist[[1]][which(unlist(refseqlist[[i]][,1]) %in% sckos),4]))))
}

refseqlist <- refseqlist[-which(medians > 1.1 | medians < 1 | is.na(medians))]

# FUNCTIONS FOR PREDICTION TOOL

# read fasta to list

GoFasta <- function(file) {
  inputF <- readLines(file)
  headerlines <- grep(">", inputF)
  headerstring <- gsub(">", "", inputF[headerlines])
  seqlines <- data.frame(name=headerstring, from=headerlines+1, to=c((headerlines-1)[-1], length(inputF)))
  seqstring <- rep(NA, length(headerlines))
  for(i in 1:length(headerlines)) {
    seqstring[i] <- paste(inputF[seqlines$from[i]:seqlines$to[i]], collapse="")
  }
  output <- as.list(seqstring)
  names(output) <- headerstring
  return(output)
}

# function to return averaged metagenome of list of genomes
# x is input as a list, commoncol1 is grouping of gene/ko/ec and assumed to be the first column
# commoncol2 will grep for column names containing this string (ie copy number)

AveraGenome <- function(x, commoncol1, commoncol2){
  # merge list and group by gene name
  x <- reduce(x, full_join, by = commoncol1)
  # store first column
  storecol <- x[,1]
  # grep to copy number columns
  x <- x[, grep(commoncol2, colnames(x))]
  # store number of columns for final division
  divideby <- ncol(x)
  # replace nas with 0, sum up copy nos, divide by number
  x[is.na(x)] <- 0
  x <- as.data.frame(rowSums(x))
  x <- x/divideby
  x[,2] <- x[,1]
  x[,1] <- storecol
  colnames(x) <- c(commoncol1, commoncol2)
  return(x)
}

SumGenome <- function(x, commoncol1, commoncol2){
  # merge list and group by gene name
  x <- reduce(x, full_join, by = commoncol1)
  # store first column
  storecol <- x[,1]
  # grep to copy number columns
  x <- x[, grep(commoncol2, colnames(x))]
  # replace nas with 0, sum up copy nos, divide by number
  x[is.na(x)] <- 0
  x <- as.data.frame(rowSums(x))
  x[,2] <- x[,1]
  x[,1] <- storecol
  colnames(x) <- c(commoncol1, commoncol2)
  return(x)
}

allbac <- read_tsv("C:/16s/blasttest/allbound.txt")
#copynotable <- readRDS("C:/16s/blasttest/copynumbertable.rds")

# TO DO: makedb function call from within r

# WIP: generate fasta for makeblastdb input from filtered list
# Needs to: search for genome id in list, find number of corresponding gene ids from 1-any
# output sequences to DF, first column sequence name (gene id), second column genome id (taxid), third sequence?

j <- 1
for(i in 1:10){
  if(length(which(allbac$genome_id == "637000011")) >1){
    testDf[j+1:(j+length(which(allbac$genome_id == "637000011"))),] <- NA
    testDf[(j+1):(j+length(which(allbac$genome_id == "637000011"))),1] <- t(t(which(allbac$genome_id == "637000011")))
    j <- j+length(which(allbac$genome_id == "637000011"))
  }
  j <- j+1
}

testlist <- list()
twocols <- c("genome_id", "sequence")
genomeid_list <- as.character(unlist(distinct(as.data.frame(allbac$genome_id))))
for(i in 1:10){
  cat(paste("\r", i))
  if(length(which(allbac$genome_id == genomeid_list[i])) >1){
    tempDf <- as.data.frame(allbac[which(allbac$genome_id == genomeid_list[i]),twocols])
    templist <- lapply( split(tempDf,seq_along(tempDf[,1])), as.list)
    names(templist) <- allbac[which(allbac$genome_id == genomeid_list[i]),1]
    
    testlist <- c(testlist, templist)
    templist <- list()
  }
  else{
    templist <- list(as.list(allbac[(which(allbac$genome_id == genomeid_list[i])),twocols]))
    names(templist) <- names(templist) <- allbac[which(allbac$genome_id == genomeid_list[i]),1]
    testlist <- c(testlist, templist)
    templist <- list()
  }
}


MakeDB <- function(path_to_makeblastdb, dbname, genomesfasta){
  path_to_makeblastdb <- gsub("/", "\\\\", path_to_makeblastdb)
  tmpdir <- getwd()
  input <- genomesfasta
  tmpdb <- paste0(tmpdir, "/tmpdb.fa")
  # convert to fasta header format
  for (i in 1:length(input)){
    input[i] <- paste0(">", names(input[i]), "\n", input[i])
  }
  input <- as.character(input)
  write(input, tmpdb)
  
  system("cmd.exe", input = paste0(path_to_makeblastdb, " -in ", tmpdir, "/tmpdb.fa  -dbtype nucl -title \"", dbname, "\" "), intern = TRUE, wait = FALSE)
  print(paste("DB output to", tmpdb))
}

#makeblastdb -in test.fsa -parse_seqids -blastdb_version 5 -taxid_map test_map.txt -title "Cookbook demo" -dbtype prot


# main prediction function, x is your fasta file, inputlist is your list of metagenomes
# metagenome list must match the file used to generate reference database
# this function calls AveraGenome but has hard coded variables for columns named Gene and copy number, change if problems occur


BlastIt <- function(x, inputlist, path_to_blastn, path_to_db, threads, numtargets) {
  
  inputlist <- inputlist
  outlist <- list()
  
  # write fasta object to temporary fasta file for blast interface
  tmpdir <- getwd()
  input <- x
  tmpquery <- paste0(tmpdir, "/tmpquery.fa")
  # insert fasta formatting or query seq ids will be replaced with "query 1" and everything fails
  for (i in 1:length(input)){
    input[i] <- paste0(">", names(input[i]), "\n", input[i])
  }
  input <- as.character(input)
  write(input, tmpquery)
  
  # column names for parsing tabular output
  parsecols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                 "send", "evalue", "bitscore", "qcovs", "qcovhsp")
  
  # main system call to blast, should eventually make each argument module
  # db = fasta database, 
  path_to_blastn <- gsub("/", "\\\\", path_to_blastn)
  cat("\r", paste("BLASTing against reference DB, please wait"))
  out <- system("cmd.exe", input = paste0(path_to_blastn, " -db ", path_to_db, " -query ", tmpdir, "/tmpquery.fa  -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp\" -evalue 1e-06 -num_threads ", threads, " -max_target_seqs ", numtargets), intern = TRUE, wait = FALSE)
  cat("\n", paste("BLAST against reference DB complete"))
  cat("\n")
  # remove ms dos output lines and enframe character output to generate data frame
  out <- out[5:(length(out) -2)]
  out <- out %>%
    tibble::enframe() %>%
    tidyr::separate(col = value, into = parsecols,  sep = "\t", convert = TRUE)
  
  # main prediction loop
  for (i in 1:length(x)){
    cat("\r", paste("Now predicting", i, "of", length(x)))
    
    # match blast alignments for each sequence id in fasta file and subset to the equal best matches
    out2 <- out[which(out$qseqid == names(x[i])),]
    out2 <- out2[which(out2$pident == max(out2$pident)),]
    
    # store equally matched metagenomes in temporary list
    templist <- inputlist[c(which(names(inputlist) %in% out2[which(out2$pident == max(out2$pident)),3][[1]]))]
    
    #initialise dumps to store data
    tempDf <- NULL
    
    if (length(which(out2$pident == max(out2$pident))) >1){
      tempDf <- AveraGenome(templist, "Gene", "copy number")
      outlist[[i]] <- tempDf
      names(outlist)[i] <- names(x)[i]
    } else {
      tryCatch(expr = {tempDf <- cbind(templist[[1]][,1], templist[[1]][,4])
      outlist[[i]] <- tempDf
      names(outlist)[i] <- names(x)[i]},
      
      error = function(e){warning(paste("BLAST failure in sequence named", names(x)[i]))},
      
      finally={}
      )
      # removed to test trycatch tempDf <- cbind(templist[[1]][,1], templist[[1]][,4])
    }
    
    # removed to test trycatch outlist[[i]] <- tempDf
    
    
  }
  #names(outlist) <- names(x)
  
  return(outlist)
  
}

# returns a biom object minus the sequences not matched to a reference metagenome

FlushBiom <- function(biom, predictionslist){
  biom <- biom[which(rownames(biom) %in% names(predictionslist)),]
  return(biom)
}

OutputBiom <- function(biom, predictions){
  biomfileout <- list()
  newlist <- NULL
  allgenenames <- NULL
  tempDf <- NULL
  temporaryDf <- NULL
  tempvec <- NULL
  
  for (i in 1:ncol(biom)){
    # assign sequences with abundance >0 to temp df and subset temp list to only those
    temporaryDf <- biom[which(biom[,i] >0),]
    newlist <- predictions[c(row.names(temporaryDf))]
    
    #loop through list and multiply gene copy number values by abundance
    #nested loops were faster in benchmarking but leave apply approach commented just in case, may not scale
    for (x in 1:length(newlist)){
      newlist[[x]][[2]] <- newlist[[x]][[2]] * temporaryDf[x,i]
    }
    
    #newlist <- Map(c, lapply(seq_along(newlist), function(x) newlist[[x]][1]), lapply(seq_along(newlist), function(x) newlist[[x]][2] * temporaryDf[x,i]))
    #newlist <- lapply(seq_along(newlist), function(x) as.data.frame(newlist[[x]]))
    #names(newlist) <- rownames(temporaryDf)
    
    tempDf <- NULL
    if (length(newlist) >1){
      tempDf <- SumGenome(newlist, "Gene", "copy number")
    } else
      tempDf <- as.data.frame(newlist)
    colnames(tempDf) <- c("Gene", "copy number")
    
    biomfileout[i] <- list(tempDf)
    names(biomfileout)[i] <- colnames(biom)[i]
  }
  
  for (i in 1:length(biomfileout)){
    names(biomfileout[[i]]) <- c("Gene", "copy number")
  }
  
  allKOs <- distinct(data.frame(unlist(c(sapply(biomfileout, "[[", "Gene")))))
  
  tempoutDf <- as.data.frame(matrix(data = NA, nrow = nrow(allKOs), ncol = (length(biomfileout)+1)))
  
  tempDf <- allKOs
  
  tempoutDf[,1] <- unlist(lapply(1:length(allKOs[,1]), function(y){
    tempDf[y,1] <- allKOs[y,1]
  }))
  
  colnames(tempoutDf)[1] <- "KO"
  colnames(tempoutDf)[2:ncol(tempoutDf)] <- c(names(biomfileout[1:length(biomfileout)]))
  
  for (i in 1:((ncol(tempoutDf))-1)){
    tempoutDf[,i+1] <- tempoutDf[,1] %in% biomfileout[[i]][[1]]
    tempoutDf[tempoutDf[,i+1] == TRUE,i+1] <- biomfileout[[i]][[2]]
  }
  
  return(tempoutDf)
}

# TO DO
# output per sequence, and output "per sequence per sample", ie which sequences contribute to what within sample


time1 <- Sys.time()

# read in your fasta file, in this instance we are using 137 poops
# from the human microbiome project

hmp_fasta <- GoFasta("C:/16s/blasttest/hmp_16S_rep_seqs.fasta")

# read in your biom file and assign column containing sequence ids to rownames
# if your biom file has a comment line, use skip command
hmp_biom <- read_tsv(file = "C:/16s/blasttest/hmp_16S.biom.tsv", col_names = TRUE)
hmp_biom <- column_to_rownames(hmp_biom, "#OTU ID")

# input folder paths to your installation of blastn and your precomputed database
# note that the system command invoked requires DOS format only for blastn, this will
# be converted automatically if you have problems, check that first

path_to_blastn <- "C:/blast/bin/blastn.exe"
path_to_db <- "C:/16s/blasttest/refseqswithsingledistinct16s.fa"

# call the main prediction function, in this instance we have set a large number of target sequences
# the larger the number of target sequences requested the longer this will take
# however too few target sequences may cause you to not capture the true amount of 16S overlap

hmp_predictions <- BlastIt(hmp_fasta, refseqlist, path_to_blastn, path_to_db, 4, 500)

# inspecting warnings() will show you which sequences failed to be aligned
# these also need to be removed from the biom file to match

hmp_biom <- FlushBiom(hmp_biom, hmp_predictions)

# now you're ready to construct an output biom, this is massively unoptimised right now

hmp_predictions_output <- OutputBiom(hmp_biom, hmp_predictions)

time2 <- Sys.time()
time2-time1

goldstandard <- read_tsv(file = "C:/16s/blasttest/humann2_ko_unstrat.tsv")
piphillin <- read_tsv(file = "C:/16s/blasttest/hmp_piphillin_ko.tsv")
picrust2 <- read_tsv(file = "C:/16s/blasttest/hmp_picrust2_ko_nsti1.0.tsv")

shared_kos <- intersect(goldstandard$"function", hmp_predictions_output$KO)
shared_cols <- intersect(colnames(hmp_predictions_output),  colnames(goldstandard))


mypredvec <- hmp_predictions_output[which(hmp_predictions_output$KO %in% shared_kos),]
goldsvec <- goldstandard[which(goldstandard$"function" %in% shared_kos),]

mypredvecarranged <- arrange(mypredvec, KO)

mypredvecarranged <- mypredvecarranged[,which(colnames(mypredvecarranged) %in% shared_cols)]
goldsvec <- goldsvec[,which(colnames(goldsvec) %in% shared_cols)]

#cor.test(mypredvec[,1], goldsvec[,1], method = "spearman")
resultslist <- list()
for (i in 1:ncol(mypredvecarranged)){
  vec1 <- c(mypredvecarranged[,i])
  vec2 <- c(unlist(goldsvec[,i]))
  
  resultslist[[i]] <- as.list((cor.test(vec1, vec2, method = "spearman")))
}

allestimates <- distinct(data.frame(unlist(c(sapply(resultslist, "[[", "estimate")))))

shared_kos2 <- intersect(goldstandard$"function", piphillin$feature)
shared_cols2 <- intersect(colnames(piphillin),  colnames(goldstandard))

piphillinvec <- piphillin[which(piphillin$feature %in% shared_kos2),]
goldsvec2 <- goldstandard[which(goldstandard$"function" %in% shared_kos2),]

piphillinvec <- piphillinvec[,which(colnames(piphillinvec) %in% shared_cols2)]
goldsvec2 <- goldsvec2[,which(colnames(goldsvec2) %in% shared_cols2)]

resultslist2 <- list()
for (i in 1:137){
  vec1 <- c(unlist(piphillinvec[,i]))
  vec2 <- c(unlist(goldsvec2[,i]))
  
  resultslist2[[i]] <- as.list((cor.test(vec1, vec2, method = "spearman")))
}

allestimates2 <- distinct(data.frame(unlist(c(sapply(resultslist2, "[[", "estimate")))))

colnames(allestimates)[1] <- "correlation"
colnames(allestimates2)[1] <- "correlation"


boxplotdata <- rbind(allestimates, allestimates2)
boxplotdata[,2] <- boxplotdata[,1]
colnames(boxplotdata)[2] <- "tool"
boxplotdata[1:137,2] <- "Proof of Concept"
boxplotdata[138:274,2] <- "Piphillin"
tool_order <- c("Piphillin", "Proof of Concept") 

ggplot(boxplotdata, aes(x=factor(tool, level = tool_order), y=correlation, fill = tool)) + 
  geom_boxplot(
    
    # Notch?
    notch=TRUE,
    notchwidth = 0.8,
    
    # custom outliers
    #outlier.colour="red",
    #outlier.fill="red",
    #outlier.size=3
    outliers = FALSE
    
  ) + ylab("Spearman Correlation Coefficient") +xlab(NULL) + geom_jitter()

calc_accuracy_metrics <- function(df1, df2, category) {
  
  # Subset only to columns and rows that overlap between both and convert to present (TRUE) and absent (FALSE)
  cols2keep <- colnames(df1)[which(colnames(df1) %in% colnames(df2))]
  rows2keep <- rownames(df1)[which(rownames(df1) %in% rownames(df2))]
  
  
  df1 <- df1[rows2keep, cols2keep, drop=FALSE] > 0
  df2 <- df2[rows2keep, cols2keep, drop=FALSE] > 0
  
  out_df <- data.frame(matrix(NA, nrow=length(cols2keep), ncol=12))
  colnames(out_df) <- c("category", "sample", "acc", "TP", "TN", "FP", "FN", "NPV", "precision", "recall", "fpr", "F1")
  
  row_i = 1
  for(sample in colnames(df1)) {
    
    total_func <- length(rows2keep)
    
    overall_acc <- sum(df1[,sample] == df2[,sample])/total_func
    
    num_true_pos <- length(which(which(df1[,sample]) %in% which(df2[,sample])))
    num_true_neg <- length(which(which(! df1[,sample]) %in% which(! df2[,sample])))
    
    num_false_pos <- length(which(which(! df1[,sample]) %in% which(df2[,sample])))
    num_false_neg <- length(which(which(df1[,sample]) %in% which(! df2[,sample])))
    
    npv <- num_true_neg/(num_true_neg + num_false_neg)
    precision <- num_true_pos/(num_true_pos + num_false_pos)
    recall <- num_true_pos/(num_true_pos + num_false_neg)
    fpr <- num_false_pos/(num_false_pos + num_true_neg)
    F1 <-  2 * ((precision * recall)/(precision + recall))
    
    out_df[row_i, ] <- c(NA, NA, overall_acc, num_true_pos, num_true_neg, num_false_pos, num_false_neg,
                         npv, precision, recall, fpr, F1)
    
    row_i = row_i + 1
  }
  
  out_df$category <- category  
  out_df$sample <- cols2keep
  
  return(out_df)
}

hmp_poc <- calc_accuracy_metrics(goldstandard, hmp_predictions_output, "HMP_POC")
hmp_piphillin <- calc_accuracy_metrics(goldstandard, piphillin, "HMP_piphillin")
hmp_picrust <- calc_accuracy_metrics(goldstandard, picrust2, "HMP_picrust2")






















Confused <- function(Actual,Predicted){
  actual = as.data.frame(table(Actual))
  names(actual) = c("Actual","ActualFreq")
  
  #build confusion matrix
  confusion = as.data.frame(table(Actual, Predicted))
  names(confusion) = c("Actual","Predicted","Freq")
  
  #calculate percentage of test cases based on actual frequency
  
  confusion = merge(confusion, actual, by=c('Actual','Actual'))
  confusion$Percent = confusion$Freq/confusion$ActualFreq*100
  confusion$ColorScale<-confusion$Percent*-1
  confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale<-confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale*-1
  confusion$Label<-paste(round(confusion$Percent,0),"%, n=",confusion$Freq,sep="")
  
  ggplot(confusion) +
    geom_tile(aes(x=Actual, y=Predicted,fill=ColorScale), data=confusion, color="black",size=0.1) + coord_equal() +
    labs(x="Actual",y="Predicted") +
    
    geom_text(aes(x=Actual,y=Predicted, label=Label),data=confusion, size=5, colour="black") +
    scale_fill_gradient2(low="red4",high="steelblue",mid="white", midpoint = 0,guide='none')
  
  #return(confusion)
}

Actual <- NULL
Predicted <- NULL

SumAll <- sum(hmp_poc$TP, hmp_poc$TN, hmp_poc$FP, hmp_poc$FN)
SumTP <- sum(hmp_poc$TP)
SumFP <- sum(hmp_poc$FP)
SumFN <- sum(hmp_poc$FN)

Actual[1:SumAll] <- 0
Predicted[1:SumAll] <- 0
Actual[1:SumTP] <- 1
Predicted[1:sum(SumTP, SumFP)] <- 1
Actual[((SumAll-SumFN)+1):SumAll] <- 1

Confused(Actual, Predicted)












GoFasta <- function(file) {
  inputF <- readLines(file)
  headerlines <- grep(">", inputF)
  headerstring <- gsub(">", "", inputF[headerlines])
  seqlines <- data.frame(name=headerstring, from=headerlines+1, to=c((headerlines-1)[-1], length(inputF)))
  seqstring <- rep(NA, length(headerlines))
  for(i in 1:length(headerlines)) {
    seqstring[i] <- paste(inputF[seqlines$from[i]:seqlines$to[i]], collapse="")
  }
  output <- as.list(seqstring)
  names(output) <- headerstring
  return(output)
}

# function to return averaged metagenome of list of genomes
# x is input as a list, commoncol1 is grouping of gene/ko/ec and assumed to be the first column
# commoncol2 will grep for column names containing this string (ie copy number)

AveraGenome <- function(x, commoncol1, commoncol2){
  # merge list and group by gene name
  x <- reduce(x, full_join, by = commoncol1)
  # store first column
  storecol <- x[,1]
  # grep to copy number columns
  x <- x[, grep(commoncol2, colnames(x))]
  # store number of columns for final division
  divideby <- ncol(x)
  # replace nas with 0, sum up copy nos, divide by number
  x[is.na(x)] <- 0
  x <- as.data.frame(rowSums(x))
  x <- x/divideby
  x[,2] <- x[,1]
  x[,1] <- storecol
  colnames(x) <- c(commoncol1, commoncol2)
  return(x)
}

SumGenome <- function(x, commoncol1, commoncol2){
  # merge list and group by gene name
  x <- reduce(x, full_join, by = commoncol1)
  # store first column
  storecol <- x[,1]
  # grep to copy number columns
  x <- x[, grep(commoncol2, colnames(x))]
  # replace nas with 0, sum up copy nos, divide by number
  x[is.na(x)] <- 0
  x <- as.data.frame(rowSums(x))
  x[,2] <- x[,1]
  x[,1] <- storecol
  colnames(x) <- c(commoncol1, commoncol2)
  return(x)
}


MakeDB <- function(path_to_makeblastdb, dbname, genomesfasta){
  path_to_makeblastdb <- gsub("/", "\\\\", path_to_makeblastdb)
  tmpdir <- getwd()
  input <- genomesfasta
  tmpdb <- paste0(tmpdir, "/", dbname, ".fa")
  # convert to fasta header format
  for (i in 1:length(input)){
    input[i] <- paste0(">", names(input[i]), "\n", input[i])
  }
  input <- as.character(input)
  write(input, tmpdb)
  
  system("cmd.exe", input = paste0(path_to_makeblastdb, " -in ", tmpdb, "  -dbtype nucl -title \"", dbname, "\" "), intern = TRUE, wait = FALSE)
  print(paste("DB output to", tmpdb))
}

BlastIt <- function(x, inputlist, path_to_blastn, path_to_db, threads, numtargets) {
  
  inputlist <- inputlist
  outlist <- list()
  
  # write fasta object to temporary fasta file for blast interface
  tmpdir <- getwd()
  input <- x
  tmpquery <- paste0(tmpdir, "/tmpquery.fa")
  # insert fasta formatting or query seq ids will be replaced with "query 1" and everything fails
  for (i in 1:length(input)){
    input[i] <- paste0(">", names(input[i]), "\n", input[i])
  }
  input <- as.character(input)
  write(input, tmpquery)
  
  # column names for parsing tabular output
  parsecols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                 "send", "evalue", "bitscore", "qcovs", "qcovhsp")
  
  # main system call to blast, should eventually make each argument module
  # db = fasta database, 
  path_to_blastn <- gsub("/", "\\\\", path_to_blastn)
  cat("\r", paste("BLASTing against reference DB, please wait"))
  out <- system("cmd.exe", input = paste0(path_to_blastn, " -db ", path_to_db, " -query ", tmpdir, "/tmpquery.fa  -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp\" -evalue 1e-06 -num_threads ", threads, " -max_target_seqs ", numtargets), intern = TRUE, wait = FALSE)
  cat("\n", paste("BLAST against reference DB complete"))
  cat("\n")
  # remove ms dos output lines and enframe character output to generate data frame
  out <- out[5:(length(out) -2)]
  out <- out %>%
    tibble::enframe() %>%
    tidyr::separate(col = value, into = parsecols,  sep = "\t", convert = TRUE)
  
  # main prediction loop
  for (i in 1:length(x)){
    cat("\r", paste("Now predicting", i, "of", length(x)))
    
    # match blast alignments for each sequence id in fasta file and subset to the equal best matches
    out2 <- out[which(out$qseqid == names(x[i])),]
    out2 <- out2[which(out2$pident == max(out2$pident)),]
    
    # store equally matched metagenomes in temporary list
    templist <- inputlist[c(which(names(inputlist) %in% out2[which(out2$pident == max(out2$pident)),3][[1]]))]
    
    #initialise dumps to store data
    tempDf <- NULL
    
    if (length(which(out2$pident == max(out2$pident))) >1){
      tempDf <- AveraGenome(templist, "Gene", "copy number")
      outlist[[i]] <- tempDf
      names(outlist)[i] <- names(x)[i]
    } else {
      tryCatch(expr = {tempDf <- cbind(templist[[1]][,1], templist[[1]][,4])
      outlist[[i]] <- tempDf
      names(outlist)[i] <- names(x)[i]},
      
      error = function(e){warning(paste("BLAST failure in sequence named", names(x)[i]))},
      
      finally={}
      )
      # removed to test trycatch tempDf <- cbind(templist[[1]][,1], templist[[1]][,4])
    }
    
    # removed to test trycatch outlist[[i]] <- tempDf
    
    
  }
  #names(outlist) <- names(x)
  
  return(outlist)
  
}

# returns a biom object minus the sequences not matched to a reference metagenome

FlushBiom <- function(biom, predictionslist){
  biom <- biom[which(rownames(biom) %in% names(predictionslist)),]
  return(biom)
}

OutputBiom <- function(biom, predictions){
  biomfileout <- list()
  newlist <- NULL
  allgenenames <- NULL
  tempDf <- NULL
  temporaryDf <- NULL
  tempvec <- NULL
  
  for (i in 1:ncol(biom)){
    # assign sequences with abundance >0 to temp df and subset temp list to only those
    temporaryDf <- biom[which(biom[,i] >0),]
    newlist <- predictions[c(row.names(temporaryDf))]
    
    #loop through list and multiply gene copy number values by abundance
    #nested loops were faster in benchmarking but leave apply approach commented just in case, may not scale
    for (x in 1:length(newlist)){
      newlist[[x]][[2]] <- newlist[[x]][[2]] * temporaryDf[x,i]
    }
    
    #newlist <- Map(c, lapply(seq_along(newlist), function(x) newlist[[x]][1]), lapply(seq_along(newlist), function(x) newlist[[x]][2] * temporaryDf[x,i]))
    #newlist <- lapply(seq_along(newlist), function(x) as.data.frame(newlist[[x]]))
    #names(newlist) <- rownames(temporaryDf)
    
    tempDf <- NULL
    if (length(newlist) >1){
      tempDf <- SumGenome(newlist, "Gene", "copy number")
    } else
      tempDf <- as.data.frame(newlist)
    colnames(tempDf) <- c("Gene", "copy number")
    
    biomfileout[i] <- list(tempDf)
    names(biomfileout)[i] <- colnames(biom)[i]
  }
  
  for (i in 1:length(biomfileout)){
    names(biomfileout[[i]]) <- c("Gene", "copy number")
  }
  
  allKOs <- distinct(data.frame(unlist(c(sapply(biomfileout, "[[", "Gene")))))
  
  tempoutDf <- as.data.frame(matrix(data = NA, nrow = nrow(allKOs), ncol = (length(biomfileout)+1)))
  
  tempDf <- allKOs
  
  tempoutDf[,1] <- unlist(lapply(1:length(allKOs[,1]), function(y){
    tempDf[y,1] <- allKOs[y,1]
  }))
  
  colnames(tempoutDf)[1] <- "KO"
  colnames(tempoutDf)[2:ncol(tempoutDf)] <- c(names(biomfileout[1:length(biomfileout)]))
  
  for (i in 1:((ncol(tempoutDf))-1)){
    tempoutDf[,i+1] <- tempoutDf[,1] %in% biomfileout[[i]][[1]]
    tempoutDf[tempoutDf[,i+1] == TRUE,i+1] <- biomfileout[[i]][[2]]
  }
  
  return(tempoutDf)
}

refseqlist <- readRDS("C:/16s/blasttest/refseqlist.rds")

for(i in 1:length(refseqlist)){
  names(refseqlist[[i]])[1] <- "Gene"
  names(refseqlist[[i]])[2] <- "Name"
  names(refseqlist[[i]])[3] <- "Definition"
  names(refseqlist[[i]])[4] <- "copy number"
}

sckos <- c("K00057","K00059","K00075","K00088","K00133","K00134","K00384","K00525","K00554","K00566","K00600","K00604","K00609","K00615","K00655","K00762","K00789","K00790","K00791","K00820","K00858","K00859","K00927","K00939","K00942","K00948","K00954","K00962","K00981","K01000","K01056","K01265","K01358","K01409","K01462","K01465","K01491","K01591","K01689","K01714","K01756","K01775","K01783","K01803","K01810","K01866","K01867","K01868","K01869","K01870","K01872","K01873","K01874","K01875","K01876","K01881","K01883","K01887","K01889","K01890","K01892","K01915","K01921","K01924","K01925","K01929","K01937","K01939","K01945","K01951","K01955","K01956","K01972","K01990","K01992","K02003","K02004","K02078",
           "K02108","K02109","K02110","K02111","K02112","K02114","K02115","K02313","K02314","K02316","K02335","K02337","K02338","K02340","K02341","K02343","K02355","K02356","K02357","K02358","K02469","K02470","K02493","K02495","K02518","K02519","K02520","K02528","K02563","K02600","K02601","K02834","K02835","K02836","K02838","K02860","K02863","K02864","K02867","K02871","K02874","K02876","K02878","K02879","K02881","K02884","K02886","K02887","K02888","K02890","K02892","K02895","K02899","K02902","K02904","K02906","K02909","K02911","K02913","K02916","K02926","K02931","K02933","K02935","K02939","K02945","K02946","K02948","K02950","K02952","K02954","K02956","K02959","K02961","K02963","K02965","K02967","K02968","K02982","K02986","K02988","K02990","K02992","K02994","K02996","K03040","K03043","K03046","K03070","K03075","K03076","K03086","K03100","K03101","K03106","K03110","K03111","K03168","K03177","K03217","K03218","K03424","K03438","K03466","K03470","K03501","K03530","K03531","K03544","K03545","K03550","K03551","K03553","K03588","K03595","K03596","K03624","K03625","K03631","K03655","K03657","K03664","K03671","K03685","K03686","K03687","K03701","K03702","K03703","K03723","K03798","K03977",
           "K03979","K04043","K04066","K04075","K04077","K04078","K04096","K04485","K04487","K06173","K06178","K06180","K06187","K06207","K06942","K07056","K07447","K07478","K09710","K09761","K09903","K10773","K11753","K11754")

medians <- NULL
for(i in 1:length(refseqlist)){
  medians[i] <- median(unlist(refseqlist[[1]][which(unlist(refseqlist[[i]][,1]) %in% sckos),4]))
}

allbac <- read_tsv("C:/16s/blasttest/allbound.txt")
allbac2 <- allbac[-which(allbac$negative == 1),]

refseqlist <- refseqlist[-which(medians > 1.1 | medians < 1 | is.na(medians))]
allbac2 <- allbac2[which(allbac2$genome_id %in% names(refseqlist)),]

seqs <- as.vector(allbac2$sequence)
names(seqs) <- allbac2$genome_id

setwd("C:/16s/blasttest/")
MakeDB("C:/Blast/bin/makeblastdb.exe", dbname = "example", seqs)

# refseqlist object and db must have same genome ids

time1 <- Sys.time()

# read in your fasta file, in this instance we are using 137 poops
# from the human microbiome project

hmp_fasta <- GoFasta("C:/16s/blasttest/hmp_16S_rep_seqs.fasta")

# read in your biom file and assign column containing sequence ids to rownames
# if your biom file has a comment line, use skip command
hmp_biom <- read_tsv(file = "C:/16s/blasttest/hmp_16S.biom.tsv", col_names = TRUE)
hmp_biom <- column_to_rownames(hmp_biom, "#OTU ID")

# input folder paths to your installation of blastn and your precomputed database
# note that the system command invoked requires DOS format only for blastn, this will
# be converted automatically if you have problems, check that first

path_to_blastn <- "C:/blast/bin/blastn.exe"
path_to_db <- "C:/16s/blasttest/example.fa"

# call the main prediction function, in this instance we have set a large number of target sequences
# the larger the number of target sequences requested the longer this will take
# however too few target sequences may cause you to not capture the true amount of 16S overlap

hmp_predictions <- BlastIt(hmp_fasta, refseqlist, path_to_blastn, path_to_db, 4, 500)

# inspecting warnings() will show you which sequences failed to be aligned
# these also need to be removed from the biom file to match

hmp_biom <- FlushBiom(hmp_biom, hmp_predictions)

# now you're ready to construct an output biom, this is massively unoptimised right now

hmp_predictions_output <- OutputBiom(hmp_biom, hmp_predictions)

time2 <- Sys.time()
time2-time1

goldstandard <- read_tsv(file = "C:/16s/blasttest/humann2_ko_unstrat.tsv")
piphillin <- read_tsv(file = "C:/16s/blasttest/hmp_piphillin_ko.tsv")
picrust2 <- read_tsv(file = "C:/16s/blasttest/hmp_picrust2_ko_nsti1.0.tsv")

shared_kos <- intersect(goldstandard$"function", hmp_predictions_output$KO)
shared_cols <- intersect(colnames(hmp_predictions_output),  colnames(goldstandard))


mypredvec <- hmp_predictions_output[which(hmp_predictions_output$KO %in% shared_kos),]
goldsvec <- goldstandard[which(goldstandard$"function" %in% shared_kos),]

mypredvecarranged <- arrange(mypredvec, KO)

mypredvecarranged <- mypredvecarranged[,which(colnames(mypredvecarranged) %in% shared_cols)]
goldsvec <- goldsvec[,which(colnames(goldsvec) %in% shared_cols)]

#cor.test(mypredvec[,1], goldsvec[,1], method = "spearman")
resultslist <- list()
for (i in 1:ncol(mypredvecarranged)){
  vec1 <- c(mypredvecarranged[,i])
  vec2 <- c(unlist(goldsvec[,i]))
  
  resultslist[[i]] <- as.list((cor.test(vec1, vec2, method = "spearman")))
}

allestimates <- distinct(data.frame(unlist(c(sapply(resultslist, "[[", "estimate")))))

shared_kos2 <- intersect(goldstandard$"function", piphillin$feature)
shared_cols2 <- intersect(colnames(piphillin),  colnames(goldstandard))

piphillinvec <- piphillin[which(piphillin$feature %in% shared_kos2),]
goldsvec2 <- goldstandard[which(goldstandard$"function" %in% shared_kos2),]

piphillinvec <- piphillinvec[,which(colnames(piphillinvec) %in% shared_cols2)]
goldsvec2 <- goldsvec2[,which(colnames(goldsvec2) %in% shared_cols2)]

resultslist2 <- list()
for (i in 1:137){
  vec1 <- c(unlist(piphillinvec[,i]))
  vec2 <- c(unlist(goldsvec2[,i]))
  
  resultslist2[[i]] <- as.list((cor.test(vec1, vec2, method = "spearman")))
}

allestimates2 <- distinct(data.frame(unlist(c(sapply(resultslist2, "[[", "estimate")))))

colnames(allestimates)[1] <- "correlation"
colnames(allestimates2)[1] <- "correlation"


boxplotdata <- rbind(allestimates, allestimates2)
boxplotdata[,2] <- boxplotdata[,1]
colnames(boxplotdata)[2] <- "tool"
boxplotdata[1:137,2] <- "Proof of Concept"
boxplotdata[138:274,2] <- "Piphillin"
tool_order <- c("Piphillin", "Proof of Concept") 

ggplot(boxplotdata, aes(x=factor(tool, level = tool_order), y=correlation, fill = tool)) + 
  geom_boxplot(
    
    # Notch?
    notch=TRUE,
    notchwidth = 0.8,
    
    # custom outliers
    #outlier.colour="red",
    #outlier.fill="red",
    #outlier.size=3
    outliers = FALSE
    
  ) + ylab("Spearman Correlation Coefficient") +xlab(NULL) + geom_jitter()

library(splinectomeR)

dbnames <- c("2500db", "5000db", "7500db", "10000db", "12500db", "15000db", "17500db", "20000db", "22500db", "25000db", "27500db", "30000db", "30775db")
spanlist <- list()
splinelist <- list()
for(i in 1:length(dbnames)){
  path_to_blastn <- "C:/blast/bin/blastn.exe"
  path_to_db <- paste0("C:/16s/blasttest/", dbnames[i], ".fa")
  
  # call the main prediction function, in this instance we have set a large number of target sequences
  # the larger the number of target sequences requested the longer this will take
  # however too few target sequences may cause you to not capture the true amount of 16S overlap
  
  spanlist[[i]] <- BlastIt(hmp_fasta, refseqlist, path_to_blastn, path_to_db, 4, 500)
  hmp_biom <- FlushBiom(hmp_biom, spanlist[[i]])
  
  # now you're ready to construct an output biom, this is massively unoptimised right now
  
  spanlist[[i]] <- OutputBiom(hmp_biom, spanlist[[i]])
  
  splinelist[[i]] <- calc_accuracy_metrics(goldstandard, spanlist[[i]], "HMP_POC")
}

tempres <- data.frame(acc = as.vector(sapply(splinelist, "[[", 3)), TP = as.vector(sapply(splinelist, "[[", 4)), TN = as.vector(sapply(splinelist, "[[", 5)), FP = as.vector(sapply(splinelist, "[[", 6)), FN = as.vector(sapply(splinelist, "[[", 7)), NPV = as.vector(sapply(splinelist, "[[", 8)), precision = as.vector(sapply(splinelist, "[[", 9)), recall = as.vector(sapply(splinelist, "[[", 10)), fpr = as.vector(sapply(splinelist, "[[", 11)), F1 = as.vector(sapply(splinelist, "[[", 12)), dbsize = as.numeric(rep(sapply(strsplit(dbnames, "db"), "[", 1), each = 137)))

# construct ordered list by genome size

genomesize_ord <- order(unlist(lapply(refseqlist, nrow)))

allbac2 <- allbac[-which(allbac$negative == 1),]

for(i in 1:length(dbnames)){
  seqs <- allbac2$sequence[which(allbac2$genome_id %in% names(refseqlist[genomesize_ord])[1:as.numeric(sapply(strsplit(dbnames, "db"), "[", 1))[i]])]
  names(seqs) <- allbac2$genome_id[which(allbac2$genome_id %in% names(refseqlist[genomesize_ord])[1:as.numeric(sapply(strsplit(dbnames, "db"), "[", 1))[i]])]
    MakeDB("C:/Blast/bin/makeblastdb.exe", dbname = paste0("ordered", as.numeric(sapply(strsplit(dbnames, "db"), "[", 1))[i]), seqs)
  }

dbnames_ord <- c("2500", "5000", "7500", "10000", "12500", "15000", "17500", "20000", "22500", "25000", "27500", "30000", "30775")
spanlist_ord <- list()
splinelist_ord <- list()
for(i in 1:length(dbnames)){
  path_to_blastn <- "C:/blast/bin/blastn.exe"
  path_to_db <- paste0("C:/16s/blasttest/ordered", dbnames_ord[i], ".fa")
  refseqlist2 <- refseqlist[which(names(refseqlist) %in% names(refseqlist[genomesize_ord])[1:as.numeric(sapply(strsplit(dbnames, "db"), "[", 1))[i]])]
  # call the main prediction function, in this instance we have set a large number of target sequences
  # the larger the number of target sequences requested the longer this will take
  # however too few target sequences may cause you to not capture the true amount of 16S overlap
  
  spanlist_ord[[i]] <- BlastIt(hmp_fasta, refseqlist2, path_to_blastn, path_to_db, 4, 500)
  hmp_biom <- read_tsv(file = "C:/16s/blasttest/hmp_16S.biom.tsv", col_names = TRUE)
  hmp_biom <- column_to_rownames(hmp_biom, "#OTU ID")
  
  hmp_biom <- FlushBiom(hmp_biom, spanlist_ord[[i]])
  
  # now you're ready to construct an output biom, this is massively unoptimised right now
  
  spanlist_ord[[i]] <- OutputBiom(hmp_biom, spanlist_ord[[i]])
  
  splinelist_ord[[i]] <- calc_accuracy_metrics(goldstandard, spanlist_ord[[i]], "HMP_POC")
}

res_ord <- data.frame(acc = as.vector(sapply(splinelist_ord, "[[", 3)), TP = as.vector(sapply(splinelist_ord, "[[", 4)), TN = as.vector(sapply(splinelist_ord, "[[", 5)), FP = as.vector(sapply(splinelist_ord, "[[", 6)), FN = as.vector(sapply(splinelist_ord, "[[", 7)), NPV = as.vector(sapply(splinelist_ord, "[[", 8)), precision = as.vector(sapply(splinelist_ord, "[[", 9)), recall = as.vector(sapply(splinelist_ord, "[[", 10)), fpr = as.vector(sapply(splinelist_ord, "[[", 11)), F1 = as.vector(sapply(splinelist_ord, "[[", 12)), dbsize = as.numeric(rep(sapply(strsplit(dbnames, "db"), "[", 1), each = 137)))

ggplot(res_ord, aes(x = dbsize, y = F1)) +
  geom_point() +
  geom_smooth()

allbac2 <- allbac[-which(allbac$negative == 1),]
tempind <- list()
set.seed(12345)
for(i in 1:length(dbnames)){
  tempind[[i]] <- sample(names(refseqlist), as.numeric(sapply(strsplit(dbnames, "db"), "[", 1))[i])
  seqs <- allbac2$sequence[which(allbac2$genome_id %in% tempind[[i]])] 
  names(seqs) <- allbac2$genome_id[which(allbac2$genome_id %in% tempind[[i]])]
  MakeDB("C:/Blast/bin/makeblastdb.exe", dbname = paste0("random", as.numeric(sapply(strsplit(dbnames, "db"), "[", 1))[i]), seqs)
}

dbnames_ord <- c("2500", "5000", "7500", "10000", "12500", "15000", "17500", "20000", "22500", "25000", "27500", "30000", "30775")
spanlist_rand <- list()
splinelist_rand <- list()
for(i in 1:length(dbnames)){
  path_to_blastn <- "C:/blast/bin/blastn.exe"
  path_to_db <- paste0("C:/16s/blasttest/random", dbnames_ord[i], ".fa")
  refseqlist2 <- refseqlist[which(names(refseqlist) %in% names(refseqlist[genomesize_ord])[1:as.numeric(sapply(strsplit(dbnames, "db"), "[", 1))[i]])]
  # call the main prediction function, in this instance we have set a large number of target sequences
  # the larger the number of target sequences requested the longer this will take
  # however too few target sequences may cause you to not capture the true amount of 16S overlap
  
  spanlist_rand[[i]] <- BlastIt(hmp_fasta, refseqlist, path_to_blastn, path_to_db, 4, 500)
  hmp_biom <- read_tsv(file = "C:/16s/blasttest/hmp_16S.biom.tsv", col_names = TRUE)
  hmp_biom <- column_to_rownames(hmp_biom, "#OTU ID")
  
  hmp_biom <- FlushBiom(hmp_biom, spanlist_rand[[i]])
  
  # now you're ready to construct an output biom, this is massively unoptimised right now
  
  spanlist_rand[[i]] <- OutputBiom(hmp_biom, spanlist_rand[[i]])
  
  splinelist_rand[[i]] <- calc_accuracy_metrics(goldstandard, spanlist_rand[[i]], "HMP_POC")
}

res_rand <- data.frame(acc = as.vector(sapply(splinelist_rand, "[[", 3)), TP = as.vector(sapply(splinelist_rand, "[[", 4)), TN = as.vector(sapply(splinelist_rand, "[[", 5)), FP = as.vector(sapply(splinelist_rand, "[[", 6)), FN = as.vector(sapply(splinelist_rand, "[[", 7)), NPV = as.vector(sapply(splinelist_rand, "[[", 8)), precision = as.vector(sapply(splinelist_rand, "[[", 9)), recall = as.vector(sapply(splinelist_rand, "[[", 10)), fpr = as.vector(sapply(splinelist_rand, "[[", 11)), F1 = as.vector(sapply(splinelist_rand, "[[", 12)), dbsize = as.numeric(rep(sapply(strsplit(dbnames, "db"), "[", 1), each = 137)))

res_both <- rbind(res_ord, res_rand)
res_both$type <- c(rep("Ordered", nrow(res_ord)), rep("Random", nrow(res_rand)))

ggplot(res_both, aes(x = dbsize, y = F1)) +
  geom_point() + facet_wrap(vars(type)) +
  geom_smooth()

# add for splinectomer
res_both$sample <- c(rep(c(paste0(shared_cols, "ord")), 13), rep(c(paste0(shared_cols, "rand")), 13))

# repeat for all indices
result <- trendyspliner(res_both, xvar = "dbsize", yvar = "F1", category = "type", cases = "sample", perms = 999)

# plots all permuted splines
allplot <- function (data = NULL, xvar = NULL, yvar = NULL) 
{
  if (is.null(data) | is.null(xvar) | is.null(yvar)) {
    stop("Missing required arugments.")
  }
  if (is.null(data["permuted_splines"][[1]])) {
    stop("Permuted data not in results. Did you set \"retain_perms = TRUE\" in the permuspliner() run?")
  }
  require(ggplot2)
  require(reshape2)
  permsplines <- data["permuted_splines"][[1]]
  permsplines <- permsplines[, grep("perm", colnames(permsplines))]
  num_perms <- (ncol(permsplines)/2)
  permsplines$x.par <- rownames(permsplines)
  rownames(permsplines) <- NULL
  permsplines <- melt(permsplines, id.vars = "x.par", variable.name = "permutation", 
                      value.name = "y.par")
  var_1 <- as.character(data["category_1"][[1]])
  var_2 <- as.character(data["category_2"][[1]])
  true_v1 <- data["v1_interpolated"][[1]]
  true_v1$var <- var_1
  colnames(true_v1)[2] <- "y"
  true_v2 <- data["v2_interpolated"][[1]]
  true_v2$var <- var_2
  colnames(true_v2)[2] <- "y"
  true_data <- rbind(true_v1, true_v2)
  var_labels <- factor(true_data$var, ordered = T)
  true_data$var <- var_labels
  num_points <- length(true_v1$x)
  colvec <- vector(length = nrow(permsplines))
  colvec[which(lengths(strsplit(as.character(permsplines$permutation), "v2")) == 1)] <- 1
  colvec[which(lengths(strsplit(as.character(permsplines$permutation), "v2")) != 1)] <- 2
  colvec <- as.factor(colvec)
  p <- ggplot() + geom_line(data = permsplines, aes(x = as.numeric(x.par), 
                                                    y = as.numeric(y.par), group = factor(permutation), color = colvec),, size = 1, na.rm = T) + scale_color_manual(name = "", values = c("darkorange", 
                                                                                                                 "steelblue")) + theme_classic() + theme(axis.text = element_text(color = "black")) + 
    xlab(xvar) + ylab(yvar)
  
  return(p)
}

allplot(result, "dbsize", "F1")

# plots single grouped permuted splines
groupplot <- function (data = NULL, xvar = NULL, yvar = NULL) 
{
  if (is.null(data) | is.null(xvar) | is.null(yvar)) {
    stop("Missing required arugments.")
  }
  if (is.null(data["permuted_splines"][[1]])) {
    stop("Permuted data not in results. Did you set \"retain_perms = TRUE\" in the permuspliner() run?")
  }
  require(ggplot2)
  require(reshape2)
  permsplines <- data["permuted_splines"][[1]]
  permsplines <- permsplines[, grep("perm", colnames(permsplines))]
  num_perms <- (ncol(permsplines)/2)
  permsplines$x.par <- rownames(permsplines)
  rownames(permsplines) <- NULL
  permsplines <- melt(permsplines, id.vars = "x.par", variable.name = "permutation", 
                      value.name = "y.par")
  var_1 <- as.character(data["category_1"][[1]])
  var_2 <- as.character(data["category_2"][[1]])
  true_v1 <- data["v1_interpolated"][[1]]
  true_v1$var <- var_1
  colnames(true_v1)[2] <- "y"
  true_v2 <- data["v2_interpolated"][[1]]
  true_v2$var <- var_2
  colnames(true_v2)[2] <- "y"
  true_data <- rbind(true_v1, true_v2)
  var_labels <- factor(true_data$var, ordered = T)
  true_data$var <- var_labels
  num_points <- length(true_v1$x)
  p <- ggplot() + geom_line(data = true_data, 
                            aes(x = as.numeric(x), y = as.numeric(y), color = var_labels), 
                            size = 1.2) + scale_color_manual(name = "", values = c("darkorange", 
                                                                                   "steelblue")) + theme_classic() + theme(axis.text = element_text(color = "black")) + 
    xlab(xvar) + ylab(yvar)
  return(p)
}

# repeat for all indices, use patchwork to construct 4 panel plots
allplot(result, "dbsize", "F1")
# trend plot, repeat for all indices

result_p <- sliding_spliner(data = res_both, xvar = 'dbsize', yvar = 'F1', category = 'type', groups = c('Ordered','Random'), cases = 'sample', ints = 100)

sliding_spliner.plot.pvals(result_p)