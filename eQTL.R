setwd("D:/Study thingie/FACH/S5 WS2021/GOBI/eQTL/input_exercise")
library(data.table)
library(R.utils)
library(ggplot2)
library(preprocessCore)
library(reshape2)
library(tidyverse)
library(MatrixEQTL)

#1
example <-readRDS("eqtl_examples.RDS")
cor.test(example$gene1,example$snp1)
#t = -25.01, df = 443, p-value < 2.2e-16 cor = -0.7651
example$snp1 <- as.character(example$snp1)
ggplot(example, aes(snp1,gene1,color = snp1)) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                outlier.size=2, notch=FALSE)+facet_wrap(~gender)
View(example)
summary(lm(example$gene1~example$snp1 +example$gender))

#2
gd660 <- fread("GD660.GeneQuantRPKM.txt.gz")
head(gd660)
matrix1 <- gd660[,1:4]
matrix2 <- gd660[,-(1:4)]
matrix2 <- data.matrix(matrix2)

row.names(matrix2)<- matrix1$TargetID
convered_matrix2 <- matrix2[rowMeans(matrix2!=0) >0.5, ]

convered_matrix2[convered_matrix2 == 0] <- 1

convered_matrix2 <- log2(convered_matrix2)
sum(convered_matrix2==0)

#for gene plotting
cut_matrix2 <- convered_matrix2[1:50,]
cut_matrix2 <- data.frame(cut_matrix2)



unlist_matrix2 <- unlist(cut_matrix2)
unlist_matrix2 <- data.matrix(unlist_matrix2)
unlist_matrix2 <- data.frame(unlist_matrix2)
unlist_matrix2$gene_name <- rownames(unlist_matrix2) 

trans_unlist_matrix2 <- data.frame(t(cut_matrix2))
trans_unlist_matrix2 <- stack(trans_unlist_matrix2)



ggplot(unlist_matrix2,aes(x = unlist_matrix2))+
  geom_histogram()+
  labs(x = "values")

ggplot(trans_unlist_matrix2,aes(ind, values))+
  geom_boxplot()+
  labs(x = "gene_names", y = "values")+
  theme(axis.text.x = element_text(angle = 90))

#for samples plotting
samples_matrix2 <- convered_matrix2[,1:50]
samples_matrix2 <- data.frame(samples_matrix2)
stack_matrix2 <- stack(samples_matrix2)
ggplot(stack_matrix2, aes(x = ind,y = values))+ geom_boxplot()+theme(axis.text.x = element_text(angle = 90))

#quantail normarlization

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
nq_matrix2 <- quantile_normalisation(convered_matrix2)

nq <- nq_matrix2[,1:50]
nq <- data.frame(nq)
nq <- stack(nq)
ggplot(nq, aes(x = ind,y = values))+ geom_boxplot()+ theme(axis.text.x = element_text(angle = 90))


colnames(convered_matrix2) <- gsub("^([^.]+).*","\\1",colnames(convered_matrix2))
clean_name <- t(apply(convered_matrix2, 1, function(x) tapply(x, colnames(convered_matrix2), mean)))


#read table patients_edit.txt
patients <- read.table("patients_edit.txt",h=T)
#Remove all samples for which we have no sample annotation in the file patients_edit.txt
annotation_samples <- clean_name[, which(colnames(clean_name) %in% c(patients$sample))]


ranker <-  t(apply(-annotation_samples, 1, rank))
devided_ranker <- t(apply(ranker, 1, function(x) x/ncol(ranker)))
q_ranker <- qnorm(devided_ranker)


#Finally, apply gene wise quantile normalization to transform the expression values to normal distribution
gene_quantail_normal <- function(df){
  ranker <-  t(apply(-annotation_samples, 1, rank))
  devided_ranker <- t(apply(ranker, 1, function(x) x/ncol(df)))
  qnormal <- qnorm(devided_ranker)
  
  return(qnormal)
}

qnormal_gene <- gene_quantail_normal(annotation_samples)

#plot genes
qnor_value <- qnormal_gene[1:50,]
qnor_value <- data.table(qnor_value, keep.rownames = T)
plotDat <- gather(qnor_value, key = "key", value = "value", -rn)
ggplot(plotDat, aes(rn, value)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "gene IDs")

#2.1.2 Gene annotations
gene_position <- fread("gene_position.tsv")
#remove gene at position X,Y, MT from file
gene_position <- gene_position[!which(gene_position$chromosome_name == "X" |gene_position$chromosome_name == "Y"| gene_position$chromosome_name == "MT"),]



#delete version number
rownames(qnormal_gene) <- gsub("^([^.]+).*","\\1",rownames(qnormal_gene))

#remove gene not in the file
cleaned_qnormal_gene <- qnormal_gene[which(rownames(qnormal_gene) %in% gene_position$ensembl_gene_id ),]
nrow(qnormal_gene)
nrow(cleaned_qnormal_gene) #23287 gene rest
ncol(cleaned_qnormal_gene) #445 samples

#2.2 Reformatting the SNP data
runPLINK <- function(PLINKoptions = "") system(paste("D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\plink\\plink", PLINKoptions))

samples <- colnames(cleaned_qnormal_gene) 
keep_names <- data.frame(samples,samples)
head(keep_names)
write.table(keep_names, file = "D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\keep_SNPdata.tsv", row.names=FALSE, sep="\t",col.names = F,quote = F)
runPLINK("plink --vcf chr22.vcf --geno 0.9 --maf 0.05 --hwe 1e-3 --keep keep_SNPdata.tsv --snps-only --recode A-transpose --out SNP-file" )
#79341 variantes pass filters and QC

SNPdata <- read.table("SNP-file.traw", h=T)


View(SNPdata)
#2.3 Covariates
#2.3.1 Measured Covariates
patients_V2 <- data.frame(patients)

nrow(patients_V2)
View(cleaned_qnormal_gene)
patients_V2 <- patients_V2[which(patients_V2$sample %in% samples),]
patients_V2$gender <- ifelse(patients_V2$gender == "male","1","0")
nrow(patients_V2)
View(patients_V2)
sum(patients_V2$gender == "1") #male 209
sum(patients_V2$gender == "0") #female 236
#because gender could have affect to the gene expression (mask or increase)

#2.3.2 Hidden factors (PEER)
peer_table <- read.csv("Estimated_Peer_Factors.csv")
View(peer_table)
#did any of the factors capture patterns inherent in the known covariates? 
relation_test <- merge(patients_V2, peer_table, by ="sample" )
View(relation_test)

long_relation_test <- melt(relation_test, id.vars = c("sample","pop","super_pop","gender"), variable.name = "peer")
ggplot(long_relation_test, aes(pop,value, color = pop))+
  geom_boxplot()+
  facet_wrap(~peer)

ggplot(long_relation_test, aes(super_pop,value, color = super_pop))+
  geom_boxplot()+
  facet_wrap(~peer)

ggplot(relation_test, aes(pop,PEER1)) + geom_boxplot()
#ggplot(relation_test, aes(pop,PEER2)) + geom_boxplot()
#ggplot(relation_test, aes(pop,PEER3)) + geom_boxplot()
#ggplot(relation_test, aes(pop,PEER4)) + geom_boxplot()
#ggplot(relation_test, aes(pop,PEER5)) + geom_boxplot()

#ggplot(relation_test, aes(super_pop,PEER1)) + geom_boxplot()
#ggplot(relation_test, aes(super_pop,PEER2)) + geom_boxplot()
#ggplot(relation_test, aes(super_pop,PEER3)) + geom_boxplot()
#ggplot(relation_test, aes(super_pop,PEER4)) + geom_boxplot()
#ggplot(relation_test, aes(super_pop,PEER5)) + geom_boxplot()

cor.test(example$gene2,example$snp2)
#t = 5.6493, df = 443, p-value = 2.886e-08       cor = 0.2592308

summary(lm(example$snp2~example$gene2 +example$gender))

example$snp2 <- as.character(example$snp2)
ggplot(example,aes(snp2,gene2,color = snp2))+geom_boxplot()+facet_wrap(~gender)



nrow(example)
example$sample <- rownames(example)
merged_relation <- merge(example, relation_test, by = c("sample","gender"))
View(merged_relation)

summary(lm(gene2 ~gender+snp2+PEER1+PEER2+PEER3+PEER4+PEER5+PEER6+PEER7+PEER8+PEER9+PEER10, data = merged_relation))
#2 boxplot
expression_value = relation_test[,4:ncol(relation_test)]
head(expression_value)
expression_value$gender = as.numeric(expression_value$gender)
lmexpression = lm(gender ~ PEER1+PEER2+PEER3+PEER4+PEER5+PEER6+PEER7+PEER8+PEER9+PEER10, data = expression_value)
points = rowSums(expression_value)
# Intercept 0.46966

#Reformat all covariates for matrix eQTL
matrix_eQTL_format <- t(relation_test)
View(matrix_eQTL_format)
colnames(matrix_eQTL_format) <- relation_test$sample
matrix_eQTL_format <- matrix_eQTL_format[4:nrow(matrix_eQTL_format),]


#3 eQTL analysis

#View(cleaned_qnormal_gene) #gene expression
#View(SNPdata) #SNP 
#View(matrix_eQTL_format) #covariates
#View(SNPdata[,c(2,7:ncol(SNPdata))])

#genotype
genotype <- SNPdata[,c(2,7:ncol(SNPdata))]
#write.table(genotype, file = "D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\genotype.tsv", sep="\t",row.names = F,col.names = T,quote = F)

#expression
genotype_name <- colnames(genotype) 
#head(genotype_name)

expression <- subset(cleaned_qnormal_gene, colnames(cleaned_qnormal_gene) %in% colnames(genotype))
#View(expression)
#write.table(express, file = "D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\input_exercise\\ms_expression.tsv", sep="\t",col.names = NA,quote = F)
#write.table(cleaned_qnormal_gene, file = "D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\expression.tsv", sep="\t",col.names = NA,quote = F)
test <- fread("D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\cleaned_qnormal_gene.tsv")
#View(test)

#covariate
covariate <- matrix_eQTL_format
#write.table(matrix_eQTL_format, file = "D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\covariate.tsv", sep="\t",col.names = NA,quote = F)

#gene location
genelocation <- gene_position[,c(1,2,3,4)]
#write.table(gene_position[,c(1,2,3,4)], file = "D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\genelocation.tsv", sep="\t",row.names = F,col.names = T,quote = F)

#SNP location
snplocation <- SNPdata[c(2,1,4)]
#write.table(SNPdata[c(2,1,4)], file = "D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\snplocation.tsv", sep="\t",row.names = F,col.names = T,quote = F)

#View(SNPdata[c(2,1,4,4)])

#View(SNPdata[,-c(1,3,4,5,6)])

express <- readRDS("expression.RDS")
#View(express)


base.dir = "D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\input_exercise"
# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-6;
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

useModel = modelLINEAR
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the space character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile("genotype.tsv");


gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the space character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile("ms_expression.tsv");


snpspos = read.table("SNP_location.tsv", header = T, stringsAsFactors = FALSE);
genepos = read.table("genelocation.tsv", header = T, stringsAsFactors = FALSE);
output_file_name_cis ="output_cis.tsv";
output_file_name_tra = "output_trans.tsv";


cvrt= SlicedData$new();
cvrt$fileDelimiter = "\t";      # the space character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile("covariate.tsv");

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);

qq_me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);


plot(me)
plot(qq_me)
me$cis$ntests
#2570223 tests

#bonferroni corrected a = 0.05/2570223
bonferroni_threshold = 0.05/2570223
me$param

output_cis <- read.table("output_cis.tsv",h = T)
sum(output_cis$p.value < bonferroni_threshold )
# 4444 after FWER<0.05
# FWER is calculated as 1-(1- a)^number_of_tests
# If the threshold is not corrected = 0.05, with 6 tests, there is a 26.5% chance of 
# discovering one or more false-positive results.
# With Bonferroni corrected threshold (a/number_of_tests) 
# 1-(1-0.008)^6 = 4.7%
# For bonferroni correction keeps FWER at 5% 
# Bonferroni correction adjust the probability value to avoid increased risk
# of Type I error when making multiple statistical tests.
# 

#3.1
#chr1
ch1_snp_data <- fread("Ch1-SNP-file.traw",h=T)
head(ch1_snp_data[, c(2,7:451)])
#genotype
ch1_genotype <- ch1_snp_data[,c(2,7:451)]
ch1_genotype
ch1_genotype
write.table(ch1_genotype, file = "ch1_genotype.tsv", sep="\t",row.names = F,col.names = T,quote = F)
#SNP location

ch1_snplocation <- ch1_snp_data[,c(2,1,4)]
write.table(ch1_snplocation, file = "ch1_snplocation.tsv", sep="\t",row.names = F,col.names = T,quote = F)

#matrix eqtl

base.dir = "D:\\Study thingie\\FACH\\S5 WS2021\\GOBI\\eQTL\\input_exercise"
# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-6;
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

useModel = modelLINEAR
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the space character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile("ch1_genotype.tsv");


gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the space character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile("ms_expression.tsv");


snpspos = read.table("ch1_snplocation.tsv", header = T, stringsAsFactors = FALSE);
genepos = read.table("genelocation.tsv", header = T, stringsAsFactors = FALSE);
output_file_name_cis ="ch1_output_cis.tsv";
output_file_name_tra = "ch1_output_trans.tsv";


cvrt= SlicedData$new();
cvrt$fileDelimiter = "\t";      # the space character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile("covariate.tsv");

ch1_me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);

ch1_qq_me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);

plot(ch1_me)
plot(ch1_qq_me,  pch = 16, cex = 0.7)


ch1_me$cis$ntests
#bonferroni corrected a = 0.05/2570223

ch1_bonferroni_threshold = 0.05/ch1_me$cis$ntests


ch1_output_cis <- read.table("ch1_output_cis.tsv",h = T)
sum(ch1_output_cis$p.value < ch1_bonferroni_threshold )
#1269

#ch6
ch6_snp_data <- fread("Ch6-SNP-file.traw",h=T)
ncol(ch6_snp_data)

#genotype
ch6_genotype <- ch6_snp_data[,c(2,7:451)]

write.table(ch6_genotype, file = "ch6_genotype.tsv", sep="\t",row.names = F,col.names = T,quote = F)
#SNP location

ch6_snplocation <- ch6_snp_data[,c(2,1,4)]
write.table(ch6_snplocation, file = "ch6_snplocation.tsv", sep="\t",row.names = F,col.names = T,quote = F)

snps$LoadFile("ch6_genotype.tsv");
snpspos = read.table("ch6_snplocation.tsv", header = T, stringsAsFactors = FALSE)
output_file_name_cis ="ch6_output_cis.tsv";
output_file_name_tra = "ch6_output_trans.tsv"

ch6_me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);

ch6_qq_me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);

plot(ch6_me)
plot(ch6_qq_me,  pch = 16, cex = 0.7)

ch6_me$cis$ntests
#bonferroni corrected a = 0.05/7128130

ch6_me$cis$neqtls
#53267

ch6_bonferroni_threshold = 0.05/ch6_me$cis$ntests


ch6_output_cis <- read.table("ch6_output_cis.tsv",h = T)
sum(ch6_output_cis$p.value < ch6_bonferroni_threshold )
#34246

#ch11
ch11_snp_data <- fread("Ch11-SNP-file.traw",h=T)
ncol(ch6_snp_data)

#genotype
ch11_genotype <- ch11_snp_data[,c(2,7:451)]

write.table(ch11_genotype, file = "ch11_genotype.tsv", sep="\t",row.names = F,col.names = T,quote = F)
#SNP location

ch11_snplocation <- ch11_snp_data[,c(2,1,4)]
write.table(ch11_snplocation, file = "ch11_snplocation.tsv", sep="\t",row.names = F,col.names = T,quote = F)

snps$LoadFile("ch11_genotype.tsv");
snpspos = read.table("ch11_snplocation.tsv", header = T, stringsAsFactors = FALSE)
output_file_name_cis ="ch11_output_cis.tsv";
output_file_name_tra = "ch11_output_trans.tsv"

ch11_me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);

ch11_qq_me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = useModel,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);

plot(ch11_me)
plot(ch11_qq_me,  pch = 16, cex = 0.7)

ch11_me$cis$ntests
#bonferroni corrected a = 0.05/4880409
ch11_me$cis$neqtls
#11962
ch11_bonferroni_threshold = 0.05/ch11_me$cis$ntests


ch11_output_cis <- read.table("ch11_output_cis.tsv",h = T)
sum(ch11_output_cis$p.value < ch11_bonferroni_threshold )
#6483

#4.1 Genomic distance of cis eQTLs
snpspos = read.table("SNP_location.tsv", header = T, stringsAsFactors = FALSE)
genepos = read.table("genelocation.tsv", header = T, stringsAsFactors = FALSE)


gene_loc_22 <- subset(genepos, chromosome_name == 22)
ch22_tss_loc <- as.numeric(unlist(gene_loc_22[3]))

ch22_tss_loc_sort <- ch22_tss_loc[order(ch22_tss_loc)]

ch22_snp_loc <- as.numeric(snpspos$POS)

ch22_snp_loc_sort <- ch22_snp_loc[order(ch22_snp_loc)]
ch22_tss_loc_sort- ch22_snp_loc_sort[1:613]

hist(ch22_tss_loc_sort- ch22_snp_loc_sort)

#build data frams
ch22_cis_eqtl <- read.table("output_cis.tsv",h=T)
ch22_max_beta <- ch22_cis_eqtl %>% group_by(SNP) %>% top_n(1,abs(beta))

ch22_ceqtl_dis <- data.frame(ch22_cis_eqtl)


ch22_snp_loc <-as.numeric(gsub("^[^-]*_([^_]+).*", "\\1", ch22_ceqtl_dis$SNP))
ch22_ceqtl_dis$snp_loc <- ch22_snp_loc

ch22_result <- cbind(ch22_ceqtl_dis, genepos[match(ch22_ceqtl_dis$gene,genepos$ensembl_gene_id),3])
ch22_result$distance <-as.numeric(abs(ch22_result$`genepos[match(ch22_ceqtl_dis$gene, genepos$ensembl_gene_id), `- ch22_result$snp_loc ))


View(ch22_result)
split(ch22_cis_eqtl, ch22_cis_eqtl$gene)

ch22_nearest <- ch22_result %>% group_by(SNP) %>%top_n(-1,distance)

ch22_nearest <- as.data.frame.table(ch22_nearest)
#match nearest and most accosiated together

sum((ch22_max_beta$SNP == ch22_nearest$Freq.SNP)&(ch22_max_beta$gene == ch22_nearest$Freq.gene))
nrow(ch22_max_beta)

ch22_nearest$snp_gene <- paste(ch22_nearest$Freq.SNP,ch22_nearest$Freq.gene,sep = "-")
ch22_max_beta$snp_gene <- paste(ch22_max_beta$SNP, ch22_max_beta$gene, sep = "-")


sum(ch22_max_beta$snp_gene %in% ch22_nearest$snp_gene) #4629/5286 = 0.88 
sum(ch22_nearest$snp_gene %in% ch22_max_beta$snp_gene) #37032/42288 = 0.88

#4.2.1 Quantify overlap
gwas_catalog <- read.table("gwas_catalog_v1_0_filtered.tsv", sep = "\t", h = T) 

catalog_duplicate <- gwas_catalog[!duplicated(gwas_catalog$SNPS),]
sum(length(catalog_duplicate$SNP_GENE_IDS)> 1)
nrow(catalog_duplicate)

w = 0
for(x in catalog_duplicate$SNP_GENE_IDS){
  if(nchar(x) > 15){
    w = w+1
  }
}
w # 14 

nrow(catalog_duplicate) #112
nrow(gwas_catalog) #154
#154-112 = 42  

catalog_duplicate <- gwas_catalog[!duplicated(gwas_catalog$UPSTREAM_GENE_ID),]
catalog_duplicate <- gwas_catalog[!duplicated(gwas_catalog$DOWNSTREAM_GENE_ID),]
catalog_duplicate <- gwas_catalog[!duplicated(gwas_catalog$SNP_GENE_IDS),]
nrow(catalog_duplicate) #38 left

#4.2.2

x <- catalog_duplicate[grep("snp_11_838722",catalog_duplicate$DISEASE.TRAIT),]
x$DISEASE.TRAIT
gene_position <- fread("gene_position.tsv")
View(gene_position)
subset(gene_position, ensembl_gene_id == c("ENSG00000280759","ENSG00000183020","ENSG00000281385"))
y <- gene_position[grep("ENSG00000280759"|"ENSG00000183020"|"ENSG00000281385",gene_position$ensembl_gene_id),]

gene_position[grepl("ENSG00000280759|ENSG00000183020|ENSG00000281385",gene_position$ensembl_gene_id)]
