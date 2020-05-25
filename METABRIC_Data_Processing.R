library(data.table)
library(ggplot2)
library(impute)

#### Load and process METABRIC expression and related clinical data ----------------------------------------------- 
metabric <- fread("./../Data/METABRIC/data_expression_medianV2.txt", header=TRUE)
glist <- metabric$Hugo_Symbol
metabric <- as.matrix(metabric[,3:ncol(metabric)])
rownames(metabric) <- glist
metabric <- impute.knn(metabric, k = 50)$data

# load METABRIC mutations and CNV
metabric.mut <- fread("./../Data/METABRIC/data_mutations_mskcc.txt", header=T)
metabric.cnv <- fread("./../Data/METABRIC/data_CNA.txt", header=T)

y <- which(metabric.mut$Consequence == "synonymous_variant") # remove synonmous mutations
metabric.mut <- metabric.mut[-y,]

# load METABRIC clinical data
metabric.annot <- read.delim("./../Data/METABRIC/data_clinical_patient.txt", comment.char='#', fill=T, header=T)
MB.annot.extra <- read.delim("./../Data/METABRIC/brca_metabric_clinical_data.tsv", header=T) 

# Select ER+/HER2- cancers that did not receive chemotherapy and are either alive or died due to the disease
F.INDEX <- which(metabric.annot$ER_IHC == 'Positve' & 
                   metabric.annot$CHEMOTHERAPY != 'YES' & 
                   MB.annot.extra$HER2.Status == 'Negative' &
                   (metabric.annot$VITAL_STATUS == 'Living' | metabric.annot$VITAL_STATUS == 'Died of Disease')) # filters down from 2509 to 833 


### Run GSVA in ssGSEA mode with hallmark and C2 gene sets and output to file -------------------------------------------
library(GSEABase)
library(GSVA)

c2set <- GSEABase::getGmt("./../Data/GeneSets/c2.all.v7.0.symbols.gmt")
Hset <- GSEABase::getGmt("./../Data/GeneSets/H.all.v7.0.symbols.gmt")

gsvaData <- gsva(metabric, Hset, min.sz=10, verbose=TRUE, kcdf="Poisson", method="ssgsea")
write.table(gsvaData, file="./../Output/METABRIC_H.ssGSEA_scores.tsv", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
METABRIC.H <- gsvaData

c(40:60)[which(ncol(metabric) %% c(40:60) == 0)] #find factor 
ColGrps <- split(1:ncol(metabric), 1:56)
gsvaData <- gsva(metabric[, ColGrps[[1]]], c2set, min.sz=10, verbose=TRUE, kcdf="Poisson", method="ssgsea")
t1 <- proc.time()
for (i in 2:length(ColGrps)) {
  TMP <- gsva(metabric[, ColGrps[[i]]], c2set, min.sz=10, verbose=TRUE, kcdf="Poisson", method="ssgsea")
  gsvaData <- cbind(gsvaData, TMP)
  print(paste("\n Done", i, "out of", length(ColGrps), sep=" "))
}
t2 <- proc.time() - t1
write.table(gsvaData, file="./../Output/METABRIC_C2.ssGSEA_scores.tsv", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
METABRIC.C2 <- gsvaData 

### Process annotations save objects for further analysus ------------------------------------------------
metabric.annotV2 <- na.omit(metabric.annot[F.INDEX, c(1, 12, 13, 14)])
x <- match(metabric.annotV2$PATIENT_ID, MB.annot.extra$Patient.ID)
MB.annot.extraV2 <- MB.annot.extra[x,]

metabric.hall <- METABRIC.H[, na.omit(match(metabric.annotV2$PATIENT_ID, colnames(METABRIC.H)))]
metabric.c2 <- METABRIC.C2[, na.omit(match(metabric.annotV2$PATIENT_ID, colnames(METABRIC.C2)))]
metabric.exp <- metabric[, na.omit(match(metabric.annotV2$PATIENT_ID, colnames(metabric)))]

save(metabric.annotV2, MB.annot.extraV2, metabric.exp, file="./../Output/metabric_er.RData")
save(metabric.c2, metabric.hall, file="./../Output/metabric.ssGSEA.RData")



