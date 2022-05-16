#######################
###### GSE130788 ######

### Load the data --------------------
gse130788 <- read.delim("./../Data/GSE130788/GSE130788_series_matrix.txt", comment.char = "!",
                        blank.lines.skip = T, header=T, fill=T)
gse130788.genes <- read.delim("./../Data/GSE130788/GPL6480-9577.txt", comment.char = "#",
                              blank.lines.skip = T, header=T, fill=T)

gse130788.annot <- t(read.delim("./../Data/GSE130788/GSE130788_series_matrix.txt", 
                              skip = 45, nrows = 5, 
                              blank.lines.skip = T, header=F, fill=T))

### Pre-process and prep data ----------
# Expression 
gse130788.mat <- as.matrix(gse130788[, 2:ncol(gse130788)])
x <- match(gse130788[,1], gse130788.genes$ID)
gids <- gse130788.genes$GENE_SYMBOL[x]
agg.gse130788 <- aggregate(gse130788.mat, by=list(gids), FUN=mean)

library(impute)
temp <- as.matrix(agg.gse130788[,2:ncol(agg.gse130788)])
temp2 <- impute.knn(temp)
gse130788.s <- scale(t(temp2$data))
colnames(gse130788.s) <- agg.gse130788$Group.1

# Phenotype
gse130788.annot <- data.frame(gse130788.annot[-1,])
colnames(gse130788.annot) <- c("Sample_time", "ER", "PR", "Group", "Response")
 
### Calculate ENDORSE risk scores -----------
ENDORSE <- readRDS("./../Output/ENDORSE.RDS")

gse130788.ssgsea <- gsva(t(gse130788.s), list("ENDORSE"=ENDORSE, "E2"=names(E2_geneset)),
                        verbose=TRUE, kcdf="Gaussian", method='ssgsea')
gse130788.ssgsea[1,] <- ReScale(gse130788.ssgsea[1,])

# Train df 
df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 'ENDORSE'=MBsc.ssgsea[1,])

cfit = coxph(Surv(Time, Event) ~ ENDORSE, data=df)

# predict
gse130788.pred <- predict(cfit, data.frame("ENDORSE"=gse130788.ssgsea[1,]), type="risk")
gse130788.pred.grp <- MakeRiskGrp(gse130788.pred, thresh = c(1,2))

df <- data.frame(gse130788.annot,
                 "RiskGroup"=gse130788.pred.grp, 
                 "Risk"=gse130788.pred)
Response_Lab <- array(dim=nrow(df))
Response_Lab[grep("PCR", df$Response)] <- "PCR"
Response_Lab[grep("RD", df$Response)] <- "RD"
df$Response_Lab <- as.factor(Response_Lab)

BL <- grep('baseline', gse130788.annot$Sample_time)
EP <- grep('Pos', gse130788.annot$ER)
EN <- grep('Neg', gse130788.annot$ER)

SetP <- intersect(BL, EP)
SetN <- intersect(BL, EN)

summary(aov(Risk ~ Response_Lab, data=df[SetP,])) #P=0.465
summary(aov(Risk ~ Response_Lab, data=df[SetN,])) #P=0.824

P1 <- ggviolin(df[SetP,], x="Response_Lab", y="Risk", fill="Response_Lab", palette = "Purples", 
               add="boxplot", add.params = list(fill = "white", width=0.1)) +
  labs(x="Response", y="ENDORSE risk", title="Neoadvuvant trastuzumab and/or lapatinib", 
       subtitle="ER+") +
  theme_classic(base_size = 18) 

P2 <- ggviolin(df[SetN,], x="Response_Lab", y="Risk", fill="Response_Lab", palette = "Purples", 
               add="boxplot", add.params = list(fill = "white", width=0.1)) +
  labs(x="Response", y="ENDORSE risk", subtitle="ER-") +
  theme_classic(base_size = 18) 
  
#
pdf("./../Figures/Predictions_GSE130788.pdf", width=11, height=5)
egg::ggarrange(ggpar(P1, legend='none'), ggpar(P2, legend='none'), ncol=2)
dev.off()




