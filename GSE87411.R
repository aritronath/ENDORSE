#####################################
########### GSE87411 ################

#Expression 
gse87411 <- read.delim("./../Data/GSE87411/GSE87411_V2.txt", comment.char="!", header=TRUE, fill=TRUE)

#Header annotation 
gse87411.annot <- read.delim("./../Data/GSE87411/GSE87411_V2.txt", header=FALSE, fill=TRUE, nrows=6)
gse87411.annot <- t(gse87411.annot)
colnames(gse87411.annot) <- c("Patient", "BL_KI", "EOT_KI", "ER_Response", "Treatment", "ID_REF")
gse87411.annot <- gse87411.annot[-1,]

gse87411.annot[grep("BL", gse87411.annot[,1], ignore.case=T), 1] <- "BL"
gse87411.annot[grep("W2", gse87411.annot[,1], ignore.case=T), 1] <- "W2"
gse87411.annot[grep("sensitive", gse87411.annot[,4], ignore.case=T), 4] <- "sensitive"
gse87411.annot[grep("resistant", gse87411.annot[,4], ignore.case=T), 4] <- "resistant"

#Array information 
gse87411.genes <- read.delim("./../Data/GSE87411/GPL6480-9577.txt", comment.char="#", header=TRUE)

#Pre-process expression matrix
x <- match(gse87411$ID_REF, gse87411.genes$ID)
agg.gse87411 <- aggregate(gse87411, by=list(gse87411.genes$GENE_SYMBOL[x]), FUN=mean)
gse87411.mat <- as.matrix(agg.gse87411[-1, 3:ncol(agg.gse87411)])
rownames(gse87411.mat) <- agg.gse87411$Group.1[-1]

library(impute)
gse87411.mat <- impute.knn(gse87411.mat, k=50)
gvars <- apply(gse87411.mat[[1]], 1 , var)

gse87411.s <- scale(t(gse87411.mat[[1]]))

##### ENDORSE signatures -----------------------------
ENDORSE <- readRDS("./../Output/ENDORSE.RDS")

gse87411.ssgsea <- gsva(t(gse87411.exp.bl), list("ENDORSE"=ENDORSE),
                            verbose=TRUE, kcdf="Gaussian", method='ssgsea')
gse87411.ssgsea[1,] <- ReScale(gse87411.ssgsea[1,])

# Train df 
df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 'ENDORSE'=MBsc.ssgsea[1,])

cfit = coxph(Surv(Time, Event) ~ ENDORSE, data=df)

# predict
gse87411.pred <- predict(cfit, data.frame("ENDORSE"=gse87411.ssgsea[1,]), type="risk")

# compare
# plot - by risk groups
gse87411.pred.grp <- MakeRiskGrp(gse87411.pred, thresh = c(1,2))

df <- data.frame("Ki67_BL"=ki67_BL,
                 "Ki67_EOT"=ki67_EOT,
                 "Ki67_delta"=ki67_BL - ki67_EOT,
                 "RiskGroup"=gse87411.pred.grp, 
                 "Risk"=gse87411.pred)

summary(aov(ki67_BL ~ RiskGroup, data=df))
summary(aov(ki67_EOT ~ RiskGroup, data=df))

P1 <- ggviolin(df, x="RiskGroup", y="Ki67_BL", fill="RiskGroup", palette = "nejm", 
               add="boxplot", add.params = list(fill = "white", width=0.1)) +
  labs(x="", y="Ki67 (%)", title="NCT00265759: ER+ neoadjuvant aromatase inhbitor") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

P2 <- ggviolin(df, x="RiskGroup", y="Ki67_EOT", fill="RiskGroup", palette = "nejm", 
               add="boxplot", add.params = list(fill = "white", width=0.1)) +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  labs(x="", y="Ki67 (%)") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

# plot - by continous predictions
summary(lm(ki67_BL ~ gse87411.pred))
cor.test(ki67_BL, gse87411.pred)

summary(lm(ki67_EOT ~ gse87411.pred))
cor.test(ki67_EOT, gse87411.pred)

P3 <- ggscatter(df, x="Risk", y="Ki67_BL", color = "salmon",
                add="reg.line", add.params = list(color = "steelblue", width=0.5), conf.int = TRUE) + 
  labs(x="Risk", y="Ki67 (%)") +
  theme_classic(base_size = 18)

P4 <- ggscatter(df, x="Risk", y="Ki67_EOT", color = "salmon",
                add="reg.line", add.params = list(color = "steelblue", width=0.5), conf.int = TRUE) + 
  labs(x="Risk", y="Ki67 (%)") +
  theme_classic(base_size = 18)

pdf("./../Figures/Predictions_GSE87411.pdf", width=10, height=8)
egg::ggarrange(P1, P3, P2, P4, ncol=2)
dev.off()





