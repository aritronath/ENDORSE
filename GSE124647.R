#####################################
########### GSE124647 ################

#Expression 
gse124647 <- read.delim("./../Data/GSE124647/GSE124647_series_matrix.txt", comment.char="!", header=TRUE, fill=TRUE, row.names=1)

#Header annotation 
gse124647.annot <- read.delim("./../Data/GSE124647/GSE124647_series_matrix copy.txt", header=FALSE, fill=TRUE, skip=10, nrows=16)
gse124647.annot <- data.frame(t(gse124647.annot[,-1]))
rownames(gse124647.annot) <- gse124647.annot$X16

temp <- apply(gse124647.annot, 2, function (x) strsplit(as.character(x), ":"))
gse124647.annotV2 <- data.frame(matrix(data=NA, nrow=nrow(gse124647.annot), ncol=15))
for (i in 1:15) {
  gse124647.annotV2[,i] <- sapply(temp[[i]], '[', 2)
  colnames(gse124647.annotV2)[i] <- temp[[i]][[1]][1]
}
rownames(gse124647.annotV2) <- gse124647.annot$X16

#Array information 
gse124647.genes <- read.delim("./../Data/GSE124647/GPL96-57554.txt", comment.char="#", header=TRUE)

#Pre-process expression matrix
x <- match(rownames(gse124647), gse124647.genes$ID)
agg.gse124647 <- aggregate(gse124647, by=list(gse124647.genes$Gene.Symbol[x]), FUN=mean)
gse124647.mat <- agg.gse124647[-1, 2:ncol(agg.gse124647)]
rownames(gse124647.mat) <- agg.gse124647$Group.1[-1]

gvars <- apply(gse124647.mat, 1 , var)
gse124647.s <- scale(t(gse124647.mat[-which(gvars < 0.1), ]))


### ER response and  survival analysis -------
library(GSEABase)
library(GSVA)
ENDORSE <- readRDS("./../Output/ENDORSE.RDS")

gse124647.ssgsea <- gsva(t(gse124647.s), list("ENDORSE"=ENDORSE),
                         verbose=TRUE, kcdf="Gaussian", method='ssgsea', ssgsea.norm=T)
gse124647.ssgsea[1,] <- ReScale(gse124647.ssgsea[1,])

gse124647.os <- Surv(time=as.numeric(gse124647.annotV2$`oas months`), event=as.numeric(gse124647.annotV2$`os event`))

# Predict risk groups --- 
# Train df 
df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 'ENDORSE'=MBsc.ssgsea[1,])

cfit = coxph(Surv(Time, Event) ~ ENDORSE, data=df)

# predict
gse124647.pred <- predict(cfit, data.frame("ENDORSE"=gse124647.ssgsea[1,]), type="risk")

# compare
# plot - by risk groups
gse124647.pred.grp <- MakeRiskGrp(gse124647.pred, thresh=c(1,2))

df <- data.frame("OSTime"=as.numeric(gse124647.annotV2$`oas months`),
                 "OSEvent"=as.numeric(gse124647.annotV2$`os event`), 
                 "RiskGroup"=gse124647.pred.grp, 
                 "Risk"=gse124647.pred, 
                 "ENDORSE"=gse124647.ssgsea[1,])

# Survival plots ----
library(survminer)
fit.os <- survfit(Surv(OSTime, OSEvent) ~ RiskGroup, data=df)

P1 <- ggsurvplot(fit.os, data=df, palette="nejm", 
                 censor.shape="|", censor.size=2, 
                 break.time.by=12, 
                 legend.title="ENDORSE", legend.labs=c("Low risk", "Medium risk",  "High risk"),
                 surv.median.line = "hv",
                 risk.table=T, 
                 title="ER+/HER2- MBC",
                 ggtheme=theme_classic2(base_size = 18))

summary(coxph(Surv(OSTime, OSEvent) ~ RiskGroup, data=df))

pdf("./../Figures/Prediction_GSE124647_os.pdf", width=12, height=8)
print(P1)
dev.off()


