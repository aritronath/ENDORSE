#####################################
########### GSE76040 ################
library(ggplot2)
library(ggpubr)

#Expression 
gse76040 <- read.delim("./../Data/GSE76040/GSE76040_series_matrix.txt", comment.char="!", header=TRUE, fill=TRUE, row.names=1)

#Header annotation 
gse76040.annot <- read.delim("./../Data/GSE76040/GSE76040_series_matrix copy.txt", header=FALSE, fill=TRUE, nrows=7)
gse76040.annot <- data.frame(t(gse76040.annot[,-1]))

temp <- apply(gse76040.annot, 2, function (x) strsplit(as.character(x), ":"))
gse76040.annotV2 <- data.frame(matrix(data=NA, nrow=nrow(gse76040.annot), ncol=5))
for (i in 3:7) {
  gse76040.annotV2[, i-2] <- sapply(temp[[i]], '[', 2)
  colnames(gse76040.annotV2)[i-2] <- temp[[i]][[1]][1]
}
rownames(gse76040.annotV2) <- gse76040.annot$X2

#Array information 
#gse76040.genes <- read.delim("./GSE76040/GPL6480-9577.txt", comment.char="#", header=TRUE)
library(biomaRt)
gids <- strtrim(rownames(gse76040), 13)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
bm.gids <- getBM(attributes = c("affy_hta_2_0", "hgnc_symbol"),
                 filters    = "affy_hta_2_0",
                 values     = gids, 
                 mart       = mart, 
                 verbose  = FALSE)

x <- match(gids, bm.gids$affy_hta_2_0)
gse76040.mat <- gse76040[which(!is.na(x)), ]
agg.gse76040 <- aggregate(gse76040.mat, by=list(bm.gids$hgnc_symbol[na.omit(x)]), FUN=max)

gse76040.mat <- as.matrix(agg.gse76040[-1, 2:ncol(agg.gse76040)])
rownames(gse76040.mat) <- agg.gse76040$Group.1[-1]

gse76040.s <- scale(t(gse76040.mat))

#### ER Response ssGSEA ---------------------
ENDORSE <- readRDS("./../Output/ENDORSE.RDS")

temp3 <- scale(t(metabric.exp))
MBsc.ssgsea <- gsva(t(temp3), list("ENDORSE"=ENDORSE),
                    verbose=TRUE, kcdf="Gaussian", method='ssgsea', ssgsea.norm=T)
MBsc.ssgsea[1,] <- ReScale(MBsc.ssgsea[1,])

gse76040.ssgsea <- gsva(t(gse76040.s), list("ENDORSE"=ENDORSE),
                        verbose=TRUE, kcdf="Gaussian", method='ssgsea', ssgsea.norm=T)
gse76040.ssgsea[1,] <- ReScale(gse76040.ssgsea[1,])

# Train df 
df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 'ENDORSE'=MBsc.ssgsea[1,])

cfit = coxph(Surv(Time, Event) ~ ENDORSE, data=df)

# predict
gse76040.pred <- predict(cfit, data.frame("ENDORSE"=gse76040.ssgsea[1,]), type="risk")

# compare
# plot - by risk groups
gse76040.pred.grp <- MakeRiskGrp(gse76040.pred, thresh = c(1,2))

df <- data.frame("Ki67"=as.numeric(gse76040.annotV2$`ki67 (%)`), 
                 "RiskGroup"=gse76040.pred.grp, 
                 "Risk"=gse76040.pred)

summary(aov(Ki67 ~ RiskGroup, data=df))
P1 <- ggviolin(df, x="RiskGroup", y="Ki67", fill="RiskGroup", palette = "nejm", 
               add="boxplot", add.params = list(fill = "white", width=0.1)) +
  labs(x="", y="Ki67 (%)", title="CONFIRM: ER+ fulvestrant") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

# plot - by continous predictions
summary(lm(as.numeric(gse76040.annotV2$`ki67 (%)`) ~ gse76040.pred))
cor.test(as.numeric(gse76040.annotV2$`ki67 (%)`), gse76040.pred)

P2 <- ggscatter(df, x="Risk", y="Ki67", color = "salmon",
                add="reg.line", add.params = list(color = "steelblue", width=0.5), conf.int = TRUE) + 
  labs(x="Risk", y="Ki67 (%)") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

pdf("./../Figures/Predictions_GSE76040.pdf", width=11, height=5)
egg::ggarrange(P1, P2, ncol=2)
dev.off()
