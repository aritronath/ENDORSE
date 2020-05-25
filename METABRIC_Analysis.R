library(ggplot2)
library(egg)
library(ggpubr)
library(survminer)
library(ggsci)

### Load results from 50 repeat 10-fold cross-validation analyses ------------

load("./../Output/metabric_er.RData")
load("./../Output/metabric.ssGSEA.RData")

FS_C2 <- readRDS("./../Output/MB_C2_Model_cv_results_FeatureSet.RDS")
FS_hall <- readRDS("./../Output/MB_H_Model_cv_results_FeatureSet.RDS")
FS_sg <- readRDS("./../Output/MB_Cox_cv_results_FeatureSet.RDS")
FS_ENDORSE <- readRDS("./../Output/MB_CoxNoClin_cv_results_FeatureSet.RDS")
FS_ENDORSE2 <- readRDS("./../Output/MB_CoxNoClin_cv_2_results_FeatureSet.RDS")

Res_C2 <- readRDS("./../Output/MB_C2_Model_cv_results_Restable.RDS")
Res_hall <- readRDS("./../Output/MB_H_Model_cv_results_Restable.RDS")
Res_sg <- readRDS("./../Output/MB_Cox_cv_results_Restable.RDS")
Res_ENDORSE <- readRDS("./../Output/MB_CoxNoClin_cv_results_Restable.RDS")
Res_ENDORSE2 <- readRDS("./../Output/MB_CoxNoClin_cv_2_results_Restable.RDS")

FS_sg <- FS_sg[-grep("NA/NaN", FS_sg)]


#### Plot the general characterisitics of the models ------------------------------------
# 1. Features in each analysis (table barplot) 
FS <- list(FS_C2, FS_hall, FS_sg, FS_ENDORSE)
Plts <- list(); Plts2 <- list()

for (i in 1:length(FS)) {  
  df <- data.frame(sort(table(FS[[i]]), decreasing = T))
  colnames(df) <- c("Feature", "Frequency")
  
  Plts2[[i]] <-ggplot(data=df, aes(x=Frequency)) + 
    geom_histogram(color="white", fill="#000370", binwidth = 25) +
    geom_vline(xintercept = 250) +
    labs(x="Frequency", y="Count") +
    theme_classic()
  
  df <- df[which(df$Frequency > 250), ]
  print(sum(df$Frequency == 500))
  #sz <- max(10 - nrow(df)/50, 3)
  Plts[[i]] <- ggplot(data=df) + 
    geom_bar(aes(x=Feature, y=Frequency), width=0.5, color=NA, fill='#FFC67E', stat='identity') +
    labs(x="") +
    theme_classic() +
    coord_flip()

#  write.csv(df, file=paste("./../Output/Feature_Set_Table_", i, ".csv", sep=""))
}

png(file="./../Figures/Feature_bar.png", height=3000, width=3600, res=200) 
egg::ggarrange(Plts[[1]], Plts[[2]], Plts[[3]], Plts[[4]], ncol=2)
dev.off()

png(file="./../Figures/Feature_hist.png", height=1600, width=3600, res=200) 
egg::ggarrange(Plts2[[1]], Plts2[[2]], Plts2[[3]], Plts2[[4]], ncol=2)
dev.off()


# 2. Consensus geneset application - METABRIC stratification with risk score
ENDORSE <- names(which(table(FS_ENDORSE) > 250)) # Define ENDORSE as genes in atleast 50% of the cross-validations (>250) 
saveRDS(ENDORSE, file="./../Output/ENDORSE.RDS")

# Calculate ENDORSE enrichment score 
MB.ENDORSE <- gsva(metabric.exp, gset.idx.list = list("ENDORSE"=ENDORSE), 
                   method="ssgsea", kcdf="Gaussian")

MakeRiskGrp <- function (risk.pred) {
  risk.grps <- array(length(risk.pred))
  risk.grps[which(risk.pred <= 1)] <- "Low Risk"
  risk.grps[which(risk.pred > 1 & risk.pred < 2)] <- "Medium Risk"
  risk.grps[which(risk.pred >= 2)] <- "High Risk"
  risk.grps <- factor(risk.grps, levels=c("Low Risk", "Medium Risk", "High Risk"))
  return(risk.grps)
}

df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 'ENDORSE'=MB.ENDORSE[1,])

cfit = coxph(Surv(Time, Event) ~ Age + ENDORSE, data=df)
saveRDS(cfit, file="./../Output/ENDORSE.CoxFit.RDS") 

res <- summary(cfit)

# Three steps to predict and classify sample: 
# 1. Calculate ENDORSE enrichment score (gsva)
# 2. Apply ENDORSE to predict risk (predict(coxfit)
# 3. Classify risk grp (stratify)

MB.risk <- predict(cfit, type='risk')
MB.risk.grps <- MakeRiskGrp(MB.risk)

MB.Surv.df <- data.frame("Time"=metabric.annotV2$OS_MONTHS, 
                   "Event"=as.numeric(metabric.annotV2$OS_EVENT), 
                   "Age"=metabric.annotV2$AGE_AT_DIAGNOSIS, 
                   "Group"=MB.risk.grps)

fit <- survfit(Surv(Time, Event) ~ Group, data=MB.Surv.df)
Plt <- ggsurvplot(fit, data=MB.Surv.df, palette="nejm", 
                  censor.shape="|", censor.size=2, 
                  break.time.by=60, 
                  legend.title="ENDORSE", legend.labs=c("Low risk", "Medium risk",  "High risk"),
                  surv.median.line = "hv",
                  risk.table=T, 
                  ggtheme=theme_classic2(base_size = 18))

ggexport(Plt, filename ="./../Figures/METABRIC.Cox.Survival.pdf", width=8, height=9)
dev.off()

table(MB.risk.grps)
hist(MB.risk, freq=T, breaks=10)
lines(ecdf(MB.risk), lwd=2, col="skyblue")
abline(v=c(1,2), lwd=2, lty=2, col="darkred")

100 * sum(MB.risk >= 2) / length(MB.risk)
100 * sum(MB.risk <= 1) / length(MB.risk)

# Make histogram and ECDF of risk distribution 
Plt[[1]] <- ggplot(data.frame(MB.risk), aes(MB.risk)) + 
  geom_histogram(color=NA, fill='#2D004E', bins=50) +  
  geom_rug(color="#8D8D8D") +
  labs(x="ENDORSE risk", y="Frequency") + 
  #geom_vline(xintercept = c(1,2), color = "grey50", size=0.75) +
  scale_x_continuous(breaks=c(1, 2, 3, 4)) +
  theme_classic(base_size = 18)

Plt[[2]] <- ggplot(data.frame(MB.risk), aes(MB.risk)) + 
  stat_ecdf(geom = "step", size=1.5, color="#D10000") +
  labs(y = "F(risk)", x="ENDORSE risk") +
  geom_hline(yintercept = c(0.4622, 0.8367), linetype = "dotted", color="#858585") +
  #geom_segment(x=0, xend=2, y=0.8367, yend = 0.8367, size=0.5, linetype = "dashed", col="#8D8D8D") +
  #geom_segment(x=0, xend=1, y=0.4622, yend =  0.4622, size=0.5, linetype = "dashed", col="#8D8D8D") +
  geom_segment(x=2, xend=2, y=0, yend = 0.8367,linetype = "dotted", col="#858585") +
  geom_segment(x=1, xend=1, y=0, yend =  0.4622, linetype = "dotted", col="#858585") +
  scale_x_continuous(breaks=c(1, 2, 3, 4)) +
  theme_classic(base_size = 18)

pdf("./../Figures/METABRIC.CoxHistEcdf.pdf", width=8, height=4)
egg::ggarrange(Plt[[1]], Plt[[2]], ncol=2)
dev.off()


# 3. Simulate effects of dropouts ---------------------------------------
# A.Perform analysis with dropouts in signature genes (probably should do a full transcriptome dropout)
MB.DS <- list()
set.seed(1234)
for (i in 1:9) {
  print(i)
  cut1 <- sample(x = nrow(metabric.exp), size = nrow(metabric.exp) - (0.1*i*nrow(metabric.exp)),
                 replace = FALSE)
  print(length(cut1))
  MB.DS[[i]] <- gsva(metabric.exp[cut1, ], gset.idx.list = list("ENDORSE"=ENDORSE), 
                     method="ssgsea", kcdf="Gaussian")
}

# Make the plots 
MB.DS.Plts <- list(); MB.DS.res <- list()
for (i in 1:length(MB.DS)) {
  df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                   'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                   'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                   'ENDORSE'=MB.DS[[i]][1,])
  
  cfit = coxph(Surv(Time, Event) ~ Age + ENDORSE, data=df)
  TEMP.risk <- predict(cfit, type='risk')
  TEMP.risk.grps <- MakeRiskGrp(TEMP.risk)
  
  TEMP.Surv.df <- data.frame("Time"=metabric.annotV2$OS_MONTHS, 
                           "Event"=as.numeric(metabric.annotV2$OS_EVENT), 
                           "Age"=metabric.annotV2$AGE_AT_DIAGNOSIS, 
                           "Group"=TEMP.risk.grps)
  
  fit <- survfit(Surv(Time, Event) ~ Group, data=TEMP.Surv.df)
  MB.DS.Plts[[i]] <- ggsurvplot(fit, data=TEMP.Surv.df, palette="nejm", 
                    censor.shape="|", censor.size=1, 
                    break.time.by=60, 
                    legend="none",
                    #legend.title="ENDORSE", legend.labs=c("Low risk", "Medium risk",  "High risk"),
                    #surv.median.line = "hv",
                    risk.table=F, 
                    ggtheme=theme_classic2(base_size = 10))
  
  MB.DS.res[[i]] = coxph(Surv(Time, Event) ~ Age + ENDORSE, data=df)
}


# B. Perform analysis with reduced number of samples 
MB.L <- list(); MB.Lsam <- list()
set.seed(1234)
for (i in 1:9) {
  print(i)
  res <- list()
  cut2 <- sample(x = 1:ncol(metabric.exp), size  = 0.1*i*ncol(metabric.exp), 
                 replace=FALSE)
  MB.L[[i]] <- gsva(metabric.exp[, cut2], gset.idx.list = list("ENDORSE"=ENDORSE), 
                     method="ssgsea", kcdf="Gaussian")
  MB.Lsam[[i]] <- cut2 
}

# Make the plots 
MB.L.Plts <- list(); MB.L.res <- list()
for (i in 1:length(MB.L)) {
  MB.Lsam
  df <- data.frame('Time'=metabric.annotV2$OS_MONTHS[MB.Lsam[[i]]], 
                   'Event'=as.numeric(metabric.annotV2$OS_EVENT[MB.Lsam[[i]]]), 
                   'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS[MB.Lsam[[i]]],
                   'ENDORSE'=MB.L[[i]][1,])
  
  cfit = coxph(Surv(Time, Event) ~ Age + ENDORSE, data=df)
  TEMP.risk <- predict(cfit, type='risk')
  TEMP.risk.grps <- MakeRiskGrp(TEMP.risk)
  
  TEMP.Surv.df <- data.frame("Time"=metabric.annotV2$OS_MONTHS[MB.Lsam[[i]]], 
                             "Event"=as.numeric(metabric.annotV2$OS_EVENT)[MB.Lsam[[i]]], 
                             "Age"=metabric.annotV2$AGE_AT_DIAGNOSIS[MB.Lsam[[i]]], 
                             "Group"=TEMP.risk.grps)
  
  fit <- survfit(Surv(Time, Event) ~ Group, data=TEMP.Surv.df)
  MB.L.Plts[[i]] <- ggsurvplot(fit, data=TEMP.Surv.df, palette="nejm", 
                          censor.shape="|", censor.size=1, 
                          break.time.by=60, 
                          legend="none",
                          #legend.title="ENDORSE", legend.labs=c("Low risk", "Medium risk",  "High risk"),
                          #surv.median.line = "hv",
                          risk.table=F, 
                          ggtheme=theme_classic2(base_size = 10))
  
  MB.L.res[[i]] = coxph(Surv(Time, Event) ~ Age + ENDORSE, data=df)
}

pdf(file="./../Figures/MB_Simulations.pdf", width=7, height=6)
arrange_ggsurvplots(MB.DS.Plts, nrow=3, ncol=3)
arrange_ggsurvplots(rev(MB.L.Plts), nrow=3, ncol=3)
dev.off()


### 4. Table of Cox proportional hazards model ------------------------------
#Clinical vars
df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 'TumorSize'=MB.annot.extraV2$Tumor.Size,
                 'TumorStage'=MB.annot.extraV2$Tumor.Stage,
                 'MutationCount'=MB.annot.extraV2$Mutation.Count,
                 'LymphNode'=MB.annot.extraV2$Lymph.nodes.examined.positive)
METABRIC.Clinical.cfit <- coxph(Surv(time = Time, event = Event) ~ ., data = df)

getMeanCI <- function(x) { 
  temp <- t.test(x, na.rm=T)
  res <- data.frame("Mean"=signif(mean(x, na.rm=T), 3),
                    "N"=sum(!is.na(x)),
                    "CI"=paste(signif(temp$conf.int[1], 3), signif(temp$conf.int[2], 3), sep=" - "))
  return(res)
}

apply(df, 2, getMeanCI)

#Other classifiers 
df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 'TumorSize'=MB.annot.extraV2$Tumor.Size,
                 'TumorStage'=MB.annot.extraV2$Tumor.Stage,
                 'MutationCount'=MB.annot.extraV2$Mutation.Count,
                 'LymphNode'=MB.annot.extraV2$Lymph.nodes.examined.positive,
                 'PAM50'=MB.annot.extraV2$Pam50...Claudin.low.subtype,
                 'IntClust'=MB.annot.extraV2$Integrative.Cluster,
                 'ESR'=metabric.exp["ESR1",],
                 'PGR'=metabric.exp["PGR",],
                 'ENDORSE'=MBsc.ssgsea[1,],
                 'RiskGroups'=MB.risk.grps)

METABRIC.PAM.cfit <- coxph(Surv(time = Time, event = Event) ~ PAM50, data = df)
METABRIC.IC.cfit <- coxph(Surv(time = Time, event = Event) ~ IntClust, data = df)
METABRIC.Recptors.cfit <- coxph(Surv(time = Time, event = Event) ~ ESR + PGR, data = df)
METABRIC.ENDORSE.cfit <- coxph(Surv(time = Time, event = Event) ~ ENDORSE, data = df)
METABRIC.RG.cfit <- coxph(Surv(time = Time, event = Event) ~ RiskGroups, data = df)

CoxExport <- function(cfit, fname) {
  temp <- summary(cfit)
  df1 <- data.frame(temp$coefficients)
  df2 <- data.frame(temp$conf.int)
  c.table <- data.frame("Variable"=rownames(df1),
                        "HR"=signif(df1$exp.coef., 3), 
                        "95%CI"=paste(signif(df2$lower..95, 3), signif(df2$upper..95, 3), sep=" - "),
                        "P"=signif(df1$Pr...z.., 2))
  write.csv(c.table, file=fname)
}

CoxExport(METABRIC.PAM.cfit, fname="./../Output/METABRIC.PAM.cfit.csv")
CoxExport(METABRIC.IC.cfit, fname="./../Output/METABRIC.IC.cfit.csv")
CoxExport(METABRIC.Recptors.cfit, fname="./../Output/METABRIC.Recptors.cfit.csv")
CoxExport(METABRIC.ENDORSE.cfit, fname="./../Output/METABRIC.ENDORSE.cfit.csv")
CoxExport(METABRIC.RG.cfit, fname="./../Output/METABRIC.RG.cfit.csv")
CoxExport(METABRIC.Clinical.cfit, fname="./../Output/METABRIC.Clinical.cfit.csv")

# Multivariate analyses 
CoxExport(coxph(Surv(time = Time, event = Event) ~ PAM50 + RiskGroups, data = df), 
          fname = "./../Output/Multi_ENDORSE_PAM50.cfit.csv")
CoxExport(coxph(Surv(time = Time, event = Event) ~ IntClust + RiskGroups, data = df), 
          fname = "./../Output/Multi_ENDORSE_IntClust.cfit.csv")
CoxExport(coxph(Surv(time = Time, event = Event) ~ ESR + PGR + RiskGroups, data = df), 
          fname = "./../Output/Multi_ENDORSE_ESRPGR.cfit.csv")


#### Negative control: METABRIC ER- ----------------------------------
x <- match(colnames(metabric), metabric.annot$PATIENT_ID)
metabric_N.annot <- metabric.annot[x, ]

OS_EVENT <- array(dim=nrow(metabric_N.annot))
OS_EVENT[which(metabric_N.annot$OS_STATUS == "LIVING")] <- 0
OS_EVENT[which(metabric_N.annot$OS_STATUS == "DECEASED")] <- 1

metabric_N.annot$OS_EVENT <- OS_EVENT
y <- which(metabric_N.annot$ER_IHC == "Negative")
metabric_N.annot <- metabric_N.annot[y,]
metabric.N <- metabric[, y]
metabric.N <- scale(t(metabric.N))

ENDORSE <- readRDS("./../Output/ENDORSE.RDS")

MB_N.ssgsea <- gsva(t(metabric.N), list("ENDORSE"=ENDORSE),
                    verbose=TRUE, kcdf="Gaussian", method='ssgsea', ssgsea.norm=T)
MB_N.ssgsea[1,] <- ReScale(MB_N.ssgsea[1,])


# Predict risk groups 
# Train df 
df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 'ENDORSE'=MBsc.ssgsea[1,])

cfit = coxph(Surv(Time, Event) ~ ENDORSE, data=df)

# predict
MB_N.pred <- predict(cfit, data.frame("ENDORSE"=MB_N.ssgsea[1,]), type="risk")

# compare
# plot - by risk groups
MB_N.pred.grp <- MakeRiskGrp(MB_N.pred, thresh=c(1,2))

df <- data.frame("Time"=as.numeric(metabric_N.annot$OS_MONTHS),
                 "Event"=metabric_N.annot$OS_EVENT,
                 "RiskGroup"=MB_N.pred.grp, 
                 "Risk"=MB_N.pred, 
                 "ENDORSE"=MB_N.ssgsea[1,])

# Survival plots ----
library(survminer)
fit.os <- survfit(Surv(Time, Event) ~ RiskGroup, data=df)

surv_median(fit.os)

P1 <- ggsurvplot(fit.os, data=df, palette="nejm", 
                 censor.shape="|", censor.size=2, 
                 break.time.by=60, 
                 legend.title="ENDORSE", legend.labs=c("Low risk", "Medium risk",  "High risk"),
                 risk.table=F, 
                 title="ER- METABRIC breast cancers",
                 ggtheme=theme_classic2(base_size = 18))


summary(coxph(Surv(Time, Event) ~ RiskGroup, data=df))

pdf("./../Figures/METABRIC_ER_Neg.pdf", width=6, height=6)
print(P1)
dev.off()

### PATHWAYS and MUTATIONS in RISK GROUPS ------------------------------
cfit <- readRDS(file="./../Output/ENDORSE.CoxFit.RDS") 
MB.risk <- predict(cfit, type='risk')
MB.risk.grps <- MakeRiskGrp(MB.risk)
names(MB.risk.grps) <- colnames(MBsc.ssgsea)

## Somatic mutations with synonymous mutations removed ### 
x <- match(metabric.mut$Tumor_Sample_Barcode, colnames(metabric.exp))
MB.mut <- metabric.mut[which(!is.na(x)), ]
length(unique(MB.mut$Hugo_Symbol)) #169

y <- match(MB.mut$Tumor_Sample_Barcode, names(MB.risk.grps))
MB.mut$RiskGroup <- MB.risk.grps[y]

MBRisk.Mutation.Dist <- as.data.frame.matrix(table(MB.mut$Hugo_Symbol, MB.mut$RiskGroup))
RiskN <- as.matrix(table(MB.risk.grps))
MutF <- apply(MBRisk.Mutation.Dist, 1, function(x) 100*sum(x)/sum(RiskN))
z <- which(MutF > 5)

Fisher.Res <- array(dim=length(z))
Chi.Res <- array(dim=length(z))
for (i in 1:length(z)) {
  # R1 = mut, R2 = no mut, C1 = Low/Med, C2 = High
  R1.C1 = MBRisk.Mutation.Dist[z[i], 1] +  MBRisk.Mutation.Dist[z[i], 2]
  R1.C2 = MBRisk.Mutation.Dist[z[i], 3]
  R2.C1 = (RiskN[1] + RiskN[2]) - R1.C1
  R2.C2 = RiskN[3] - R1.C2
  
  mat2 <- matrix(c(R1.C1, R2.C1, R1.C2, R2.C2), nrow=2)
  Fisher.Res[i] <- fisher.test(mat2)$p.value
  Chi.Res[i] <- chisq.test(mat2)$p.value
  print(i)
}

k <- which(Chi.Res < 0.05)
MBRisk.Mutation.Dist[z[k],]

Col.Set <- list()
Col.Set[[1]] <- c("#DF6868", "#803B3B")
Col.Set[[2]] <- c("#B6DF68", "#68803B")
Col.Set[[3]] <- c("#68DFD3", "#3B8078")
Col.Set[[4]] <- c("#686FDF", "#3B4080")
Col.Set[[5]] <- c("#DA68DF", "#7C3B80")

pdf("./../Figures/METABRIC_MutF.pdf", width=8, height=4)
par(mar=c(6.5,5.5,1,1), mfrow=c(1,5), cex.lab=1.5, cex.axis=1.5)
i = 1
f1 <- 100*(MBRisk.Mutation.Dist[z[k][i], 1] + MBRisk.Mutation.Dist[z[k][i], 2]) / (RiskN[1] + RiskN[2])
f2 <- 100*(MBRisk.Mutation.Dist[z[k][i], 3]) / RiskN[3]
print(c(f1, f2))
barplot(c(f1, f2), width=0.5, space=0.25, col=Col.Set[[i]], border=NA, axes=F, ylim=c(0,60), ylab="Samples Mutated (%)")
axis(side=2, at=c(0,10,20,30,40,50,60), lwd.ticks = 2, lwd = 2)
axis(side=1, at=c(0.375,1), labels=c("Low/Med", "High"), lwd=0, las=2)

for (i in 2:length(k)) {
  f1 <- 100*(MBRisk.Mutation.Dist[z[k][i], 1] + MBRisk.Mutation.Dist[z[k][i], 2]) / (RiskN[1] + RiskN[2])
  f2 <- 100*(MBRisk.Mutation.Dist[z[k][i], 3]) / RiskN[3]
  print(c(f1, f2))
  barplot(c(f1, f2), width=0.5, space=0.25, col=Col.Set[[i]], border=NA, axes=F, ylim=c(0,60))
  #axis(side=2, at=c(0,10,20,30,40,50,60), lwd.ticks = 2, lwd = 2)
  axis(side=1, at=c(0.375,1), labels=c("Low/Med", "High"), lwd=0, las=2)
}
dev.off()

which(p.adjust(Chi.Res) < 0.05)
MBRisk.Mutation.Dist[z[which(p.adjust(Chi.Res) < 0.05)],]

write.csv(data.frame(MBRisk.Mutation.Dist[z,], "Chi.sq(P)"=Chi.Res, "FDR"=p.adjust(Chi.Res)), 
          "./../Output/MBRiskMutationTable.csv")

## CNAs ### 
y <- match(colnames(MBsc.ssgsea), colnames(metabric.cnv))
MB.cnv <- data.frame(metabric.cnv)[, which(!is.na(y))]
rownames(MB.cnv) <- metabric.cnv$Hugo_Symbol

# R1 = mut, R2 = no mut, C1 = Low/Med, C2 = High
I_lm <- which(MB.risk.grps != "High Risk")
I_h <- which(MB.risk.grps == "High Risk")

MB.cnv <- na.omit(MB.cnv)
CNV_Loss_Chi.Res <- array(dim=nrow(MB.cnv))
CNV_Gain_Chi.Res <- array(dim=nrow(MB.cnv))
for (i in 1:nrow(MB.cnv)) {
  R1.C1 = sum(MB.cnv[i, I_lm] < 0)
  R1.C2 = sum(MB.cnv[i, I_h] < 0)
  R2.C1 = (RiskN[1] + RiskN[2]) - R1.C1
  R2.C2 = RiskN[3] - R1.C2
  
  mat2 <- matrix(c(R1.C1, R2.C1, R1.C2, R2.C2), nrow=2)
  CNV_Loss_Chi.Res[i] <- chisq.test(mat2)$p.value
  print(i)
}

library(foreach)
library(doParallel)

CNV_Gain_Chi.Res <- foreach (i=1:nrow(MB.cnv), .combine='c') %dopar% {
  R1.C1 = sum(MB.cnv[i, I_lm] > 0)
  R1.C2 = sum(MB.cnv[i, I_h] > 0)
  R2.C1 = (RiskN[1] + RiskN[2]) - R1.C1
  R2.C2 = RiskN[3] - R1.C2
  
  mat2 <- matrix(c(R1.C1, R2.C1, R1.C2, R2.C2), nrow=2)
  print(i)
  return(chisq.test(mat2)$p.value)
}

length(which(CNV_Loss_Chi.Res < 0.05)) #831
length(which(p.adjust(CNV_Loss_Chi.Res) < 0.05)) #0

length(which(CNV_Gain_Chi.Res < 0.05)) #831
length(which(p.adjust(CNV_Gain_Chi.Res) < 0.05)) #0


## Curated and Hallmark pathways ### -------------------------------
y <- match(colnames(MBsc.ssgsea), colnames(metabric.c2))
metabric.c2 <- metabric.c2[,y]

y <- match(colnames(MBsc.ssgsea), colnames(metabric.c2))
metabric.hall <- metabric.hall[,y]

I_lm <- which(MB.risk.grps != "High Risk")
I_h <- which(MB.risk.grps == "High Risk")

MBRisk.hall.p <- apply(metabric.hall, 1, function (x) t.test(x[I_lm], x[I_h])$p.value)
MBRisk.c2.p  <- apply(metabric.c2, 1, function (x) t.test(x[I_lm], x[I_h])$p.value)

MBRisk.hall.e <- apply(metabric.hall, 1, function (x) {
  (mean(x[I_h]) - mean(x[I_lm]))/sd(x)
})

MBRisk.c2.e  <- apply(metabric.c2, 1, function (x)  {
  (mean(x[I_h]) - mean(x[I_lm]))/sd(x)
})

write.csv(data.frame("Effect"=MBRisk.hall.e, 
                     "P-value"=MBRisk.hall.p,
                     "FDR"=p.adjust(MBRisk.hall.p)),
          file="./../Output/MBRisk.hall.csv")

write.csv(data.frame("Effect"=MBRisk.c2.e, 
                     "P-value"=MBRisk.c2.p,
                     "FDR"=p.adjust(MBRisk.c2.p)),
          file="./../Output/MBRisk.c2.csv")


# Important pathways 
HRC2Dn <- names(which(MBRisk.c2.e < -0.8 & p.adjust(MBRisk.c2.p) < 0.05)) #12
HRC2Up <- names(which(MBRisk.c2.e > 0.8 & p.adjust(MBRisk.c2.p) < 0.05)) #229

write.csv(MBRisk.c2.p[HRC2Up], "./../Output/HRC2Up.csv")
write.csv(MBRisk.c2.p[HRC2Dn], "./../Output/HRC2Dn.csv")

## Heatmaps of some pathways based on common themes 
library(RColorBrewer)
library(stringr)
coul <- colorRampPalette(c("navy", "white", "darkred"))(15)

pdf("./../Figures/MBRisk_Pathway_Heatmaps.pdf", height=8, width=12)
heatmap(metabric.c2[grep("P53", HRC2Up, value=T), c(I_lm, I_h)], 
        scale = "row",
        col = coul, 
        cexRow = 1.25,
        Colv = NA, 
        margins = c(25,25),
        ColSideColors = c(rep("darkred", length(I_lm)), rep("salmon", length(I_h))),
        labCol= NA,
        labRow = str_wrap(gsub("_"," ", grep("P53", HRC2Up, value=T)), width=50))


l1 <- grep("RB", HRC2Up, value=T)
l1 <- l1[-2]
heatmap(metabric.c2[l1, c(I_lm, I_h)], 
        scale = "row",
        col = coul, 
        cexRow = 1.25,
        Colv = NA, 
        margins = c(35,25),
        ColSideColors = c(rep("darkred", length(I_lm)), rep("salmon", length(I_h))),
        labCol= NA,
        labRow = str_wrap(gsub("_"," ", l1), width=50))

l1 <- grep("CELL_CYCLE", HRC2Up, value=T)
l1 <- l1[-c(5,6,15,16,17)]
heatmap(metabric.c2[l1, c(I_lm, I_h)], 
        scale = "row",
        col = coul, 
        cexRow = 1.25,
        Colv = NA, 
        margins = c(25,25),
        ColSideColors = c(rep("darkred", length(I_lm)), rep("salmon", length(I_h))),
        labCol= NA,
        labRow = str_wrap(gsub("_"," ", l1), width=50))


l1 <- c(grep("RECOMBINATION", HRC2Up, value=T), 
        grep("REPAIR", HRC2Up, value=T),
        grep("ATR", HRC2Up, value=T),
        grep("ATM", HRC2Up, value=T))
l1 <- l1[-11]
heatmap(metabric.c2[l1, c(I_lm, I_h)], 
        scale = "row",
        col = coul, 
        cexRow = 1.25,
        Colv = NA, 
        margins = c(25,25),
        ColSideColors = c(rep("darkred", length(I_lm)), rep("salmon", length(I_h))),
        labCol= NA,
        labRow = str_wrap(gsub("_"," ", l1), width=50))

heatmap(metabric.c2[HRC2Up, c(I_lm, I_h)], 
        scale = "row",
        col = coul, 
        cexRow = 0.25,
        Colv = NA, 
        margins = c(5,25),
        ColSideColors = c(rep("darkred", length(I_lm)), rep("salmon", length(I_h))),
        labCol= NA,
        labRow = NA)
#labRow = str_wrap(gsub("_"," ", HRC2Up), width=150))

legend_image <- as.raster(matrix(coul, nrow=1))
plot.new()
rasterImage(legend_image, 0, 0.5, 1, 1)
text(y=0.35, x= seq(0, 1, l=5), cex=2,
     labels = round(quantile(scale(t(metabric.c2[HRC2Up, c(I_lm, I_h)]))), 1))

heatmap(metabric.c2[HRC2Dn, c(I_lm, I_h)], 
        scale = "row",
        col = coul, 
        cexRow = 1.25,
        Colv = NA, 
        margins = c(25,25),
        ColSideColors = c(rep("darkred", length(I_lm)), rep("salmon", length(I_h))),
        labCol= NA,
        labRow = str_wrap(gsub("_"," ", HRC2Dn), width=50))

legend_image <- as.raster(matrix(coul, nrow=1))
plot.new()
rasterImage(legend_image, 0, 0.5, 1, 1)
text(y=0.35, x= seq(0, 1, l=5), cex=2,
     labels = round(quantile(scale(t(metabric.c2[HRC2Dn, c(I_lm, I_h)]))), 1))

dev.off()


#### Make individual plots 
library(stringr)
risk.labs <- array(dim=length(MB.risk.grps))
risk.labs[I_lm] <- "Low/Med"
risk.labs[I_h] <- "High"

# P53,  RB, DNA damage, cell cycle
Plts <- list()
l1 <- c("SCIAN_CELL_CYCLE_TARGETS_OF_TP53_AND_TP73_DN",
        "MARKEY_RB1_ACUTE_LOF_UP",
        "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
        "REACTOME_CELL_CYCLE_MITOTIC")
df <- data.frame(t(metabric.c2[l1,]), risk.labs)
for (i in 1:4) {
  Plts[[i]] <- ggviolin(df, x="risk.labs", y=l1[i], color="risk.labs", 
                        add="jitter") +
    labs(y="ssGSEA score", title=str_wrap(gsub("_"," ", l1[i]), width=30))
}

pdf(file="./../Figures/MBRisk_P53RBDNACC.pdf", width=12, height=4)
ggpubr::ggarrange(plotlist=ggpar(Plts, xlab = FALSE,  font.main=c(9,"bold","black"), legend="none"), nrow=1)
dev.off()

# Metastasis and others 
Plts <- list()
l1 <- c("RAMASWAMY_METASTASIS_UP",
        "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP",
        "YU_MYC_TARGETS_UP",
        "STEIN_ESRRA_TARGETS_RESPONSIVE_TO_ESTROGEN_DN")

df <- data.frame(t(metabric.c2[l1,]), risk.labs)
for (i in 1:4) {
  Plts[[i]] <- ggviolin(df, x="risk.labs", y=l1[i], color="risk.labs", 
                        add="jitter") +
    labs(y="ssGSEA score", title=str_wrap(gsub("_"," ", l1[i]), width=30))
}

pdf(file="./../Figures/MBRisk_MetaSign.pdf", width=12, height=4)
ggpubr::ggarrange(plotlist=ggpar(Plts, xlab = FALSE,  font.main=c(9,"bold","black"), legend="none"), nrow=1)
dev.off()

# Tumor-ECM interactions
Plts <- list()
l1 <- c("REACTOME_TYPE_I_HEMIDESMOSOME_ASSEMBLY",
        "PID_INTEGRIN4_PATHWAY",
        "REACTOME_LAMININ_INTERACTIONS",
        "NABA_BASEMENT_MEMBRANES")

df <- data.frame(t(metabric.c2[l1,]), risk.labs)
for (i in 1:4) {
  Plts[[i]] <- ggviolin(df, x="risk.labs", y=l1[i], color="risk.labs", 
                        add="jitter") +
    labs(y="ssGSEA score", title=str_wrap(gsub("_"," ", l1[i]), width=30))
}

pdf(file="./../Figures/MBRisk_ECM.pdf", width=12, height=4)
ggpubr::ggarrange(plotlist=ggpar(Plts, xlab = FALSE,  font.main=c(9,"bold","black"), legend="none"), nrow=1)
dev.off()

#### Top 5 pathways by database 
# REACTOME
Plts <- list()
df <- data.frame(t(metabric.c2[l1,]), risk.labs)
for (i in 1:5) {
  Plts[[i]] <- ggviolin(df, x="risk.labs", y=l1[i], color="risk.labs", 
                        add="jitter") +
    labs(y="ssGSEA score", title=str_wrap(gsub("_"," ", l1[i]), width=25))
}

pdf(file="./../Figures/MBRisk_REACTOME.pdf", width=15, height=5)
ggpubr::ggarrange(plotlist=ggpar(Plts, xlab = FALSE,  font.main=c(10,"bold","black"), legend="none"), nrow=1)
dev.off()

# BIOCARTA
Plts <- list()
df <- data.frame(t(metabric.c2[l2,]), risk.labs)
for (i in 1:5) {
  Plts[[i]] <- ggviolin(df, x="risk.labs", y=l2[i], color="risk.labs", 
                        add="jitter") +
    labs(y="ssGSEA score", title=str_wrap(gsub("_"," ", l2[i]), width=40))
}

pdf(file="./../Figures/MBRisk_BIOCARTA.pdf", width=15, height=5)
ggpubr::ggarrange(plotlist=ggpar(Plts, xlab = FALSE,  font.main=c(10,"bold","black"), legend="none"), nrow=1)
dev.off()

# KEGG
Plts <- list()
df <- data.frame(t(metabric.c2[l3,]), risk.labs)
for (i in 1:5) {
  Plts[[i]] <- ggviolin(df, x="risk.labs", y=l3[i], color="risk.labs", 
                        add="jitter") +
    labs(y="ssGSEA score", title=str_wrap(gsub("_"," ", l3[i]), width=25))
}

pdf(file="./../Figures/MBRisk_KEGG.pdf", width=15, height=5)
ggpubr::ggarrange(plotlist=ggpar(Plts, xlab = FALSE,  font.main=c(10,"bold","black"), legend="none"), nrow=1)
dev.off()

# PID
Plts <- list()
df <- data.frame(t(metabric.c2[l4,]), risk.labs)
for (i in 1:5) {
  Plts[[i]] <- ggviolin(df, x="risk.labs", y=l4[i], color="risk.labs", 
                        add="jitter") +
    labs(y="ssGSEA score", title=str_wrap(gsub("_"," ", l4[i]), width=40))
}

pdf(file="./../Figures/MBRisk_PID.pdf", width=15, height=5)
ggpubr::ggarrange(plotlist=ggpar(Plts, xlab = FALSE, legend="none", font.main=c(10,"bold","black")), nrow=1)
dev.off()


######### Compare meta-PCNA proliferation signature -----------
install.packages("ProliferativeIndex")
library(ProliferativeIndex)

temp <- readDataForPI(data.frame(metabric.exp), modelIDs = ENDORSE)
metabric.PI <- calculatePI(temp)

summary(lm(MB.risk ~ metabric.PI))
plot(MB.risk, metabric.PI)
abline(lm(metabric.PI ~ MB.risk))

df <- data.frame("Time"=as.numeric(metabric.annotV2$OS_MONTHS),
                 "Event"=as.numeric(metabric.annotV2$OS_EVENT),
                 "RiskGroup"=MB.risk.grps, 
                 "Risk"=MB.risk, 
                 "PI"=metabric.PI)

F1 <- coxph(Surv(Time, Event) ~ PI + RiskGroup, data=df)
F2 <- coxph(Surv(Time, Event) ~ PI + Risk, data=df)

CoxExport(F1, fname="./../Output/METABRIC.ProlIndex1.cfit.csv")
CoxExport(F2, fname="./../Output/METABRIC.ProlIndex2.cfit.csv")

CoxExport <- function(cfit, fname) {
  temp <- summary(cfit)
  df1 <- data.frame(temp$coefficients)
  df2 <- data.frame(temp$conf.int)
  c.table <- data.frame("Variable"=rownames(df1),
                        "HR"=signif(df1$exp.coef., 3), 
                        "95%CI"=paste(signif(df2$lower..95, 3), signif(df2$upper..95, 3), sep=" - "),
                        "P"=signif(df1$Pr...z.., 2))
  write.csv(c.table, file=fname)
}


summary(coxph(Surv(Time, Event) ~ PI, data=df))
summary(coxph(Surv(Time, Event) ~ PI + RiskGroup, data=df))


