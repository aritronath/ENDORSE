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


### Match matrices and train the model on batch-corrected data ####-------------------------
library(sva)
temp <- DoCombat(metabric.s, t(gse124647.s))
train <- match(colnames(metabric.s), colnames(temp))
test <- match(rownames(gse124647.s), colnames(temp))

MBGSE124647.ENDORSE <- gsva(temp[, train], gset.idx.list = list("ENDORSE"=ENDORSE, 
                                                                "HALLMARK_ESTROGEN_RESPONSE_EARLY" = H_ESTR_EARLY, 
                                                                "HALLMARK_ESTROGEN_RESPONSE_LATE" = H_ESTR_LATE), 
                            method="ssgsea", kcdf="Gaussian")



df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 t(MBGSE124647.ENDORSE))

cfit.0 = coxph(Surv(Time, Event) ~ ENDORSE, data=df)
cfit.5 = coxph(Surv(Time, Event) ~ HALLMARK_ESTROGEN_RESPONSE_EARLY, data=df)
cfit.6 = coxph(Surv(Time, Event) ~ ENDORSE + HALLMARK_ESTROGEN_RESPONSE_EARLY, data=df, x = T)

### ER response and  survival analysis -------
gse124647.ssgsea <- gsva(temp[, test], list("ENDORSE"=ENDORSE,
                                              "HALLMARK_ESTROGEN_RESPONSE_EARLY" = H_ESTR_EARLY, 
                                              "HALLMARK_ESTROGEN_RESPONSE_LATE" = H_ESTR_LATE),
                         verbose=TRUE, kcdf="Gaussian", method='ssgsea', ssgsea.norm=T)

gse124647.ssgsea2 <- rbind(gsva(temp[, test], Hset, method='ssgsea', kcdf='Gaussian'), 
                          gsva(temp[, test], c6set, method='ssgsea', kcdf='Gaussian'))

gse124647.os <- Surv(time=as.numeric(gse124647.annotV2$`oas months`), event=as.numeric(gse124647.annotV2$`os event`))
gse124647.pfs <- Surv(time=as.numeric(gse124647.annotV2$`pfs months`), event=as.numeric(gse124647.annotV2$pfs.event))

gse124647.ssgsea2 <- rbind(gsva(temp[, test], Hset, method='ssgsea', kcdf='Gaussian'), 
                           gsva(temp[, test], c6set, method='ssgsea', kcdf='Gaussian'))

# predict
gse124647.pred <- predict(cfit.6, data.frame(t(gse124647.ssgsea)), type="risk")
gse124647.pred.grp <- MakeRiskGrp(gse124647.pred)

############# Apply TransCONFIRM signature #########################
gse124647.TransCONFIRM <- TransCONFIRM(gse124647.s)

############# Apply SET signature #########################
gse124647.SET <- SET_Score(gse124647.s)

x <- match(names(gse124647.SET), gse124647.annot$X16)
annot.SET <- as.numeric(sapply(strsplit(as.character(gse124647.annot[x,15]), " "), "[", 4))

cor(gse124647.SET, annot.SET, use='complete.obs')
cor.test(gse124647.SET, annot.SET, use='complete.obs')
plot(gse124647.SET, annot.SET)


df.gse124647 <- data.frame("OSTime"=as.numeric(gse124647.annotV2$`oas months`),
                           "OSEvent"=as.numeric(gse124647.annotV2$`os event`), 
                           "PFSTime"=as.numeric(gse124647.annotV2$`pfs months`),
                           "PFSEvent"=as.numeric(gse124647.annotV2$pfs.event),
                           "RiskGroup"=gse124647.pred.grp, 
                           "Risk"=gse124647.pred, 
                           "SET"=gse124647.SET,
                           "TransCONFIRM"=gse124647.TransCONFIRM)

gse124647.SET_Cat <- array(data=0, dim = length(gse124647.SET))
gse124647.SET_Cat[which(gse124647.SET > median(gse124647.SET))] <- 1
df.gse124647$SET_Cat <- gse124647.SET_Cat

cpOS.E2.gse124647 <- coxph(Surv(OSTime, OSEvent) ~ RiskGroup, data=df.gse124647, x=T)
cpOS.SET.gse124647 <- coxph(Surv(OSTime, OSEvent) ~ SET, data=df.gse124647, x=T)
cpOS.TC.gse124647 <- coxph(Surv(OSTime, OSEvent) ~ TransCONFIRM, data=df.gse124647, x=T)

cpPFS.E2.gse124647 <- coxph(Surv(PFSTime, PFSEvent) ~ RiskGroup, data=df.gse124647, x=T)
cpPFS.SET.gse124647 <- coxph(Surv(PFSTime, PFSEvent) ~ SET, data=df.gse124647, x=T)
cpPFS.TC.gse124647 <- coxph(Surv(PFSTime, PFSEvent) ~ TransCONFIRM, data=df.gse124647, x=T)

summary(cpOS.E2.gse124647)
cpOS.SET.gse124647
cpOS.TC.gse124647

summary(cpPFS.E2.gse124647)
cpPFS.SET.gse124647
cpPFS.TC.gse124647

plrtest(cpOS.E2.gse124647, cpOS.SET.gse124647, nested=F)
plrtest(cpOS.E2.gse124647, cpOS.TC.gse124647, nested=F)
        
plrtest(cpPFS.E2.gse124647, cpPFS.SET.gse124647, nested=F)
plrtest(cpPFS.E2.gse124647, cpPFS.TC.gse124647, nested=F)

##### Survival plots --------------------------------------------------
library(survminer)

fitOS.E2.gse124647 <- survfit(Surv(OSTime, OSEvent) ~ RiskGroup, data=df.gse124647)
fitOS.SET.gse124647 <- survfit(Surv(OSTime, OSEvent) ~ SET_Cat, data=df.gse124647)
fitOS.TC.gse124647 <- survfit(Surv(OSTime, OSEvent) ~ TransCONFIRM, data=df.gse124647)

fitPFS.E2.gse124647 <- survfit(Surv(PFSTime, PFSEvent) ~ RiskGroup, data=df.gse124647)
fitPFS.SET.gse124647 <- survfit(Surv(PFSTime, PFSEvent) ~ SET_Cat, data=df.gse124647)
fitPFS.TC.gse124647 <- survfit(Surv(PFSTime, PFSEvent) ~ TransCONFIRM, data=df.gse124647)


P1 <- ggsurvplot(fitOS.E2.gse124647, data=df.gse124647, palette="jco", 
                 censor.shape="|", censor.size=2, 
                 break.time.by=12, 
                 legend.title="ENDORSE", legend.labs=c("Low risk", "Medium risk", "High risk"),
                 surv.median.line = "hv",
                 risk.table=T, 
                 #title="ER+/HER2- mBC OS",
                 ggtheme=theme_classic2(base_size = 14))

P2 <- ggsurvplot(fitOS.SET.gse124647, data=df.gse124647, palette="jco", 
                 censor.shape="|", censor.size=2, 
                 break.time.by=12, 
                 legend.title="SET", legend.labs=c("Low SET", "High SET"),
                 surv.median.line = "hv",
                 risk.table=T, 
                 #title="ER+/HER2- mBC OS",
                 ggtheme=theme_classic2(base_size = 14))


P3 <- ggsurvplot(fitOS.TC.gse124647, data=df.gse124647, palette="jco", 
                 censor.shape="|", censor.size=2, 
                 break.time.by=12, 
                 legend.title="TransCONFIRM", legend.labs=c("TC Low", "TC High"),
                 surv.median.line = "hv",
                 risk.table=T, 
                 #title="ER+/HER2- mBC PFS",
                 ggtheme=theme_classic2(base_size = 14))

P4 <- ggsurvplot(fitPFS.E2.gse124647, data=df.gse124647, palette="jco", 
                 censor.shape="|", censor.size=2, 
                 break.time.by=12, 
                 legend.title="ENDORSE", legend.labs=c("Low risk", "Medium risk", "High risk"),
                 surv.median.line = "hv",
                 risk.table=T, 
                 #title="ER+/HER2- mBC PFS",
                 ggtheme=theme_classic2(base_size = 14))

P5 <- ggsurvplot(fitPFS.SET.gse124647, data=df.gse124647, palette="jco", 
                 censor.shape="|", censor.size=2, 
                 break.time.by=12, 
                 legend.title="SET", legend.labs=c("Low SET", "High SET"),
                 surv.median.line = "hv",
                 risk.table=T, 
                 #title="ER+/HER2- mBC PFS",
                 ggtheme=theme_classic2(base_size = 14))


P6 <- ggsurvplot(fitPFS.TC.gse124647, data=df.gse124647, palette="jco", 
                 censor.shape="|", censor.size=2, 
                 break.time.by=12, 
                 legend.title="TransCONFIRM", legend.labs=c("TC Low", "TC High"),
                 surv.median.line = "hv",
                 risk.table=T, 
                 #title="ER+/HER2- mBC PFS",
                 ggtheme=theme_classic2(base_size = 14))



pdf("./../Figures/Prediction_GSE124647_os_pfs_ENDORSEv2.pdf", width=6, height=6)
print(P1)
print(P2)
print(P3)
print(P4)
print(P5)
print(P6)
dev.off()


survdiff(Surv(OSTime, OSEvent) ~ RiskGroup, data=df.gse124647)
survdiff(Surv(OSTime, OSEvent) ~ SET_Cat, data=df.gse124647)
survdiff(Surv(OSTime, OSEvent) ~ TransCONFIRM, data=df.gse124647)

survdiff(Surv(PFSTime, PFSEvent) ~ RiskGroup, data=df.gse124647)
survdiff(Surv(PFSTime, PFSEvent) ~ SET_Cat, data=df.gse124647)
survdiff(Surv(PFSTime, PFSEvent) ~ TransCONFIRM, data=df.gse124647)

#### Scatter plot 
gse124647.gams <- apply(t(gse124647.ssgsea2), 2, GamApp, gse124647.pred)
write.csv(t(gse124647.gams), file="./../Output/gse124647.gams.csv")

tail(sort(gse124647.gams[3,]), 20)

df <- data.frame(cbind(t(gse124647.ssgsea2), "Risk"=gse124647.pred))  
PPL <- list()
PPL[[1]] <- ggplot(df, aes(x=HALLMARK_G2M_CHECKPOINT, y=Risk, color=HALLMARK_G2M_CHECKPOINT))  +
  geom_point(size=3) +
  scale_colour_gradient2(low = "blue", mid = "grey80", high = "red", name="GES") +
  geom_smooth(method='gam', se=T, color='navy', size=3) +
  labs(x="Cell cycle progression", y="ENDORSE") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(size=14, face="bold"),
        title = element_text(size=16, face="bold", color="grey25"))

PPL[[2]] <- ggplot(df, aes(x=E2F1_UP.V1_UP, y=Risk, color=E2F1_UP.V1_UP))  +
  geom_point(size=3) +
  scale_colour_gradient2(low = "blue", mid = "grey80", high = "red", name="GES") +
  geom_smooth(method='gam', se=T, color='navy', size=3) +
  labs(x="E2F1 upregulation", y="ENDORSE") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(size=14, face="bold"),
        title = element_text(size=16, face="bold", color="grey25"))

PPL[[3]] <- ggplot(df, aes(x=MTOR_UP.V1_UP, y=Risk, color=MTOR_UP.V1_UP))  +
  geom_point(size=3) +
  scale_colour_gradient2(low = "blue", mid = "grey80", high = "red", name="GES") +
  geom_smooth(method='gam', se=T, color='navy', size=3) +
  labs(x="MTOR upregulation", y="ENDORSE") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(size=14, face="bold"),
        title = element_text(size=16, face="bold", color="grey25"))

pdf("./../Figures/GSE124647.pathways.pdf", height =  4, width=15)
egg::ggarrange(plots=PPL, nrow=1)
dev.off()


######### Compare with the empirical and estrogen response signatures ############
gse124647.pred.c0 <- predict(cfit.0, data.frame(t(gse124647.ssgsea)), type="risk")
gse124647.pred.grp.c0 <- MakeRiskGrp(gse124647.pred.c0)

gse124647.pred.c5 <- predict(cfit.5, data.frame(t(gse124647.ssgsea)), type="risk")
gse124647.pred.grp.c5 <- MakeRiskGrp(gse124647.pred.c5)

df.gse124647.2 <- data.frame("OSTime"=as.numeric(gse124647.annotV2$`oas months`),
                           "OSEvent"=as.numeric(gse124647.annotV2$`os event`), 
                           "PFSTime"=as.numeric(gse124647.annotV2$`pfs months`),
                           "PFSEvent"=as.numeric(gse124647.annotV2$pfs.event),
                           "RiskGroup"=gse124647.pred.grp,
                           "Risk"=gse124647.pred,
                           "RiskGroup0"=gse124647.pred.grp.c0, 
                           "Risk0"=gse124647.pred.c0, 
                           "RiskGroup5"=gse124647.pred.grp.c5, 
                           "Risk5"=gse124647.pred.c5) 
                           

cpPFS.c0.gse124647 <- coxph(Surv(PFSTime, PFSEvent) ~ Risk0, data=df.gse124647.2, x=T) #7.735e-05
cpPFS.c5.gse124647 <- coxph(Surv(PFSTime, PFSEvent) ~ Risk5, data=df.gse124647.2, x=T) #0.01597
cpPFS.c6.gse124647 <- coxph(Surv(PFSTime, PFSEvent) ~ Risk, data=df.gse124647.2, x=T) #4.184e-06

cpOS.c0.gse124647 <- coxph(Surv(OSTime, OSEvent) ~ Risk0, data=df.gse124647.2, x=T) #0.0003096
cpOS.c5.gse124647 <- coxph(Surv(OSTime, OSEvent) ~ Risk5, data=df.gse124647.2, x=T) #0.0254
cpOS.c6.gse124647 <- coxph(Surv(OSTime, OSEvent) ~ Risk, data=df.gse124647.2, x=T) #4.782e-05

plrtest(cpOS.c6.gse124647, cpOS.c0.gse124647, nested=T) #H1 p = 0.0208 #H1A p = 0.036
plrtest(cpOS.c6.gse124647, cpOS.c5.gse124647, nested=T) #H1 p = 0.000199 #H1A p = 0.000472

plrtest(cpPFS.c6.gse124647, cpPFS.c0.gse124647, nested=T) #H1 p = 0.00484 #H1A p = 0.0109
plrtest(cpPFS.c6.gse124647, cpPFS.c5.gse124647, nested=T) #H1 p = 6e-04 #H1A p = 0.00026

fitOS.c0.gse124647 <- survfit(Surv(OSTime, OSEvent) ~ RiskGroup0, data=df.gse124647.2)
fitOS.c5.gse124647 <- survfit(Surv(OSTime, OSEvent) ~ RiskGroup5, data=df.gse124647.2)

fitPFS.c0.gse124647 <- survfit(Surv(PFSTime, PFSEvent) ~ RiskGroup0, data=df.gse124647.2)
fitPFS.c5.gse124647 <- survfit(Surv(PFSTime, PFSEvent) ~ RiskGroup5, data=df.gse124647.2)

pdf("./../Figures/GSE124647_Emp_Hstr.pdf", width=6, height=6)
ggsurvplot(fitOS.c0.gse124647, data=df.gse124647.2, palette="jco", 
                 censor.shape="|", censor.size=2, 
                 break.time.by=12, 
                 legend.title="Empirical", legend.labs=c("Low risk", "Medium risk", "High risk"),
                 surv.median.line = "hv",
                 ggtheme=theme_classic2(base_size = 14))

ggsurvplot(fitOS.c5.gse124647, data=df.gse124647.2, palette="jco", 
           censor.shape="|", censor.size=2, 
           break.time.by=12, 
           legend.title="Estrogen response", legend.labs=c("Low risk", "Medium risk"),
           surv.median.line = "hv",
           ggtheme=theme_classic2(base_size = 14))

ggsurvplot(fitPFS.c0.gse124647, data=df.gse124647.2, palette="jco", 
           censor.shape="|", censor.size=2, 
           break.time.by=12, 
           legend.title="Empirical", legend.labs=c("Low risk", "Medium risk", "High risk"),
           surv.median.line = "hv",
           ggtheme=theme_classic2(base_size = 14))

ggsurvplot(fitPFS.c5.gse124647, data=df.gse124647.2, palette="jco", 
           censor.shape="|", censor.size=2, 
           break.time.by=12, 
           legend.title="Estrogen response", legend.labs=c("Low risk", "Medium risk"),
           surv.median.line = "hv",
           ggtheme=theme_classic2(base_size = 14))
dev.off()

