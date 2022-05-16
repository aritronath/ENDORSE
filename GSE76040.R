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

### Match matrices and train the model on batch-corrected data ####-------------------------
library(sva)
temp <- DoCombat(metabric.s, t(gse76040.s))
train <- match(colnames(metabric.s), colnames(temp))
test <- match(rownames(gse76040.s), colnames(temp))

##### ENDORSE signatures -----------------------------

gse76040.ssgsea <- gsva(temp[, test], gset.idx.list = list("ENDORSE"=ENDORSE, 
                                                           "HALLMARK_ESTROGEN_RESPONSE_EARLY" = H_ESTR_EARLY, 
                                                           "HALLMARK_ESTROGEN_RESPONSE_LATE" = H_ESTR_LATE), 
                        method="ssgsea", kcdf="Gaussian")

gse76040.ssgsea2 <- rbind(gsva(temp[, test], Hset, method='ssgsea', kcdf='Gaussian'), 
                           gsva(temp[, test], c6set, method='ssgsea', kcdf='Gaussian'))

df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 t(MBGSE76040.ENDORSE))

cfit.0 = coxph(Surv(Time, Event) ~ ENDORSE, data=df)
cfit.5 = coxph(Surv(Time, Event) ~ HALLMARK_ESTROGEN_RESPONSE_EARLY, data=df)
cfit.6 = coxph(Surv(Time, Event) ~ ENDORSE + HALLMARK_ESTROGEN_RESPONSE_EARLY, data=df)

########### predict with E2 
gse76040.pred <- predict(cfit.6, data.frame(t(gse76040.ssgsea)), type="risk")
gse76040.pred.grp <- MakeRiskGrp(gse76040.pred)
gse76040.pred2.grp <- MakeRiskGrp(gse76040.pred, ngrps = 2)

############# Apply TransCONFIRM signature #########################
gse76040.TransCONFIRM <- TransCONFIRM(t(temp[, test]))

############# Apply SET signature #########################
gse76040.SET <- SET_Score(t(temp[, test]))

df.gse76040 <- data.frame("Ki67"=as.numeric(gse76040.annotV2$`ki67 (%)`), 
                          "RiskGroup"=gse76040.pred.grp, 
                          "RiskGroup2"=gse76040.pred2.grp, 
                          "Risk"=gse76040.pred,
                          "SET"=gse76040.SET,
                          "TransCONFIRM"=gse76040.TransCONFIRM)

summary(glm(Ki67 ~ Risk, data=df.gse76040)) #P=2.53e-05
TukeyHSD(aov(Ki67 ~ RiskGroup, data=df.gse76040)) #Medium Risk-Low Risk P=0.0405270 #High Risk-Low Risk P=0.0001753
summary(aov(Ki67 ~ RiskGroup2, data=df.gse76040)) #P=0.02
summary(glm(Ki67 ~ SET, data=df.gse76040)) #P=0.307
summary(aov(Ki67 ~ TransCONFIRM, data=df.gse76040)) #P=0.0516 

geom_point(size=3) +
  scale_color_continuous(type = 'viridis', name="GES") +
  geom_smooth(method='glm', se=T, color='navy', size=3) +
  labs(x="ENDORSE", y="GES", title=gsub("_", " ", colnames(df)[i])) +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(size=14, face="bold"),
        title = element_text(size=16, face="bold", color="grey25"))

P1 <- ggplot(df.gse76040, aes(x=Risk, y=Ki67, color = Risk))  +
  geom_point(size=3) +
  scale_color_continuous(name="ENDORSE") +
  geom_smooth(method='glm', se=T, color='darkgrey', size=3) +
  labs(x="ENDORSE", y="Ki67 (%)") +
  theme_classic(base_size = 18)

P2 <-  ggboxplot(df.gse76040, x="RiskGroup", y="Ki67", fill="RiskGroup", palette = "jco", add="jitter") +
  labs(x="ENDORSE", y="Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

P3 <- ggplot(df.gse76040, aes(x=SET, y=Ki67, color = SET))  +
  geom_point(size=3) +
  scale_color_continuous(name="SET") +
  geom_smooth(method='glm', se=T, color='darkgrey', size=3) +
  labs(x="SET", y="Ki67 (%)") +
  theme_classic(base_size = 18)

P4 <-  ggboxplot(df.gse76040, x="TransCONFIRM", y="Ki67", fill="TransCONFIRM", palette = "jco", add="jitter") +
  labs(x="TransCONFIRM", y="Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

pdf("./../Figures/Predictions_GSE76040_ENDORSE_V2.pdf", width=6.5, height=6)
print(P1)
print(P2)
print(P3)
print(P4)
dev.off()


#### Scatter plot 
gse76040.gams <- apply(t(gse76040.ssgsea2), 2, GamApp, gse76040.pred)
write.csv(t(gse76040.gams), file="./../Output/gse76040.gams.csv")

head(sort(gse76040.gams[1,]), 20)

df <- data.frame(cbind(t(gse76040.ssgsea2), df.gse76040))  
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
  
PPL[[4]] <-  ggplot(df, aes(x=ATF2_S_UP.V1_UP, y=Risk, color=ATF2_S_UP.V1_UP  ))  +
  geom_point(size=3) +
  scale_colour_gradient2(low = "blue", mid = "grey80", high = "red", name="GES") +
  geom_smooth(method='gam', se=T, color='navy', size=3) +
  labs(x="ATF2 upregulation", y="ENDORSE") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(size=14, face="bold"),
        title = element_text(size=16, face="bold", color="grey25"))

pdf("./../Figures/GSE76040.pathways.pdf", height =  4, width=7.5)
egg::ggarrange(plots=PPL, nrow=2)
dev.off()


######### Compare with the empirical and estrogen response signatures ############
gse76040.pred.c0 <- predict(cfit.0, data.frame(t(gse76040.ssgsea)), type="risk")
gse76040.pred.grp.c0 <- MakeRiskGrp(gse76040.pred.c0)

gse76040.pred.c5 <- predict(cfit.5, data.frame(t(gse76040.ssgsea)), type="risk")
gse76040.pred.grp.c5 <- MakeRiskGrp(gse76040.pred.c5)

df.gse76040.2 <- data.frame("Ki67"=as.numeric(gse76040.annotV2$`ki67 (%)`), 
                          "RiskGroup0"=gse76040.pred.grp.c0, 
                          "Risk0"=gse76040.pred.c0,
                          "RiskGroup5"=gse76040.pred.grp.c5, 
                          "Risk5"=gse76040.pred.c5)

summary(glm(Ki67 ~ Risk0, data=df.gse76040.2)) #7.08e-09
summary(aov(Ki67 ~ RiskGroup0, data=df.gse76040.2)) #1.11e-08

summary(glm(Ki67 ~ Risk5, data=df.gse76040.2)) #0.2
summary(aov(Ki67 ~ RiskGroup5, data=df.gse76040.2)) #0.85

pdf("./../Figures/GSE76040_Emp_Hstr.pdf", height=6.5, width=6)
ggplot(df.gse76040.2, aes(x=Risk0, y=Ki67, color = Risk0))  +
  geom_point(size=3) +
  scale_color_continuous(name="Empirical") +
  geom_smooth(method='glm', se=T, color='darkgrey', size=3) +
  labs(x="Empirical", y="Ki67 (%)") +
  theme_classic(base_size = 18)

ggboxplot(df.gse76040.2, x="RiskGroup0", y="Ki67", fill="RiskGroup0", palette = "jco", add="jitter") +
  labs(x="Empirical", y="Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

ggplot(df.gse76040.2, aes(x=Risk5, y=Ki67, color = Risk5))  +
  geom_point(size=3) +
  scale_color_continuous(name="Estrogen respomnse") +
  geom_smooth(method='glm', se=T, color='darkgrey', size=3) +
  labs(x="Estrogen respomnse", y="Ki67 (%)") +
  theme_classic(base_size = 18)

ggboxplot(df.gse76040.2, x="RiskGroup5", y="Ki67", fill="RiskGroup5", palette = "jco", add="jitter") +
  labs(x="Estrogen response", y="Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())
dev.off()
