#####################################
########### GSE87411 ################

#Expression 
gse87411 <- read.delim("./../Data/GSE87411/GSE87411_V2.txt", comment.char="!", header=TRUE, fill=TRUE)

#Header annotation 
gse87411.annot <- read.delim("./../Data/GSE87411/GSE87411_V2.txt", header=FALSE, fill=TRUE, nrows=6)
gse87411.annot <- t(gse87411.annot)
colnames(gse87411.annot) <- c("Patient", "BL_KI", "EOT_KI", "ER_Response", "Treatment", "ID_REF")
gse87411.annot <- gse87411.annot[-1,]
gse87411.annot <- data.frame(gse87411.annot)

gse87411.annot$Time <- "EOT"
gse87411.annot$Time[grep("BL", gse87411.annot$Patient, ignore.case=T)] <- "BL"
gse87411.annot$Response <- "sensitive"
gse87411.annot$Response[grep("resistant", gse87411.annot$ER_Response, ignore.case=T)] <- "resistant"

# Extract KI67 and response info 
ki67_BL <- sapply(strsplit(as.character(gse87411.annot$BL_KI), " "), "[", 4)
ki67_BL <- gsub("%", "", ki67_BL)
ki67_BL <- gsub("~", "", ki67_BL)
ki67_BL[grep("-", ki67_BL)] <- c(65, 65)
ki67_BL <- as.numeric(ki67_BL)

ki67_EOT <- sapply(strsplit(as.character(gse87411.annot$EOT_KI), " "), "[", 3)
ki67_EOT <- gsub("%", "", ki67_EOT)
ki67_EOT <- gsub("<", "", ki67_EOT)
ki67_EOT <- gsub(">", "", ki67_EOT)
ki67_EOT <- as.numeric(ki67_EOT)

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

### Match matrices and train the model on batch-corrected data ####-------------------------
library(sva)
temp <- DoCombat(metabric.s, t(gse87411.s))
train <- match(colnames(metabric.s), colnames(temp))
test <- match(rownames(gse87411.s), colnames(temp))

MBGSE87411.ENDORSE <- gsva(temp[, train], gset.idx.list = list("ENDORSE"=ENDORSE, 
                                                                "HALLMARK_ESTROGEN_RESPONSE_EARLY" = H_ESTR_EARLY, 
                                                                "HALLMARK_ESTROGEN_RESPONSE_LATE" = H_ESTR_LATE), 
                            method="ssgsea", kcdf="Gaussian")

df <- data.frame('Time'=metabric.annotV2$OS_MONTHS, 
                 'Event'=as.numeric(metabric.annotV2$OS_EVENT), 
                 'Age'=metabric.annotV2$AGE_AT_DIAGNOSIS,
                 t(MBGSE87411.ENDORSE))

cfit.0 = coxph(Surv(Time, Event) ~ ENDORSE, data=df)
cfit.5 = coxph(Surv(Time, Event) ~ HALLMARK_ESTROGEN_RESPONSE_EARLY, data=df)
cfit.6 = coxph(Surv(Time, Event) ~ ENDORSE + HALLMARK_ESTROGEN_RESPONSE_EARLY, data=df)

##### ENDORSE signatures -----------------------------

gse87411.ssgsea <- gsva(temp[, test], gset.idx.list = list("ENDORSE"=ENDORSE,
                                           "HALLMARK_ESTROGEN_RESPONSE_EARLY" = H_ESTR_EARLY, 
                                           "HALLMARK_ESTROGEN_RESPONSE_LATE" = H_ESTR_LATE), 
                        method="ssgsea", kcdf="Gaussian")

gse87411.ssgsea2 <- rbind(gsva(temp[, test], Hset, method='ssgsea', kcdf='Gaussian'), 
                           gsva(temp[, test], c6set, method='ssgsea', kcdf='Gaussian'))

# predict
gse87411.pred <- predict(cfit.6, data.frame(t(gse87411.ssgsea)), type="risk")
gse87411.pred.grp <- MakeRiskGrp(gse87411.pred)

############# Apply TransCONFIRM signature #########################
gse87411.TransCONFIRM <- TransCONFIRM(t(temp[, test]))

############# Apply SET signature #########################
gse87411.SET <- SET_Score(t(temp[, test]))

df.gse87411 <- data.frame("Time"=gse87411.annot$Time,
                          "Ki67_BL"=ki67_BL,
                          "Ki67_EOT"=ki67_EOT,
                          "RiskGroup"=gse87411.pred.grp, 
                          "Risk"=gse87411.pred, 
                          "Response"=gse87411.annot$Response, 
                          "SET"=gse87411.SET,
                          "TransCONFIRM"=gse87411.TransCONFIRM)

BL <- which(gse87411.annot$Time == "BL")  
EOT <- which(gse87411.annot$Time == "EOT")  

summary(aov(Ki67_BL ~ RiskGroup, data=df.gse87411[BL,]))
TukeyHSD(aov(Ki67_BL ~ RiskGroup, data=df.gse87411[BL,]))
P1 <- ggboxplot(df.gse87411[BL,], x="RiskGroup", y="Ki67_BL", fill="RiskGroup", palette = "jco", add="jitter") +
  labs(x="ENDORSE", y="BL Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

summary(aov(Ki67_EOT ~ RiskGroup, data=df.gse87411[EOT,]))[[1]][[5]]
TukeyHSD(aov(Ki67_EOT ~ RiskGroup, data=df.gse87411[EOT,]))
P2 <- ggboxplot(df.gse87411[EOT,], x="RiskGroup", y="Ki67_EOT", fill="RiskGroup", palette = "jco", add="jitter") +
  labs(x="ENDORSE", y="EOT Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

summary(lm(Ki67_BL ~ Risk, data=df.gse87411[BL,]))[[4]][2,4]
P3 <- ggplot(df.gse87411[BL,], aes(x=Risk, y=Ki67_BL, color = Risk))  +
  geom_point(size=3) +
  scale_color_continuous(name="ENDORSE") +
  geom_smooth(method='glm', se=T, color='darkgrey', size=3) +
  labs(x="ENDORSE", y="BL Ki67 (%)") +
  theme_classic(base_size = 18)

summary(lm(Ki67_EOT ~ Risk, data=df.gse87411[EOT,]))[[4]][2,4]
P4 <- ggplot(df.gse87411[EOT,], aes(x=Risk, y=Ki67_EOT, color = Risk))  +
  geom_point(size=3) +
  scale_color_continuous(name="ENDORSE") +
  geom_smooth(method='glm', se=T, color='darkgrey', size=3) +
  labs(x="ENDORSE", y="EOT Ki67 (%)") +
  theme_classic(base_size = 18)

summary(lm(gse87411.pred[BL] ~ gse87411.annot$Response[BL]))
P5 <- ggboxplot(df.gse87411[BL,], x="Response", y="Risk", fill="Response", palette = "npg", add="jitter") +
  labs(x="Clinical Response", y="ENDORSE") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

summary(lm(Ki67_BL ~ SET, data=df.gse87411[BL,]))
P6 <- ggplot(df.gse87411[BL,], aes(x=SET, y=Ki67_BL, color = Risk))  +
  geom_point(size=3) +
  scale_color_continuous(name="SET") +
  geom_smooth(method='glm', se=T, color='darkgrey', size=3) +
  labs(x="SET", y="BL Ki67 (%)") +
  theme_classic(base_size = 18)

summary(lm(Ki67_EOT ~ SET, data=df.gse87411[EOT,]))
P7 <- ggplot(df.gse87411[EOT,], aes(x=SET, y=Ki67_EOT, color = Risk))  +
  geom_point(size=3) +
  scale_color_continuous(name="SET") +
  geom_smooth(method='glm', se=T, color='darkgrey', size=3) +
  labs(x="SET", y="EOT Ki67 (%)") +
  theme_classic(base_size = 18)

summary(lm(gse87411.SET[BL] ~ gse87411.annot$Response[BL]))
P8 <- ggboxplot(df.gse87411[BL,], x="Response", y="SET", fill="Response", palette = "npg", add="jitter") +
  labs(x="Clinical Response", y="SET") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

summary(aov(Ki67_BL ~ TransCONFIRM, data=df.gse87411[BL,]))
P9 <- ggboxplot(df.gse87411[BL,], x="TransCONFIRM", y="Ki67_BL", fill="TransCONFIRM", palette = "jco", add="jitter") +
  labs(x="TransCONFIRM", y="BL Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

summary(aov(Ki67_EOT ~ TransCONFIRM, data=df.gse87411[EOT,]))
P10 <- ggboxplot(df.gse87411[EOT,], x="TransCONFIRM", y="Ki67_EOT", fill="TransCONFIRM", palette = "jco", add="jitter") +
  labs(x="TransCONFIRM", y="EOT Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

summary(aov(gse87411.TransCONFIRM[BL] ~ gse87411.annot$Response[BL]))
P11 <- ggboxplot(df.gse87411[BL,], x="Response", y="TransCONFIRM", fill="Response", palette = "npg", add="jitter") +
  labs(x="Clinical Response", y="TransCONFIRM") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())


pdf("./../Figures/Predictions_GSE87411_ENDORSE_V2.pdf", width=6.5, height=6)
print(P1); print(P2); print(P3); print(P4); print(P5)
print(P6); print(P7); print(P8)
print(P9); print(P10); print(P11)
dev.off()


########## Scatter plots 
gse87411_BL.gams <- apply(t(gse87411.ssgsea2[, BL]), 2, GamApp, gse87411.pred[BL])
write.csv(t(gse87411_BL.gams), file="./../Output/gse87411_BLgams.csv")

gse87411_EOT.gams <- apply(t(gse87411.ssgsea2[, EOT]), 2, GamApp, gse87411.pred[EOT])
write.csv(t(gse87411_EOT.gams), file="./../Output/gse87411_EOTgams.csv")

tail(sort(gse87411_BL.gams[3,]), 30)
tail(sort(gse87411_EOT.gams[3,]), 30)

df <- data.frame(cbind(t(gse87411.ssgsea2[,BL]), "Risk"=gse87411.pred[BL])) 

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

pdf("./../Figures/GSE87411_BL.pathways.pdf", height =  4, width=15)
egg::ggarrange(plots=PPL, nrow=1)
dev.off()



df <- data.frame(cbind(t(gse87411.ssgsea2[,EOT]), "Risk"=gse87411.pred[EOT])) 

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

pdf("./../Figures/GSE87411_EOT.pathways.pdf", height =  4, width=15)
egg::ggarrange(plots=PPL, nrow=1)
dev.off()



######### Compare with the empirical and estrogen response signatures ############
gse87411.pred.c0 <- predict(cfit.0, data.frame(t(gse87411.ssgsea)), type="risk")
gse87411.pred.grp.c0 <- MakeRiskGrp(gse87411.pred.c0)

gse87411.pred.c5 <- predict(cfit.5, data.frame(t(gse87411.ssgsea)), type="risk")
gse87411.pred.grp.c5 <- MakeRiskGrp(gse87411.pred.c5)

df.gse87411.2 <- data.frame("Ki67_BL"=ki67_BL,
                          "Ki67_EOT"=ki67_EOT,
                          "RiskGroup0"=gse87411.pred.grp.c0, 
                          "Risk0"=gse87411.pred.c0, 
                          "RiskGroup5"=gse87411.pred.grp.c5, 
                          "Risk5"=gse87411.pred.c5, 
                          "Response"=gse87411.annot$Response)

pdf("./../Figures/GSE87411_Emp_Hstr.pdf", height=6.5, width=6)
summary(aov(Ki67_BL ~ RiskGroup0, data=df.gse87411.2[BL,])) #3.8e-11
ggboxplot(df.gse87411.2[BL,], x="RiskGroup0", y="Ki67_BL", fill="RiskGroup0", palette = "jco", add="jitter") +
  labs(x="Empirical", y="BL Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

summary(aov(Ki67_EOT ~ RiskGroup0, data=df.gse87411.2[EOT,])) #<2e-16
ggboxplot(df.gse87411.2[EOT,], x="RiskGroup0", y="Ki67_EOT", fill="RiskGroup0", palette = "jco", add="jitter") +
  labs(x="Empirical", y="EOT Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())


summary(aov(Ki67_BL ~ RiskGroup5, data=df.gse87411.2[BL,])) #0.09
ggboxplot(df.gse87411.2[BL,], x="RiskGroup5", y="Ki67_BL", fill="RiskGroup5", palette = "jco", add="jitter") +
  labs(x="Estrogen response", y="BL Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

summary(aov(Ki67_EOT ~ RiskGroup5, data=df.gse87411.2[EOT,])) #0.88
ggboxplot(df.gse87411.2[EOT,], x="RiskGroup5", y="Ki67_EOT", fill="RiskGroup5", palette = "jco", add="jitter") +
  labs(x="Estrogen response", y="EOT Ki67 (%)") +
  geom_hline(yintercept=10, colour="grey25", linetype="dotted") +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

dev.off()
