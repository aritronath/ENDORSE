
#### Extract glmnet coefficients from cv.glmnet object ---------------------
getCoefs <- function(cv.obj, option=c("1se")) {
  if (option == "min") lamb <- cv.obj$lambda.min
  if (option == "1se") lamb <- cv.obj$lambda.1se
  lN <- match(lamb, cv.obj$glmnet.fit$lambda)
  TheCoefs <- cv.obj$glmnet.fit$beta[, lN]
  TheCoefs <- TheCoefs[-which(TheCoefs == 0 | names(TheCoefs) == "")]
  return(TheCoefs)
}

### Get genes highly correlated with selected genes -------------------------
GetCorGenes <- function(gs, exp, cof=0.75) { 
  C_UP <- NULL
  
  gs <- intersect(gs, colnames(exp)) # some co-variates might not be genes
  
  if (length(gs) == 1) {
    CORS <- apply(exp, 2, function (x) cor(x, exp[, gs]))
    C_UP <- append(C_UP, which(CORS > cof))
  }
  
  if (length(gs) > 1) {
    pb = txtProgressBar(min = 1, max = length(gs), initial = 1) 
    for (i in 1:length(gs)) {
      CORS <- apply(exp, 2, function (x) cor(x, exp[, gs[i]]))
      C_UP <- append(C_UP, which(CORS > cof))
      setTxtProgressBar(pb,i)
    }
  }
  
  return(C_UP)
}

##### Make risk groups from predicted risk ----------------------
MakeRiskGrp <- function (risk.pred, ngrps=3, thresh=c(1,2)) {
  risk.grps <- array(length(risk.pred))
  
  if (ngrps == 2) {
    risk.grps[which(risk.pred <= 0.5*(thresh[1] + thresh[2]))] <- "Low Risk"
    risk.grps[which(risk.pred > 0.5*(thresh[1] + thresh[2]))] <- "High Risk"
    risk.grps <- factor(risk.grps, levels=c("Low Risk", "High Risk"))
  }
  
  if (ngrps == 3) {
  risk.grps[which(risk.pred <= thresh[1])] <- "Low Risk"
  risk.grps[which(risk.pred > thresh[1] & risk.pred < thresh[2])] <- "Medium Risk"
  risk.grps[which(risk.pred >= thresh[2])] <- "High Risk"
  risk.grps <- factor(risk.grps, levels=c("Low Risk", "Medium Risk", "High Risk"))
  }
  
  return(risk.grps)
}

###### Rescale values between 0 and 1 -------------------------
ReScale <- function (V) {
  maxV = max(V)
  minV = min(V)
  Vs   = (V - minV) / (maxV - minV)
  return(Vs)
}

#### Export a Cox model to csv file -----------------------
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

#### Match matrices by rownames and perform ComBat batch correction ---------------------
DoCombat <- function(mat1, mat2) {
  
  comG <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[comG, ]
  mat2 <- mat2[comG, ]
  
  combined_mat <- cbind(mat1, mat2)
  batch_labels <- c(rep("A", ncol(mat1)), rep("B", ncol(mat2)))
  
  corrected_mat <- ComBat(combined_mat, batch = batch_labels)
  
  return(corrected_mat)
}


##### SET ER/PR Score ############# ------------------------------------
SET_Score <- function(mat) {
  
  # mat: genes should be in columns, samples in rows 
  
  T_genes <- c("SLC39A6", "STC2", "CA12", "ESR1", "PDZK1", "NPY1R", "CD2", "MAPT", "QDPR",
               "AZGP1", "ABAT", "ADCY1", "CD3D", "NAT1", "MRPS30", "DNAJC12", "SCUBE2", "KCNE4")
  
  R_genes <- c("LDHA", "ATP5J2", "VDAC2", "DARS", "UGP2", "UBE2Z", "AK2", "WIPF2", "APPBP2", "TRIM2")
  
  Mean_T <- apply(mat[, intersect(T_genes, colnames(mat))], 1, mean)
  Mean_R <- apply(mat[, intersect(R_genes, colnames(mat))], 1, mean)
  
  SET <- Mean_T - Mean_R + 2
  
  return(SET)
}

##### TransCONFIRM gene set clusters ############# ------------------------------------
TransCONFIRM <- function(mat) {
  
  # mat: genes should be in columns, samples in rows
  
  TC_genes <- c("CPNE1", "TFAP2C", "NOTCH2", "CTTNBP2NL", "PLA2G10", "NCLN", "SLC25A41", "USP5", "TANC1", "CSDA", 
                "ARRDC3", "KIAA0913", "C2orf68", "GLTP", "SP1", "ILF2", "C10orf137", "TTPAL", "ADD3", "BMPR1A",
                "ATP11B", "PPP2R3A", "COQ2", "TRIQK", "DEFB104A", "KRTAP13-3", "ZNF483", "TCEB3", "SERINC5",
                "MBD4", "SEC23A", "CCDC93", "DIPAH3-AS2", "LIN28B", "EGFLAM-AS", "OR4K5", "TMEM81")
  
  TC_dist <- dist(mat[, intersect(TC_genes, colnames(mat))])
  
  TC_cuts <- cutree(hclust(TC_dist), k = 2)
  
  return(TC_cuts)
}


##### Oncotype Dx scores  ############# ------------------------------------
ODX_RS <- function (mat) {
  # Define list of genes
  Ref_grp <- c("ACTB", "GAPDH", "RPLP0", "GUSB", "TFRC")
  ER_grp <- c("ESR1", "PGR", "BCL2", "SCUBE2")
  HER2_grp <- c("GRB7", "ERBB2")
  Prol_grp <- c("MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2")
  #Inv_grp <- c("MMP11", "CTSV")
  Inv_grp <- c("MMP11", "CTSL2")
  
  All_grp <- c(Ref_grp, ER_grp, HER2_grp, Prol_grp, Inv_grp, "GSTM1", "CD68") #,"BAG1")

  # Filter for 21 genes  
  mat <- data.frame(mat[, match(All_grp, colnames(mat))])
  
  # Rescale expression values to lie in the range of 1:20
  #mat <- (((mat - min(mat)) / (max(mat) - min(mat))) * 30) + 20
  
  # Normalize by reference genes 
  mat <- mat - rowMeans(mat[, Ref_grp])
  
  # Calculate group scores 
  Her2_score <- 0.9*mat$"GRB7" + 0.1*mat$"ERBB2"
  Her2_score[which(Her2_score < 8)]  <- 8 
    
  ER_score <- (0.8*mat$"ESR1" + 1.2*(mat$"PGR" + mat$"BCL2" + mat$"SCUBE2")) / 4
    
  Prol_score <- rowMeans(mat[ ,Prol_grp])
  Prol_score[which(Prol_score < 6.5)] <- 6.5
    
  Inv_score <- rowMeans(mat[, Inv_grp])
    
  # Calculate unscaled recurrence scores
  RSu <- 0.47*Her2_score + -0.34*ER_score + 1.04*Prol_score + 0.1*Inv_score + 0.05*mat$"CD68" + -0.08*mat$"GSTM1" #+ -0.07*mat$"BAG1"
    
  # Calculate scaled recurrence scores
  #if (RSu > 100) RS <- 100
  #if (RSu < 0) RS <- 0
  #if (RSu > 0 && RSu < 100) RS <- 20*(RSu - 6.7)
  return(RSu)
}


##### Functions for returning results of generalized additve models ############# ------------------------------------
GamApp <- function(X, Y) {
  temp <- data.frame("X"=X, "Y"=Y)
  gmod <- gam(Y ~ s(X), data=temp)
  gmod_sum <- summary(gmod)
  return(c(gmod_sum$s.pv, gmod_sum$r.sq,  gmod_sum$dev.expl))
}

GamApp2 <- function(X1, X2, Y) {
  temp <- data.frame("X1"=X1, "X2"=X2, "Y"=Y)
  gmod <- gam(Y ~ s(X1) + s(X2), data=temp)
  gmod_sum <- summary(gmod)
  return(c(gmod_sum$s.pv, gmod_sum$r.sq,  gmod_sum$dev.expl))
}



