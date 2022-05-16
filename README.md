# ENDORSE v1.0: Endocrine resistance signature model for estrogen receptor-positive (ER+) breast cancers

## The ENDORSE model predicts risk of death on endocrine therapy using gene expression profiles from tumors. 

#### This repository contains scripts for developing the ENDORSE model using METABRIC training data. Additional scripts for independent validation of the model are also included. Follow the instructions below to apply ENDORSE to calculate risk: 


[Click here to view and download models and training files](https://github.com/aritronath/ENDORSE/releases/tag/v1.0)


## Using ENDORSE:  

#### 1. Load required libraries and functions 
```{r}
library(GSVA)
library(survival)
library(sva)

##### Make risk groups from predicted risk scores based on specificed threshold ----------------------
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

#### Match matrices by rownames and perform batch correction ---------------------
DoCombat <- function(mat1, mat2) {
  
  comG <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[comG, ]
  mat2 <- mat2[comG, ]
  
  combined_mat <- cbind(mat1, mat2)
  batch_labels <- c(rep("A", ncol(mat1)), rep("B", ncol(mat2)))
  
  corrected_mat <- ComBat(combined_mat, batch = batch_labels)
  
  return(corrected_mat)
}

```

#### 2. Integrate sample gene expression matrix with training data and perform batch correction  
```{r}
# load required training data and model files 
metabric.s <- readRDS("metabric.s.RDS")
metabric.os <- readRDS("metabric.os.RDS")
metabric.os.event <- readRDS("metabric.os.event.RDS")
ENDORSE <- readRDS("ENDORSE.RDS")
H_ESTR_EARLY <- readRDS("H_ESTR_EARLY.RDS")

# perform batch correction
temp <- DoCombat(metabric.s, t(NEWDAT))
train <- match(colnames(metabric.s), colnames(temp))
test <- match(rownames(NEWDAT), colnames(temp))
```

#### 3. Obtain ENDORSE model parameters fron training data 
```{r}
MBS.ENDORSE <- gsva(temp[, train], 
                    gset.idx.list = list("ENDORSE"=ENDORSE, "HALLMARK_ESTROGEN_RESPONSE_EARLY" = H_ESTR_EARLY), 
                    method="ssgsea", kcdf="Gaussian")

df <- data.frame('Time'=metabric.os, 
                 'Event'=as.numeric(metabric.os.event), 
                 t(MBS.ENDORSE))

cfit.6 = coxph(Surv(Time, Event) ~ ENDORSE + HALLMARK_ESTROGEN_RESPONSE_EARLY, data=df)
```

#### 3. Apply ENDORSE model to sample gene expression data 
```{r}
NEWDAT.ssgsea <- gsva(temp[, test], gset.idx.list = list("ENDORSE"=ENDORSE,
                                                           "HALLMARK_ESTROGEN_RESPONSE_EARLY" = H_ESTR_EARLY), 
                        method="ssgsea", kcdf="Gaussian")

# predicted continuous risk score
NEWDAT.pred <- predict(cfit.6, data.frame(t(NEWDAT.ssgsea)), type="risk")

# predicted risk groups (low, medium, high risk)
NEWDAT.pred.grp <- MakeRiskGrp(NEWDAT.pred)

```
