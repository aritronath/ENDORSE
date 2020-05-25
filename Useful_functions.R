
#### Extract glmnet coefficients from cv.glmnet object ---------------------
getCoefs <- function(cv.obj, option=c("min","1se")) {
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
MakeRiskGrp <- function (risk.pred, thresh=c(1,2)) {
  risk.grps <- array(length(risk.pred))
  risk.grps[which(risk.pred <= thresh[1])] <- "Low Risk"
  risk.grps[which(risk.pred > thresh[1] & risk.pred < thresh[2])] <- "Medium Risk"
  risk.grps[which(risk.pred >= thresh[2])] <- "High Risk"
  risk.grps <- factor(risk.grps, levels=c("Low Risk", "Medium Risk", "High Risk"))
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
