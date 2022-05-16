library(caret)
library(glmnet)
library(survival)
library(GSVA)
library(doParallel)

no_cores <- detectCores()  
print(no_cores)
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl)  

load("./../Output/metabric_er.RData")
source("Useful_functions.R")

# Create data splits for 50 repeats 10-fold CV 
testindex <- list()
for (i in 1:50) {
  set.seed(2020 + i)
  testindex[[i]] <- createFolds(metabric.annotV2$OS_EVENT, k = 10)
}

# 50 rep 10-fold CV
MB_CoxNoClin_cv_results <- list()

MB_CoxNoClin_cv_results <- foreach (i=1:50) %dopar% { 
  sink("LOG_METABRIC_noclinCV.txt", append=TRUE)
  cat(paste("Starting iteration", i, "\n")) 
  
  foreach (k=1:10) %do% {  
    cat(paste("Fold", k, "\n"))
    temp <- list()
    
    # Define TRAINing and test sample index 
    test_s <- testindex[[i]][[k]]
    TRAIN_s <- c(1:nrow(metabric.annotV2))[-test_s]
    
    # Survival objects 
    TRAIN.surv = Surv(metabric.annotV2$OS_MONTHS[TRAIN_s], as.numeric(metabric.annotV2$OS_EVENT[TRAIN_s]))
    test.surv = Surv(metabric.annotV2$OS_MONTHS[test_s], as.numeric(metabric.annotV2$OS_EVENT[test_s]))
    
    # Make the signatures 
    TRAIN.cv.exp <- cv.glmnet(x=cbind(t(metabric.exp[, TRAIN_s]), 
                                      metabric.annotV2$AGE_AT_DIAGNOSIS[TRAIN_s]),
                              y=TRAIN.surv, family = 'cox', alpha=1, parallel = T)
    
    TRAIN_gseed <- names(which(getCoefs(TRAIN.cv.exp, option='1se') > 0))
    TRAIN_geneset <- names(GetCorGenes(TRAIN_gseed, t(metabric.exp), cof=0.75))
    
    # Calculate ssGSEA score
    TRAIN.ssgsea <- gsva(metabric.exp[, TRAIN_s], gset.idx.list = list("ENDORSE"=TRAIN_geneset), 
                         verbose=TRUE, kcdf="Gaussian", method='ssgsea')
    
    test.ssgsea <- gsva(metabric.exp[, test_s], gset.idx.list = list("ENDORSE"=TRAIN_geneset), 
                        verbose=TRUE, kcdf="Gaussian", method='ssgsea')
    
    # Create Cox proportional hazards model
    TRAIN.df <- data.frame("Age"=metabric.annotV2$AGE_AT_DIAGNOSIS[TRAIN_s], 
                           "ENDORSE"=TRAIN.ssgsea[1,],
                           "Time"=metabric.annotV2$OS_MONTHS[TRAIN_s], 
                           "Event"=as.numeric(metabric.annotV2$OS_EVENT[TRAIN_s]))
    
    test.df <- data.frame("Age"=metabric.annotV2$AGE_AT_DIAGNOSIS[test_s], 
                          "ENDORSE"=test.ssgsea[1,],
                          "Time"=metabric.annotV2$OS_MONTHS[test_s], 
                          "Event"=as.numeric(metabric.annotV2$OS_EVENT[test_s]))
    
    TRAIN.model <- coxph(Surv(Time, Event) ~ ., data=TRAIN.df)
    
    # Model test datset with TRAINing gene signature
    test.model <- coxph(Surv(Time, Event) ~ ., data=test.df)
    
    # Predict on test datset
    test.actual <- predict(coxph(Surv(Time, Event) ~ ., data=test.df), type="risk")
    test.pred <- predict(TRAIN.model, newdata = test.df, type="risk")
    
    plot(test.actual, test.pred)
    cor(test.actual, test.pred, method='spearman')
    
    # Export objects to save 
    temp[[1]] <- TRAIN_geneset # geneset that can be used to calculate ssGSEA scores and get consensus geneset 
    temp[[2]] <- TRAIN.model # the Cox PH model for TRAINing data
    temp[[3]] <- test.model # the Cox PH model for test data with geneset identified from TRAINing split
    temp[[4]] <- test.actual # actual risk for test data
    temp[[5]] <- test.pred # risk for test data predicted using TRAINing model with test age and ssGSEA score
    
    MB_CoxNoClin_cv_results[[k]] <- temp
  }
}
sink()

save(MB_CoxNoClin_cv_results, file="./../Output/MB_CoxNoClin_cv_results.RData")

