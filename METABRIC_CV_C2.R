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
load("./../Output/metabric.ssGSEA.RData")
source("Useful_functions.R")

# Create data splits for 50 repeats 10-fold CV 
testindex <- list()
for (i in 1:50) {
  set.seed(2020 + i)
  testindex[[i]] <- createFolds(metabric.annotV2$OS_EVENT, k = 10)
}

# 50 rep 10-fold CV
MB_C2_Model_cv_results <- list()

MB_C2_Model_cv_results <- foreach (i=1:50) %dopar% { 
  sink("LOG_MB_C2_CV.txt", append=TRUE)
  cat(paste("Starting iteration", i, "\n")) 
  
  foreach (k=1:10, .errorhandling='pass') %do% {  
    cat(paste("Fold", k, "\n"))
    temp <- list()
    
    # Define TRAINing and test sample index 
    test_s <- testindex[[i]][[k]]
    TRAIN_s <- c(1:nrow(metabric.annotV2))[-test_s]
    
    # Survival objects 
    TRAIN.surv = Surv(metabric.annotV2$OS_MONTHS[TRAIN_s], as.numeric(metabric.annotV2$OS_EVENT[TRAIN_s]))
    test.surv = Surv(metabric.annotV2$OS_MONTHS[test_s], as.numeric(metabric.annotV2$OS_EVENT[test_s]))
    
    # Make the signatures 
    TRAIN.cv.c2 <- cv.glmnet(x=cbind(t(metabric.c2[, TRAIN_s]), 
                                     metabric.annotV2$AGE_AT_DIAGNOSIS[TRAIN_s]), 
                             y=TRAIN.surv, family = 'cox', alpha=1, parallel = T)
    
    
    n.mbcv <- names(getCoefs(TRAIN.cv.c2, "min"))
    
    # Create Cox proportional hazards model
    TRAIN.df <- data.frame("Age"=metabric.annotV2$AGE_AT_DIAGNOSIS[TRAIN_s], 
                           t(metabric.c2[n.mbcv, TRAIN_s]),
                           "Time"=metabric.annotV2$OS_MONTHS[TRAIN_s], 
                           "Event"=as.numeric(metabric.annotV2$OS_EVENT[TRAIN_s]))
    
    test.df <- data.frame("Age"=metabric.annotV2$AGE_AT_DIAGNOSIS[test_s], 
                          t(metabric.c2[n.mbcv, test_s]),
                          "Time"=metabric.annotV2$OS_MONTHS[test_s], 
                          "Event"=as.numeric(metabric.annotV2$OS_EVENT[test_s]))
    
    # Model TRAINing datset with TRAINing gene signature
    TRAIN.model <- coxph(Surv(Time, Event) ~ ., data=TRAIN.df)
    
    # Model test datset with TRAINing gene signature
    test.model <- coxph(Surv(Time, Event) ~ ., data=test.df)
    
    # Predict on test datset
    test.actual <- predict(coxph(Surv(Time, Event) ~ ., data=test.df), type='risk')
    test.pred <- predict(TRAIN.model, newdata = test.df, type='risk')
    
    # Export objects to save 
    temp[[1]] <- n.mbcv # geneset that can be used to calculate ssGSEA scores and get consensus geneset 
    temp[[2]] <- TRAIN.model # the Cox PH model for TRAINing data
    temp[[3]] <- test.model # the Cox PH model for test data with geneset identified from TRAINing split
    temp[[4]] <- test.actual # actual risk for test data
    temp[[5]] <- test.pred # risk for test data predicted using TRAINing model with test age and ssGSEA score
    
    return(temp)
  }
}
sink()
save(MB_C2_Model_cv_results, file="./../Output/MB_C2_Model_cv_results.RData")

