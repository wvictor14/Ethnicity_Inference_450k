---
title: "05_Deploy_classifier"
author: "Victor Yuan"
date: "November 14, 2018"
output:
  html_document:
    keep_md: true
    toc: true
    toc_float: true
    toc_depth: 2
editor_options: 
  chunk_output_type: console
---

This script is for extracting the final model and saving it for publication. To make this tool more 
user-friendly, I create a function that wraps the prediction to output the following variables:

* Probabilities that a sample belongs to a specific class (1 for each ethnicity)
* A class label determined by the highest class-specific probability
* A class label determined after applying a user-supplied threshold function for 'ambiguous' samples
* The highest class-specific probability, used to determine the threshold

# Libraries and data


```r
library(dplyr)
```

```
## Warning: package 'dplyr' was built under R version 3.5.1
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(impute)
library(caret)
```

```
## Loading required package: lattice
```

```
## Loading required package: ggplot2
```

```
## Warning: package 'ggplot2' was built under R version 3.5.1
```

```r
library(glmnet)
```

```
## Loading required package: Matrix
```

```
## Loading required package: foreach
```

```
## Loaded glmnet 2.0-16
```

```r
library(broom)
```

```
## Warning: package 'broom' was built under R version 3.5.1
```

```r
# Load model and data
GLM_cv <- readRDS('../../Robjects_final/02_GLM_cv_logloss.rds')
betas <- readRDS('../../Robjects_final/01_processed_betas_EPIC.rds') 
dim(betas) #  319625    510
```

```
## [1] 319625    510
```

```r
pDat <- readRDS('../../Robjects_final/01_pDat.rds')
dim(pDat) # 510 22
```

```
## [1] 510  22
```

```r
# knn impute
sum(is.na(betas))
```

```
## [1] 340
```

```r
set.seed(1)
betas <- impute.knn(as.matrix(betas), maxp = 15000)$data
```

```
## Cluster size 319625 broken into 163356 156269 
## Cluster size 163356 broken into 106218 57138 
## Cluster size 106218 broken into 64684 41534 
## Cluster size 64684 broken into 35473 29211 
## Cluster size 35473 broken into 32576 2897 
## Cluster size 32576 broken into 17376 15200 
## Cluster size 17376 broken into 13853 3523 
## Done cluster 13853 
## Done cluster 3523 
## Done cluster 17376 
## Cluster size 15200 broken into 9274 5926 
## Done cluster 9274 
## Done cluster 5926 
## Done cluster 15200 
## Done cluster 32576 
## Done cluster 2897 
## Done cluster 35473 
## Cluster size 29211 broken into 13385 15826 
## Done cluster 13385 
## Cluster size 15826 broken into 12305 3521 
## Done cluster 12305 
## Done cluster 3521 
## Done cluster 15826 
## Done cluster 29211 
## Done cluster 64684 
## Cluster size 41534 broken into 21426 20108 
## Cluster size 21426 broken into 9507 11919 
## Done cluster 9507 
## Done cluster 11919 
## Done cluster 21426 
## Cluster size 20108 broken into 8321 11787 
## Done cluster 8321 
## Done cluster 11787 
## Done cluster 20108 
## Done cluster 41534 
## Done cluster 106218 
## Cluster size 57138 broken into 28359 28779 
## Cluster size 28359 broken into 16239 12120 
## Cluster size 16239 broken into 11002 5237 
## Done cluster 11002 
## Done cluster 5237 
## Done cluster 16239 
## Done cluster 12120 
## Done cluster 28359 
## Cluster size 28779 broken into 15560 13219 
## Cluster size 15560 broken into 7090 8470 
## Done cluster 7090 
## Done cluster 8470 
## Done cluster 15560 
## Done cluster 13219 
## Done cluster 28779 
## Done cluster 57138 
## Done cluster 163356 
## Cluster size 156269 broken into 63395 92874 
## Cluster size 63395 broken into 32426 30969 
## Cluster size 32426 broken into 15319 17107 
## Cluster size 15319 broken into 4468 10851 
## Done cluster 4468 
## Done cluster 10851 
## Done cluster 15319 
## Cluster size 17107 broken into 10970 6137 
## Done cluster 10970 
## Done cluster 6137 
## Done cluster 17107 
## Done cluster 32426 
## Cluster size 30969 broken into 15437 15532 
## Cluster size 15437 broken into 2943 12494 
## Done cluster 2943 
## Done cluster 12494 
## Done cluster 15437 
## Cluster size 15532 broken into 6949 8583 
## Done cluster 6949 
## Done cluster 8583 
## Done cluster 15532 
## Done cluster 30969 
## Done cluster 63395 
## Cluster size 92874 broken into 65459 27415 
## Cluster size 65459 broken into 43091 22368 
## Cluster size 43091 broken into 17456 25635 
## Cluster size 17456 broken into 9317 8139 
## Done cluster 9317 
## Done cluster 8139 
## Done cluster 17456 
## Cluster size 25635 broken into 14548 11087 
## Done cluster 14548 
## Done cluster 11087 
## Done cluster 25635 
## Done cluster 43091 
## Cluster size 22368 broken into 21552 816 
## Cluster size 21552 broken into 9176 12376 
## Done cluster 9176 
## Done cluster 12376 
## Done cluster 21552 
## Done cluster 816 
## Done cluster 22368 
## Done cluster 65459 
## Cluster size 27415 broken into 12321 15094 
## Done cluster 12321 
## Cluster size 15094 broken into 7648 7446 
## Done cluster 7648 
## Done cluster 7446 
## Done cluster 15094 
## Done cluster 27415 
## Done cluster 92874 
## Done cluster 156269
```

```r
#subset out south asians
betas_SA <- betas[,which(pDat$Ethnicity == 'South_Asian')]
pDat_SA <- pDat[which(pDat$Ethnicity == 'South_Asian'),]
dim(betas_SA);dim(pDat_SA) # 7 samples
```

```
## [1] 319625      7
```

```
## [1]  7 22
```

```r
glm_fit <- GLM_cv$finalModel
features <- glm_fit$beta$African@Dimnames[[1]]
```

Here I manually create the outputs mentioned above. Later I implement a function that automatically
does this.

# Infer ethnicity 


```r
# obtain predictions
preds <- as.data.frame(predict(glm_fit, t(betas_SA), s = glm_fit$lambdaOpt, type = 'class'))
probs <- as.data.frame(predict(glm_fit, t(betas_SA), s = glm_fit$lambdaOpt, type = 'response'))

#combine
pred_prob <- cbind(preds, probs)
colnames(pred_prob) <- c('Predicted_ethnicity_nothresh', paste0('Prob_', glm_fit$classnames))

# calculate highest probability
pred_prob$Highest_Prob <- apply(pred_prob[,2:4], 1, max)

# call ambiguous if below threshold
pred_prob$Predicted_ethnicity <- ifelse(pred_prob$Highest_Prob < 0.75, 'Ambiguous', 
                                        as.character(pred_prob$Predicted_ethnicity_nothresh))
pred_prob$Sample_ID <- rownames(pred_prob)
pred_prob <- pred_prob[,c(7,1,6, 2:5)]
```

# Function: Attempt 1

Now I wrap the above code into a function:


```r
pl_infer_ethnicity <- function(betas, threshold = 0.75){ # betas need to be in the form of samples in columns
  if(!all(features %in% rownames(betas))) {
    stop('Rownames of betas df must include all 319625 features used to fit classifier.')
  }
  
  betas <- t(betas[features,])
  
  # obtain predictions
  preds <- as.data.frame(predict(glm_fit, betas, s = glm_fit$lambdaOpt, type = 'class'))
  probs <- as.data.frame(predict(glm_fit, betas, s = glm_fit$lambdaOpt, type = 'response'))
  p <- cbind(preds, probs)
  colnames(p) <- c('Predicted_ethnicity_nothresh', paste0('Prob_', glm_fit$classnames))

  p$Highest_Prob <- apply(p[,2:4], 1, max)
  
  p$Predicted_ethnicity <- ifelse(p$Highest_Prob < threshold, 'Ambiguous', 
                                        as.character(p$Predicted_ethnicity_nothresh))
  p$Sample_ID <- rownames(p)
  p <- p[,c(7,1,6, 2:5)]
  
  return(p)
}

pl_infer_ethnicity(betas_SA, threshold = 0.75)
```

```
##             Sample_ID Predicted_ethnicity_nothresh Predicted_ethnicity
## PM263           PM263                    Caucasian           Caucasian
## PM272           PM272                    Caucasian           Ambiguous
## PM29             PM29                    Caucasian           Ambiguous
## COX_6987     COX_6987                    Caucasian           Caucasian
## COX_7646     COX_7646                        Asian           Ambiguous
## GSM1947112 GSM1947112                    Caucasian           Caucasian
## GSM1947297 GSM1947297                      African             African
##            Prob_African Prob_Asian Prob_Caucasian Highest_Prob
## PM263        0.01433442 0.09107262      0.8945930    0.8945930
## PM272        0.03455324 0.23050156      0.7349452    0.7349452
## PM29         0.01881591 0.26736345      0.7138206    0.7138206
## COX_6987     0.01740569 0.09146101      0.8911333    0.8911333
## COX_7646     0.03549401 0.66372451      0.3007815    0.6637245
## GSM1947112   0.02761343 0.01941917      0.9529674    0.9529674
## GSM1947297   0.82631224 0.01353996      0.1601478    0.8263122
```

Unfortunately the above code requires that new samples have all the features used for training 
(n=319625), when only 1862 are required for the final prediction. Below I extract the coefficients
from the final model and see if I can create the same predictions using a cross product with the 
sample vector.

# Function: Attempt 2


```r
# extract coefficients
coef <- coef(glm_fit, glm_fit$lambdaOpt)

# bind feature coefficients for each cpg and intercept
out <- do.call("cbind", lapply(coef, function(x) x[,1])) 
out <- as.data.frame(out) %>% mutate(feature = rownames(out)) %>% as_tibble()
out$feature[1] <- 'Intercept'

# should be 1862 + 1 intercept
out <- out %>% filter(abs(African) > 0 | abs(Asian) > 0 | abs(Caucasian) > 0)
out %>% arrange(desc(Asian), desc(African), desc(Caucasian))

#filter data to features
newDat <- betas_SA[out$feature[2:nrow(out)],2]

#cross product, adding 1 for intercept term
af <- out$African %*% c(1, newDat)
as <- out$Asian %*% c(1, newDat)
ca <- out$Caucasian %*% c(1, newDat)

af/sum(af,as,ca)
as/sum(af,as,ca)
ca/sum(af,as,ca)
```

After taking out the coefficients and trying the cross product, I get values that I can't make sense
of. I think I need to implement a loglink function on this output, but I'm not sure how to do this.

Instead, I create a workaround, where given a new samples with the final 1862 features, I add on
'fake' data of the remaining 319625-1862 features so that the predict() function accepts it.

# Function: Attempt 3


```r
# create 'new' data of only the necessary features
pl_features <- predictors(GLM_cv)
newDat <- betas_SA[pl_features,]
dim(newDat) #1862 features, 7 samples
```

```
## [1] 1862    7
```

```r
pl_infer_ethnicity <- function(betas, threshold = 0.75) {
  
  # find all final predictors in new data
  present_features <- intersect(rownames(betas), pl_features)
  print(paste(length(present_features), 'of 1862 predictors present.'))
  
  dat <- betas[present_features,]
  dim(dat)
  
  # find missing non-predictors used for training
  train_features <- glm_fit$beta$African@Dimnames[[1]]
  missing_features <- setdiff(train_features, present_features)
  
  # Create a matrix of zeros for these missing features
  zeros <- data.frame(
    matrix(data = 0, nrow = length(missing_features), ncol = ncol(betas), 
           dimnames = list(missing_features, colnames(betas))))
  
  # add back into the data
  dat_in <- t(rbind(dat, zeros)[train_features,])
  
  # run prediction
  preds <- as.data.frame(predict(glm_fit, dat_in, s = glm_fit$lambdaOpt, type = 'class'))
  probs <- as.data.frame(predict(glm_fit, dat_in, s = glm_fit$lambdaOpt, type = 'response'))
  p <- cbind(preds, probs)
  colnames(p) <- c('Predicted_ethnicity_nothresh', paste0('Prob_', glm_fit$classnames))
  p$Highest_Prob <- apply(p[,2:4], 1, max)
  p$Predicted_ethnicity <- ifelse(p$Highest_Prob < threshold, 'Ambiguous', 
                                        as.character(p$Predicted_ethnicity_nothresh))
  p$Sample_ID <- rownames(p)
  p <- p[,c(7,1,6, 2:5)]
  
  return(p)
}

#now compare:
pl_infer_ethnicity(newDat)
```

```
## [1] "1862 of 1862 predictors present."
```

```
##             Sample_ID Predicted_ethnicity_nothresh Predicted_ethnicity
## PM263           PM263                    Caucasian           Caucasian
## PM272           PM272                    Caucasian           Ambiguous
## PM29             PM29                    Caucasian           Ambiguous
## COX_6987     COX_6987                    Caucasian           Caucasian
## COX_7646     COX_7646                        Asian           Ambiguous
## GSM1947112 GSM1947112                    Caucasian           Caucasian
## GSM1947297 GSM1947297                      African             African
##            Prob_African Prob_Asian Prob_Caucasian Highest_Prob
## PM263        0.01433442 0.09107262      0.8945930    0.8945930
## PM272        0.03455324 0.23050156      0.7349452    0.7349452
## PM29         0.01881591 0.26736345      0.7138206    0.7138206
## COX_6987     0.01740569 0.09146101      0.8911333    0.8911333
## COX_7646     0.03549401 0.66372451      0.3007815    0.6637245
## GSM1947112   0.02761343 0.01941917      0.9529674    0.9529674
## GSM1947297   0.82631224 0.01353996      0.1601478    0.8263122
```

```r
pred_prob
```

```
##             Sample_ID Predicted_ethnicity_nothresh Predicted_ethnicity
## PM263           PM263                    Caucasian           Caucasian
## PM272           PM272                    Caucasian           Ambiguous
## PM29             PM29                    Caucasian           Ambiguous
## COX_6987     COX_6987                    Caucasian           Caucasian
## COX_7646     COX_7646                        Asian           Ambiguous
## GSM1947112 GSM1947112                    Caucasian           Caucasian
## GSM1947297 GSM1947297                      African             African
##            Prob_African Prob_Asian Prob_Caucasian Highest_Prob
## PM263        0.01433442 0.09107262      0.8945930    0.8945930
## PM272        0.03455324 0.23050156      0.7349452    0.7349452
## PM29         0.01881591 0.26736345      0.7138206    0.7138206
## COX_6987     0.01740569 0.09146101      0.8911333    0.8911333
## COX_7646     0.03549401 0.66372451      0.3007815    0.6637245
## GSM1947112   0.02761343 0.01941917      0.9529674    0.9529674
## GSM1947297   0.82631224 0.01353996      0.1601478    0.8263122
```

yay

# Save data

saveRDS(glm_fit, '../../Robjects_final/05_glm_fit.rds')

For south asian samples amy is using:


```r
# add predictions to pDat
pDat_SA <- pDat_SA %>% left_join(pred_prob, by = c('sampleNames' = 'Sample_ID'))

# combine with other samples pDat
pDat_final <- readRDS('../../Robjects_final/03_final_pData.rds')
pDat_save <- pDat_final %>% 
  bind_rows(pDat_SA %>% select(-matAge, -CH, -GA_2) %>% mutate(GA = as.numeric(GA)) %>%
              rename(Self_reported_ethnicity = Ethnicity, Cohort_owner = Dataset))
```

```
## Warning in bind_rows_(x, .id): binding factor and character vector,
## coercing into character vector
```

```
## Warning in bind_rows_(x, .id): binding character and factor vector,
## coercing into character vector
```

```r
dim(pDat_save) # 506 33
```

```
## [1] 506  33
```

```r
# recall ambiguous samples based on 0.75 threshold
pDat_save$Predicted_ethnicity <- ifelse(pDat_save$Highest_Prob < 0.75, 'Ambiguous', 
                                        as.character(pDat_save$Predicted_ethnicity))
table(pDat_save$Predicted_ethnicity)
```

```
## 
##   African Ambiguous     Asian Caucasian 
##        50        33        35       388
```

```r
saveRDS(pDat_save, '../../Robjects_final/05_pDat_506samps.rds')

ggplot(pDat_save, aes(x = Prob_Caucasian, y = Prob_African, col = Predicted_ethnicity)) +
  geom_point()
```

![](05_Deploy_classifier_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
ggplot(pDat_save, aes(x = Prob_Caucasian, y = Prob_African, col = Self_reported_ethnicity)) +
  geom_point()
```

![](05_Deploy_classifier_files/figure-html/unnamed-chunk-6-2.png)<!-- -->
