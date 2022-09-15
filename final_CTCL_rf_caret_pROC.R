
#### ---InstallPackages----####
install.packages("randomForest")
install.packages("scales")

#load libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(magrittr)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(pROC)


####-----Visit 1- Lesional-----####
#lesional 

# How many OTUs do we currently have? 
ntaxa(ps_v1_les.nbUVB.rarefied)



#Random forests can handle sparse matrices, but we still want to prune out lots of our rare OTUs which are just contributing noise. We will do this by eliminating OTUs with an average relative abundance below 0.0001

# Set minlib
minlib = 10000

prunescale = 0.0001


# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean <- taxa_sums(ps_v1_les.nbUVB.rarefied)/nsamples(ps_v1_les.nbUVB.rarefied)
sites.prune <- prune_taxa(tax.mean > prunescale*minlib, ps_v1_les.nbUVB.rarefied)

sites.prune

# Create the training and test datasets

dataMatrix <- data.frame(response = sample_data(sites.prune)$Response, otu_table(sites.prune))
View (dataMatrix)


###
library(caret)
set.seed(3456)

#Split data into training and test set- set to 70: 30
trainIndex <-  caret::createDataPartition(dataMatrix$response, p=0.70, list=FALSE)

# Step 2: Create the training  dataset

Train <- dataMatrix[ trainIndex,]
Test <- dataMatrix[-trainIndex,]

library (randomForest)
#set seed
set.seed(10)

# Performing Random Forest- create model based on train set - set number of tree as 1000

model = randomForest(as.factor(response)~.,data = Train,ntree = 1500,importance = T)
print(model)

#make predictions 

pred <- predict(model, Test, type = "prob", probability =TRUE)[,2]
preds <- prediction(pred, Test$response)

####---plots----####

#make ROC curve
library (ROCR)
#option1 
test_roc = roc(Test$response, pred, plot = TRUE, print.auc = TRUE)
#option 2
pROC_obj <- roc(Test$response, pred,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

# Make a data frame with predictor names and their importance
imp <- importance(model)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 10 predictors
imp.20 <- imp.sort[1:30, ]


# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying pre-nbUVB lesional samples\n into responders or non-responders")

####----Visit 1-Nonles----
ntaxa(ps_v1_non.nbUVB.rarefied)

# Set minlib
minlib = 10000 #adjust from 15000

prunescale = 0.0001


# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean <- taxa_sums(ps_v1_non.nbUVB.rarefied)/nsamples(ps_v1_non.nbUVB.rarefied)
sites.prune <- prune_taxa(tax.mean > prunescale*minlib, ps_v1_non.nbUVB.rarefied)

# Create the training and test datasets

dataMatrix <- data.frame(response = sample_data(sites.prune)$Response, otu_table(sites.prune))
View (dataMatrix)


###
library(caret)
set.seed(3456)

#Split data into training and test set- set to 70: 30
trainIndex <-  caret::createDataPartition(dataMatrix$response, p=0.70, list=FALSE)

# Step 2: Create the training  dataset

Train <- dataMatrix[ trainIndex,]
Test <- dataMatrix[-trainIndex,]

library (randomForest)
#set seed
set.seed(10)

# Performing Random Forest- create model based on train set - set number of tree as 1000

model = randomForest(as.factor(response)~.,data = Train,ntree = 1500,importance = T)
print(model)

#make predictions 

pred <- predict(model, Test, type = "prob", probability =TRUE)[,2]
preds <- prediction(pred, Test$response)

####---plots----####

#make ROC curve
library (ROCR)
#option1 
test_roc = roc(Test$response, pred, plot = TRUE, print.auc = TRUE)
#option 2
pROC_obj <- roc(Test$response, pred,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

# Make a data frame with predictor names and their importance
imp <- importance(model)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 10 predictors
imp.20 <- imp.sort[1:20, ]


# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying pre-nbUVB non-lesional samples\n into responders or non-responders")


####----Visit 2- Lesional----####

#lesional
# How many OTUs do we currently have? 
ntaxa(ps_v2_lesional_nbUVB.rarefied)

# Set minlib
minlib = 10000

prunescale = 0.0001


# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean <- taxa_sums(ps_v2_lesional_nbUVB.rarefied)/nsamples(ps_v2_lesional_nbUVB.rarefied)
sites.prune <- prune_taxa(tax.mean > prunescale*minlib, ps_v2_lesional_nbUVB.rarefied)

sites.prune

# Create the training and test datasets

dataMatrix <- data.frame(response = sample_data(sites.prune)$Response, otu_table(sites.prune))
View (dataMatrix)


###
library(caret)
set.seed(3456)

#Split data into training and test set- set to 70: 30
trainIndex <-  caret::createDataPartition(dataMatrix$response, p=0.70, list=FALSE)

# Step 2: Create the training  dataset

Train <- dataMatrix[ trainIndex,]
Test <- dataMatrix[-trainIndex,]

library (randomForest)
#set seed
set.seed(10)

# Performing Random Forest- create model based on train set - set number of tree as 1000

model = randomForest(as.factor(response)~.,data = Train,ntree = 1500,importance = T)
print(model)

#make predictions 

pred <- predict(model, Test, type = "prob", probability =TRUE)[,2]
preds <- prediction(pred, Test$response)

####---plots----####

#make ROC curve
library (ROCR)
#option1 
test_roc = roc(Test$response, pred, plot = TRUE, print.auc = TRUE)
#option 2
pROC_obj <- roc(Test$response, pred,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

# Make a data frame with predictor names and their importance
imp <- importance(model)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 10 predictors
imp.20 <- imp.sort[1:25, ]


# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying post-nbUVB lesional samples\n into responders or non-responders")


####---- Visit 2-Nonlesional-----
# How many OTUs do we currently have? 
ntaxa(ps_v2_non_nbUVB.rarefied)

# Set minlib
minlib = 10000

prunescale = 0.0001


# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean <- taxa_sums(ps_v2_non_nbUVB.rarefied)/nsamples(ps_v2_non_nbUVB.rarefied)
sites.prune <- prune_taxa(tax.mean > prunescale*minlib, ps_v2_non_nbUVB.rarefied)

sites.prune

# Create the training and test datasets

dataMatrix <- data.frame(response = sample_data(sites.prune)$Response, otu_table(sites.prune))
View (dataMatrix)


###
library(caret)
set.seed(3456)

#Split data into training and test set- set to 70: 30
trainIndex <-  caret::createDataPartition(dataMatrix$response, p=0.70, list=FALSE)

# Step 2: Create the training  dataset

Train <- dataMatrix[ trainIndex,]
Test <- dataMatrix[-trainIndex,]

library (randomForest)
#set seed
set.seed(3)

# Performing Random Forest- create model based on train set - set number of tree as 1000

model = randomForest(as.factor(response)~.,data = Train,ntree = 1500,importance = T)
print(model)

#make predictions 

pred <- predict(model, Test, type = "prob", probability =TRUE)[,2]
preds <- prediction(pred, Test$response)

####---plots----####

#make ROC curve
library (ROCR)
#option1 
test_roc = roc(Test$response, pred, plot = TRUE, print.auc = TRUE)
#option 2
pROC_obj <- roc(Test$response, pred,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)

# Make a data frame with predictor names and their importance
imp <- importance(model)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top 10 predictors
imp.20 <- imp.sort[1:30, ]


# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying post-nbUVB non-lesional samples\n into responders or non-responders")
