InitialData <- read.csv(file="gene-expression-invasive-vs-noninvasive-cancer (1).csv")
str(InitialData)
dim(InitialData) # To determine the shape of the dataset

dimnames(InitialData)[[2]][4947:4949]
table(InitialData[4949])

set.seed(2312122)

team.gene.subset <- rank(runif(1:4948))[1:2000] 
raw.gene.subset <- InitialData[team.gene.subset]

# Replace the class values of having 1 with 0's and having 2 with 1's
InitialData$Class[InitialData$Class == 1] <- 0
InitialData$Class[InitialData$Class == 2] <- 1

# Check null values
sum(is.na(raw.gene.subset))

# Find row number which contains missing values
rowSums(is.na(raw.gene.subset))

# Find out the names of the columns having missing values

mis_col <- names(which(colSums(is.na(raw.gene.subset))>0))
mis_col

# Plotting the box plot 
par(mfrow = c(1,3))
boxplot(raw.gene.subset[,'AB029013'])
boxplot(raw.gene.subset[,'Contig53370_RC'])
boxplot(raw.gene.subset[,'NM_004774'])

library(DMwR2)
library(tidyr)

# Filling null values using KNN Imputation method and Median values

new_ds <- knnImputation(raw.gene.subset)

#new_ds<- new_ds %>%
 # mutate_all(~replace_na(.,median(.,na.rm = TRUE)))
sum(is.na(new_ds))

# Part 1
# Dimensional Reduction Methods

# Supervised Learning

# To perform Two Sample T-test Method

library(dplyr)
t_test_results <- new_ds %>%
  
  #select(-gene_id) %>%

summarise_all(~t.test(. ~ InitialData$Class)$p.value)

t <- t_test_results[,t_test_results < 0.05]

sig.genes <- colnames(t)

reduced.data <- new_ds[sig.genes]
reduced.data$Class <- InitialData$Class    # Initializing class to the datasets
new_ds$Class <- InitialData$Class
ncol(reduced.data)
nrow(reduced.data)

# Unsupervised Learning 

# Variance Method

var_data <- apply(new_ds, 2, var)
var_data_col <- names(var_data[var_data > 0.25])
var_data_subset <- new_ds[, var_data_col]
ncol(var_data_subset)
nrow(var_data_subset)

# R-tsne Method

install.packages('Rtsne')
library(Rtsne)
tsne <- Rtsne(new_ds, perplexity=10)
tsne_data<-data.frame(tsne[['Y']])

# Part 2
# Groups of Patients

# PCA 
library(stats)
pca_mygenesubset <- princomp(var_data_subset, cor = TRUE) # Performing PCA using variance method dataset

summary(pca_mygenesubset)

pca_comp <- pca_mygenesubset$scores[,1:26]
ncol(pca_comp)

# plotting the bar graph
plot(pca_mygenesubset)


# K-Means Clustering

library(ggpubr)
library(factoextra)

#Scaling the patients data

reduced.data$Class <- NULL
scaled_data <- scale(reduced.data)


#cluster visualization 
cluster_result <- mygene_clust$cluster

#Elbow method
fviz_nbclust(reduced.data, # dataset 
             kmeans, # clustering algorithm 
             nstart = 25, 
             iter.max = 200, # the number of iterations allowed
             method = "wss")+geom_vline(xintercept = 4, linetype = 3 )

#plotting the cluster

mygene_clust <- kmeans(reduced.data, centers = 4, nstart = 20)
table(mygene_clust$cluster)

fviz_cluster(mygene_clust, data = reduced.data,
             palette = c("#EF731E","#3F6E9A","#AD2765", '#BC6213'), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

table(mygene_clust$cluster)


# Hierarchical clustering

dist_mat <- dist(scaled_data, method = "euclidean")
hclust_result <- hclust(dist_mat, method = "complete")

#visualizing hierarchial clustering
plot(hclust_result, main = "Hirarchial clustering method")

clusters_ <- cutree(hclust_result, k=2)
table(clusters_)

# Groups of Genes
transpose_t_data_subset <- t(reduced.data)
table(transpose_t_data_subset)

# PCA 
library(stats)
pca_t_mygenesubset <- princomp(transpose_t_data_subset, cor = TRUE)  # Performing PCA using variance method Transposed dataset

summary(pca_t_mygenesubset)

pca_t_comp <- pca_t_mygenesubset$scores[,1:37]
table(pca_t_comp)

#plotting the bar graph

plot(pca_t_mygenesubset)


# K-Means Clustering

library(ggpubr)
library(factoextra)

#Scaling the gene data

scaled_t_data <- scale(reduced.data)

mygene_t_clust <- kmeans(reduced.data, centers = 2, nstart = 50)

#visualizing the cluster 
cluster_t_result <- mygene_t_clust$cluster

#Elbow method
fviz_nbclust(reduced.data, # data  
             kmeans, # clustering algorithm 
             nstart = 25, 
             iter.max = 200, # the number of iterations allowed
             method = "wss")+ geom_vline(xintercept = 4, linetype = 3 )

#plotting the cluster

mygene_t_clust <- kmeans(reduced.data, centers = 4, nstart = 20)
table(mygene_t_clust$cluster)

fviz_cluster(mygene_t_clust, data = reduced.data,
             palette = c("#EF731E","#3F6E9A","#AD2765", '#BC6213'), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

table(mygene_t_clust$cluster)

# Hierarchial clustering

dist_t_mat <- dist(scaled_t_data, method = "euclidean")
hclust_t_result <- hclust(dist_t_mat, method = "complete")

#visualizing hierarchical clustering
plot(hclust_t_result, main = "Hierarchical clustering method")

clusters_t <- cutree(hclust_t_result, k=2)
table(clusters_t)

# Part 3

################################### Sampling using Two sample t-test data #############################################

# Logistic Regression

#train_data$Class<-new_ds$Class
library(caTools)

# Splitting the data to train and test
split <- sample.split(reduced.data, SplitRatio = 0.7, set.seed(2312122))
train_data <- subset(reduced.data, split == TRUE)
test_data <- subset(reduced.data, split == FALSE)

# Checking the null values
sum(is.na(train_data))
train_data<- na.omit(train_data)
test_data<- na.omit(test_data)
table(train_data$Class)
table(test_data$Class)

# Performing the Logistic regression model
log_model <- glm(Class ~ ., data = train_data, family = "binomial")
summary(log_model)

# Predicting the model
log_model.predict <- predict(log_model,newdata = test_data, type = "response")
summary(log_model.predict)
str(log_model.predict)

# Converting the predicted values to numeric values
log_numeric_data <- as.numeric(log_model.predict)
log_numeric_data
log_numeric_data <- na.omit(log_numeric_data)


predicted_class_log <- ifelse(log_numeric_data > 0.5,1,0)
predicted_class_log
test_data <- na.omit(test_data)
test_data

# Constructing confusion matrix
log_conf_matrix <- table(list(predicted_class_log, test_data$Class))
log_conf_matrix

# Calculating misclassification error
log_misclassification_error <- 1 - sum(diag(conf_matrix)) / sum(conf_matrix)
log_misclassification_error

# LDA

library(MASS)

# Performing the LDA model
lda_model_1 <- lda(Class ~ ., data = train_data)
summary(lda_model_1)

# Predicting the model
model_lda.predict <- predict(lda_model_1,newdata = test_data)$class
summary(model_lda.predict)

# Constructing confusion matrix
conf_matrix_lda <- table(model_lda.predict, test_data$Class)
conf_matrix_lda

# Calculating misclassification error
misclassification_error_lda <- 1 - sum(diag(conf_matrix_lda)) / sum(conf_matrix_lda)
misclassification_error_lda
 

# Random Forest
train_data <- na.omit(train_data)

library(randomForest)
e_rf = list() # Creating the list to store the errors

for (i in 1:15){
  
  # Performing the Random Forest model
  rf_model_1 <- randomForest(Class ~ ., data = train_data, ntree = i,set.seed(2312122)) 
  
  # Predicting the model
  predict.rf<- predict(rf_model_1, newdata = test_data, type = "response")
  
  str(predict.rf)
  
  # Converting the predicted values to numeric values
  numeric_rf <- as.numeric(predict.rf)
  predicted_class_rf <- ifelse(numeric_rf > 0.5, 1, 0)
  
  # Constructing confusion matrix
  rf_conf_matrix <- table(list(predicted_class_rf, test_data$Class))
  rf_conf_matrix
  
  # Calculating misclassification error
  rf_misclassification_error <- 1 - sum(diag(rf_conf_matrix)) / sum(rf_conf_matrix)  
  rf_misclassification_error
  e_rf <- append(e_rf,rf_misclassification_error)
}
e_rf


# SVM
library(caret)
library(e1071)

# Performing the SVM model
svm_mod <- svm(Class ~ ., data = train_data, kernel = "radial")

# Predicting the model
predict_svm <- predict(svm_mod, newdata = test_data)

# Summarize the predictions
summary(predict_svm)

# Converting predictions to numeric data
numeric_svm <- as.numeric(predict_svm)

# Assign class labels based on a threshold
svm_predicted_class <- ifelse(numeric_svm > 0.5, 1, 0)

conf_matrix_svm <- confusionMatrix(as.factor(svm_predicted_class), as.factor(test_data$Class))

# Constructing confusion matrix
conf_matrix_svm <- table(svm_predicted_class, test_data$Class)

conf_matrix_svm

# Calculating misclassification error
misclassification_svm <- 1 - sum(diag(conf_matrix_svm)) / sum(conf_matrix_svm)
misclassification_svm

# Naive Bayes

# Load the necessary library
library(e1071)

# performing the Naive Bayes model
model_nvb <- naiveBayes(Class ~ ., data = train_data)

# Make predictions on the model
predict_nvb <- predict(model_nvb, newdata = test_data)

# Summarize the predictions
summary(predict_nvb)

# Converting predicted values to numeric data
numeric_nvb <- as.numeric(predict_nvb)

# Assign class labels based on a threshold value
predicted_class_nvb <- ifelse(numeric_nvb > 0.5, 1, 0)

# generate confusion matrix
nvb_conf_matrix <- table(predicted_class_nvb, test_data$Class)

nvb_conf_matrix

# Calculate misclassification error
nvb_error <- 1 - sum(diag(nvb_conf_matrix)) / sum(nvb_conf_matrix)
nvb_error



# KNN

library(class)
test_data = na.omit(test_data)
e_knn <- list() # Creating list to store error values
for (i in 1:10){
  
  # performing the KNN model
  knn_m<- knn(train = train_data, test = test_data, cl=train_data$Class, k = i)
  summary(knn_model)
  
  # generate confusion matrix
  matrix_knn <- table(knn_m, test_data$Class)
  matrix_knn
  
  #Calculate misclassification error
  error_knn <- 1 - sum(diag(matrix_knn)) / sum(matrix_knn)  
  error_knn
  e_knn <- append(e_knn,error_knn)
}
e_knn
sum(is.na(train_data))
#na.omit(train_data)

# QDA

# Initializing Class Variable
tsne_data$Class<-new_ds$Class

# Split data to train and test
split_ts <- sample.split(tsne_data, SplitRatio = 0.7, set.seed(2312122))
train_data_ts <- subset(tsne_data, split == TRUE)
test_data_ts <- subset(tsne_data, split == FALSE)
sum(is.na(train_data_ts))
train_data_ts<- na.omit(train_data_ts)
test_data_ts<- na.omit(test_data_ts)


library(MASS)

# performing the Naive Bayes model
qda_model_1 <- qda(Class ~ ., data = train_data_ts)
summary(qda_model_1)

# Make predictions on the model
model_qda.predict <- predict(qda_model_1,newdata = test_data_ts)$class
summary(model_qda.predict)

# generate confusion matrix
conf_matrix_qda <- table(model_qda.predict, test_data_ts$Class)
conf_matrix_qda

#Calculate misclassification error
misclassification_error_qda <- 1 - sum(diag(conf_matrix_qda)) / sum(conf_matrix_qda)
misclassification_error_qda

########################################### Re-sampling K-Fold Technique #################################


# Logistic Regression

library(caret)

k <- 5
error <- list()

for (i in 1:(k-1)) {
  set.seed(2312122)
  indices <- sample(1:nrow(train_data), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Perform logistic regression model
  logistic_m_2 <- glm(Class ~ ., data = train_data[split != i, ], family = "binomial")
  
  # Prediction of the validation set
  predicted_log_2 <- predict(logistic_m_2, newdata = train_data[split == i, ], type = "response")
  
  # Convert predicted probabilities to class labels
  predicted_class_log_2 <- ifelse(predicted_log_2 > 0.5, 1, 0)
  
  # Create confusion matrix
  conf_matrix_log_2 <- confusionMatrix(as.factor(predicted_class_log_2), as.factor(train_data$Class[split == i]))
  
  # Calculate misclassification error
  log_Mis.error <- 1 - conf_matrix_log_2$overall["Accuracy"]
  
  error <- append(error, log_Mis.error)
}

error
e <- unlist(error)

# mean error
mean(e)

# KNN

library(caret)

# Initialize minimum error to 1
min_error <- 1
error <- list()
k_fold <- 5
# Sample indices
set.seed(2312122)
indices <- sample(1:nrow(train_data), replace = FALSE)
for (j in 1:50){
  for (i in 1:(k_fold-1)) {
    split <- cut(indices, breaks = k, labels = FALSE)
  
  # Perform KNN model
    knn_model_2 <- knn(train = train_data[split != i, -ncol(train_data)], 
                    test = train_data[split == i, -ncol(train_data)], 
                    cl = train_data$Class[split != i], 
                    k = j)  # Adjust k as needed
  
  # Prediction of the model
    predicted_knn_2 <- knn_model_2
  
  # Create confusion matrix
    conf_matrix_knn_2 <- confusionMatrix(as.factor(predicted_knn_2), 
                                             as.factor(train_data$Class[split == i]))
  
  # Calculate misclassification error
    knn_error2 <- 1 - conf_matrix_knn_2$overall["Accuracy"]
  
    error <- append(error, knn_error2)
  }

  error
  e <- unlist(error)
  mean_e<-mean(e)
  mean_e
  # Store minimum error in a variable
  if (mean_e < min_error){
    min_error<-mean_e
    v <- j
  }
}
min_error
v

# LDA

library(caret)

k <- 5
error <- list()

# Sampling the indices 
set.seed(2312122)
indices <- sample(1:nrow(train_data), replace = FALSE)

for (i in 1:(k-1)) {
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Perform LDA model
  lda_model_2 <- lda(Class ~ ., data = train_data[split != i, ])
  
  # Predicting the model
  predicted_lda_2 <- predict(lda_model_2, newdata = train_data[split == i, ])
  
  # Convert predicted class to labels
  class_lda_2 <- predicted_lda_2$class
  
  # Generating confusion matrix
  conf_matrix_lda_2 <- confusionMatrix(as.factor(class_lda_2), 
                                             as.factor(train_data$Class[split == i]))
  
  # Calculate misclassification error
  lda_error_2 <- 1 - conf_matrix_lda_2$overall["Accuracy"]
  
  error <- append(error, lda_error_2)
}

# Print errors
error
e <- unlist(error)

# Calculate mean error
mean(e)

#SVM

library(caret)
library(e1071)  

k <- 5
error <- list()

for (i in 1:(k-1)) {
  set.seed(2312122)
  indices <- sample(1:nrow(train_data), replace = FALSE)
  
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train SVM model
  model_2_svm <- svm(Class ~ ., data = train_data[split != i, ], kernel = "radial")
  
  # Predicting the validation set
  svm_predicted_2 <- predict(model_2_svm, newdata = train_data[split == i, ])
  
  predicted_class_svm_2 <- ifelse(svm_predicted_2 > 0.5, 1, 0)
  
  # generate confusion matrix
  conf_matrix_svm_2 <- confusionMatrix(as.factor(predicted_class_svm_2), 
                                             as.factor(train_data$Class[split == i]))
  
  # Calculate misclassification error
  svm_2_error <- 1 - conf_matrix_svm_2$overall["Accuracy"]
  
  # Append error to list
  error <- append(error, svm_2_error)
}

error
e <- unlist(error)

# Calculating mean error
mean(e)

# Random Forest

library(caret)

min_error <- 1
error <- list()
k_fold <- 5
# Sampling the indices 
set.seed(2312122)
indices <- sample(1:nrow(train_data), replace = FALSE)
for (j in 1:50) {
  for (i in 1:(k_fold-1)) {
    split <- cut(indices, breaks = k_fold, labels = FALSE)
    
    # Performing Random Forest model
    model_rf_2 <- randomForest(Class ~ ., data = train_data[split != i, ], ntree = j)  # Number of trees in the forest
    # Adjust parameters as needed
    
    # Predicting the data
    predicted_2_rf <- predict(model_rf_2, newdata = train_data[split == i,],type ="response")
    predicted_class_2_rf <- ifelse(predicted_2_rf > 0.5, 1, 0)
    
    # Creating confusion matrix
    rf_2_conf_matrix <- confusionMatrix(as.factor(predicted_class_2_rf), 
                                              as.factor(train_data$Class[split == i]))
    
    # Calculating the misclassification error
    rf_2_error <- 1 - rf_2_conf_matrix$overall["Accuracy"]
    
    # Appending error to list
    error <- append(error, rf_2_error)
  }
  
  
  e <- unlist(error)
  mean_e <- mean(e)
  # Storing the minimum value
  if (mean_e < min_error) {
    min_error <- mean_e
    var <- j 
  }
}

min_error
var


# Naive Bayes

library(caret)
library(e1071)

k <- 5
error <- list()

for (i in 1:(k-1)) {
  set.seed(2312122)
  indices <- sample(1:nrow(train_data), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Train Naive Bayes model
  model_2_nb <- naiveBayes(Class ~ ., data = train_data[split != i, ])
  
  # Predicting the validation set
  nb_2_predicted <- predict(model_2_nb, newdata = train_data[split == i, ])
  
  # Create the confusion matrix
  conf_matrix_nb_2 <- confusionMatrix(as.factor(nb_2_predicted), 
                                             as.factor(train_data$Class[split == i]))
  
  # Calculate the misclassification error
  nb_error_2 <- 1 - conf_matrix_nb_2$overall["Accuracy"]
  error <- append(error, nb_error_2)
}

error
e <- unlist(error)

# Calculate mean error
mean(e)


# QDA

library(caret)

k <- 5
error <- list()

# Sample indices once before the loop
set.seed(2312122)
indices <- sample(1:nrow(train_data_ts), replace = FALSE)

for (i in 1:(k-1)) {
  split_ts <- cut(indices, breaks = k, labels = FALSE)
  
  # Train QDA model
  qda_model_2 <- qda(Class ~ ., data = train_data_ts[split_ts != i, ])
  
  # Predict on the validation set
  predicted_qda_2 <- predict(qda_model_2, newdata = train_data_ts[split_ts == i, ])
  
  # Convert predicted class to labels
  class_qda_2 <- predicted_qda_2$class
  
  # Create confusion matrix
  conf_matrix_qda_2 <- confusionMatrix(as.factor(class_qda_2), 
                                       as.factor(train_data_ts$Class[split_ts == i]))
  
  # Calculate misclassification error
  qda_error_2 <- 1 - conf_matrix_qda_2$overall["Accuracy"]
  
  # Append error to list
  error <- append(error, qda_error_2)
}

# Print errors
error
e <- unlist(error)

# Calculate mean error
mean(e)

  ############################################# 2000 Gene Data #########################################################

# Logistic Regression

library(caTools)

# Split the data into train and test
split <- sample.split(new_ds, SplitRatio = 0.7, set.seed(2312122))
train_data_2000 <- subset(new_ds, split == TRUE)
test_data_2000 <- subset(new_ds, split == FALSE)
#train_data_2000$Class <- InitialData$Class

table(train_data_2000$Class)
table(test_data_2000$Class)

# Perform the Logistic Regression
model_log_3 <- glm(Class ~ ., data = train_data_2000, family = "binomial")
summary(model_log_3)

# Predicting the model
model_3_predict_log <- predict(model_log_3,newdata = test_data_2000, type = "response")
summary(model_3_predict_log)
str(model_3_predict_log)

# Converting the predicted values to numeric values
numeric_data_log<- as.numeric(model_3_predict_log)
numeric_data_log
numeric_data_log <- na.omit(numeric_data)
predicted_class_log_3 <- ifelse(numeric_data_log > 0.5,1,0)
predicted_class_log_3
test_data_2000 <- na.omit(test_data_2000)
test_data_2000

# Generating the confusion Matrix
conf_matrix_log_3 <- table(list(predicted_class_log_3, test_data_2000$Class))
conf_matrix_log_3

# Calculating the error
misclassification_error_log_3 <- 1 - sum(diag(conf_matrix_log_3)) / sum(conf_matrix_log_3)
misclassification_error_log_3

# LDA
#install.packages('MASS')
library(MASS)

# Perform the LDA model
model_lda_3 <- lda(Class ~ ., data = train_data_2000)
summary(model_lda_3)

# Predicting the model
model_3_predict_lda <- predict(model_lda_3,newdata = test_data_2000)$class
summary(model_3_predict_lda)
model_3_predict_lda
test_data_2000$Class

# Creating the confusion Matrix
conf_matrix_3_lda <- table(model_3_predict_lda, test_data_2000$Class)
conf_matrix_3_lda

# Calculating the misclassification error
error_3_lda <- 1 - sum(diag(conf_matrix_3_lda)) / sum(conf_matrix_3_lda)
error_3_lda

# Random Forest
train_data_2000 <- na.omit(train_data_2000)
library(randomForest)
e_rf = list() # Creating a list to store errors

for (i in 1:100){
  
  # Performing the random forest
  rf_3_model <- randomForest(Class ~ ., data = train_data_2000, ntree = i,set.seed(2312122))
  
  # Predicting the model
  model_3_rf<- predict(rf_3_model, newdata = test_data_2000, type = "response")  # Corrected: rf_model instead of model.predict
  
  str(model.predict.rf)
  
  # Converting the predicted values to numeric
  numeric_data_3_rf <- as.numeric(model_3_rf)
  predicted_class_3_rf <- ifelse(numeric_data_3_rf > 0.5, 1, 0)
  
  # Generating the confusion matrix
  conf_matrix_3_rf <- table(list(predicted_class_3_rf, test_data_2000$Class))
  conf_matrix_3_rf
  
  # Calculating the misclassification error
  error_3_rf <- 1 - sum(diag(conf_matrix_3_rf)) / sum(conf_matrix_3_rf)  
  error_3_rf
  e_rf <- append(e_rf,error_3_rf)
}
e_rf   

# SVM

library(e1071)

# Perform the SVM model
svm_model_3 <- svm(Class ~ ., data = train_data_2000, kernel = "radial")

# predicting the test data
predict_svm_3 <- predict(svm_model_3, newdata = test_data_2000)
summary(predict_svm_3)

# Convert predicted values to numeric data
numeric_data_svm_3 <- as.numeric(predict_svm_3)
predicted_class_svm_3 <- ifelse(numeric_data_svm_3 > 0.5, 1, 0)

# Creating the confusion matrix
conf_matrix_3_svm <- confusionMatrix(as.factor(predicted_class_svm_3), as.factor(test_data_2000$Class))
conf_matrix_3_svm <- table(predicted_class_svm_3, test_data_2000$Class)

conf_matrix_3_svm

# Calculating the misclassification error
error_svm_3 <- 1 - sum(diag(conf_matrix_3_svm)) / sum(conf_matrix_3_svm)
error_svm_3

# Naive Bayes

library(e1071)

# Performing the Naive Bayes model
model_3_nb <- naiveBayes(Class ~ ., data = train_data_2000)

# Predictions on validation set
model_pred_nb_3 <- predict(model_3_nb, newdata = test_data_2000)
summary(model_pred_nb_3)

# Convert predicted values to numeric data
numeric_nb_3 <- as.numeric(model_pred_nb_3)

# Assign class labels based on a threshold
class_nb_3 <- ifelse(numeric_nb_3 > 0.5, 1, 0)

# Creating the confusion matrix
conf_matrix_nb_3 <- table(class_nb_3, test_data_2000$Class)
conf_matrix_nb_3

# Calculating the misclassification error
error_nb_3 <- 1 - sum(diag(conf_matrix_nb_3)) / sum(conf_matrix_nb_3)
error_nb_3


# KNN

library(class)
train_data_2000 = na.omit(train_data_2000)
test_data_2000 = na.omit(test_data_2000)
e_knn_3 <- list() # Creating the list
for (j in (1:10)){
  # Performing the KNN Model
  knn_3_model <- knn(train = train_data_2000, test = test_data_2000, cl= train_data_2000$Class, k=j)
  
  # Creating the matrix
  conf_matrix_3_knn <- table(list(knn_3_model, test_data_2000$Class))
  conf_matrix_3_knn
  
  # Calculating the error
  error_3_knn <- 1 - sum(diag(conf_matrix_3_knn)) / sum(conf_matrix_3_knn)  
  error_3_knn
  e_knn_3 <- append(e_knn_3,error_3_knn)
}
e_knn_3


################################# Re-Sampling of 2000 Gene Data ################################################

# Logistic Regression

library(caret)

k <- 5
error <- list()

for (i in 1:(k-1)) {
  set.seed(2312122)
  indices <- sample(1:nrow(train_data_2000), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Performing the logistic regression model
  logistic_m_2 <- glm(Class ~ ., data = train_data_2000[split != i, ], family = "binomial")
  
  # Predicting the model
  predicted_log_2 <- predict(logistic_m_2, newdata = train_data_2000[split == i, ], type = "response")
  predicted_class_log_2 <- ifelse(predicted_log_2 > 0.5, 1, 0)
  
  # Creating the confusion matrix
  conf_matrix_log_2 <- confusionMatrix(as.factor(predicted_class_log_2), as.factor(train_data_2000$Class[split == i]))
  
  # Calculating the misclassification error
  log_Mis.error <- 1 - conf_matrix_log_2$overall["Accuracy"]
  error <- append(error, log_Mis.error)
}

error
e <- unlist(error)

# Calculate mean error
mean(e)

# KNN

library(caret)


min_error <- 1
error <- list()
k_fold <- 5
# Sampling the indices 
set.seed(2312122)
sum(is.na(train_data_2000))
indices <- sample(1:nrow(train_data_2000), replace = FALSE)
for (j in 1:50){
  for (i in 1:(k_fold-1)) {
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Perform KNN model
    knn_model_2 <- knn(train = train_data_2000[split != i, -ncol(train_data_2000)], 
                       test = train_data_2000[split == i, -ncol(train_data_2000)], 
                       cl = train_data_2000$Class[split != i], 
                       k = j)  # Adjust k value
    
    # Predict the model
    predicted_knn_2 <- knn_model_2
    
    # Creating the confusion matrix
    conf_matrix_knn_2 <- confusionMatrix(as.factor(predicted_knn_2), 
                                         as.factor(train_data_2000$Class[split == i]))
    
    # Calculating the misclassification error
    knn_error2 <- 1 - conf_matrix_knn_2$overall["Accuracy"]
    error <- append(error, knn_error2)
  }
  
  error
  e <- unlist(error)
  mean_e<-mean(e)
  mean_e  # Store minimum error in a variable
  if (mean_e < min_error){
    min_error<-mean_e
    v <- j
  }
}
min_error
v

# LDA

library(caret)

k <- 5
error <- list()

# Sample the indices
set.seed(2312122)
indices <- sample(1:nrow(train_data_2000), replace = FALSE)

for (i in 1:(k-1)) {
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Training the LDA model
  lda_model_2 <- lda(Class ~ ., data = train_data_2000[split != i, ])
  
  # Predicting the validation set
  predicted_lda_2 <- predict(lda_model_2, newdata = train_data_2000[split == i, ])
  
  class_lda_2 <- predicted_lda_2$class
  
  # Create the confusion matrix
  conf_matrix_lda_2 <- confusionMatrix(as.factor(class_lda_2), 
                                       as.factor(train_data_2000$Class[split == i]))
  
  # Calculate the misclassification error
  lda_error_2 <- 1 - conf_matrix_lda_2$overall["Accuracy"]
  error <- append(error, lda_error_2)
}

# Print errors
error
e <- unlist(error)

# Calculate mean error
mean(e)

#SVM
library(caret)
library(e1071)  

k <- 5
error <- list()

for (i in 1:(k-1)) {
  set.seed(2312122) # Sampling the indices
  indices <- sample(1:nrow(train_data_2000), replace = FALSE)
  
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Training the SVM model
  model_2_svm <- svm(Class ~ ., data = train_data_2000[split != i, ], kernel = "radial")
  
  # Predicting the validation set
  svm_predicted_2 <- predict(model_2_svm, newdata = train_data_2000[split == i, ])
  
  predicted_class_svm_2 <- ifelse(svm_predicted_2 > 0.5, 1, 0)
  
  # Creating the confusion matrix
  conf_matrix_svm_2 <- confusionMatrix(as.factor(predicted_class_svm_2), 
                                       as.factor(train_data_2000$Class[split == i]))
  
  # Calculate misclassification error
  svm_2_error <- 1 - conf_matrix_svm_2$overall["Accuracy"]
  error <- append(error, svm_2_error)
}

error
e <- unlist(error)

# Calculate mean error
mean(e)

# Random Forest

library(caret)

min_error <- 1
error <- list()
k_fold <- 5
# Sample the indices values
set.seed(2312122)
indices <- sample(1:nrow(train_data_2000), replace = FALSE)
for (j in 1:50) {
  for (i in 1:(k_fold-1)) {
    split <- cut(indices, breaks = k_fold, labels = FALSE)
    
    # Performing the Random Forest model
    model_rf_2 <- randomForest(Class ~ ., data = train_data_2000[split != i, ], ntree = j)
    
    # Predicting the validation set
    predicted_2_rf <- predict(model_rf_2, newdata = train_data_2000[split == i,],type ="response")
    predicted_class_2_rf <- ifelse(predicted_2_rf > 0.5, 1, 0)
    
    # Creating the confusion matrix
    rf_2_conf_matrix <- confusionMatrix(as.factor(predicted_class_2_rf), 
                                        as.factor(train_data_2000$Class[split == i]))
    
    # Calculating the misclassification error
    rf_2_error <- 1 - rf_2_conf_matrix$overall["Accuracy"]
  
    error <- append(error, rf_2_error)
  }
  
  
  e <- unlist(error)
  mean_e <- mean(e) # Store the minimum error in a variable
  if (mean_e < min_error) {
    min_error <- mean_e
    var <- j 
  }
}

min_error
var


# Naive Bayes

library(caret)
library(e1071)  

k <- 5
error <- list()

for (i in 1:(k-1)) {
  set.seed(2312122) # Sampling the indices
  indices <- sample(1:nrow(train_data_2000), replace = FALSE)
  split <- cut(indices, breaks = k, labels = FALSE)
  
  # Training the Naive Bayes model
  model_2_nb <- naiveBayes(Class ~ ., data = train_data_2000[split != i, ])
  
  # Predicting the trained model
  nb_2_predicted <- predict(model_2_nb, newdata = train_data_2000[split == i, ])
  
  # Creating the confusion matrix
  conf_matrix_nb_2 <- confusionMatrix(as.factor(nb_2_predicted), 
                                      as.factor(train_data_2000$Class[split == i]))
  
  # Calculating the misclassification error
  nb_error_2 <- 1 - conf_matrix_nb_2$overall["Accuracy"]
  
  error <- append(error, nb_error_2)
}

# Print errors
error
e <- unlist(error)
mean(e)

# Part 4

reduced.data$Class <- new_ds$Class
reduced.data$Cluster<- mygene_clust$cluster

# KNN Model

library(caret)

min_error <- 1
error <- list()
k_fold <- 5
# Sample the indices
set.seed(2312122)
indices <- sample(1:nrow(reduced.data), replace = FALSE)
for (j in 1:50){
  for (i in 1:(k_fold-1)) {
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Training the KNN model
    knn_model_2 <- knn(train = reduced.data[split != i, -ncol(reduced.data)], 
                       test = reduced.data[split == i, -ncol(reduced.data)], 
                       cl = reduced.data$Class[split != i], 
                       k = j)  # Adjust k the value
    
    # Predicting the validation set
    predicted_knn_2 <- knn_model_2
    
    # Creating the confusion matrix
    conf_matrix_knn_2 <- confusionMatrix(as.factor(predicted_knn_2), 
                                         as.factor(reduced.data$Class[split == i]))
    
    # Calculating the misclassification error
    knn_error2 <- 1 - conf_matrix_knn_2$overall["Accuracy"]
    
    error <- append(error, knn_error2)
  }
  
  error
  e <- unlist(error)
  mean_e<-mean(e)
  mean_e    # Storing the minimum error to a variable
  if (mean_e < min_error){
    min_error<-mean_e
    v <- j
  }
}
min_error
v

# Hierarchical Clustering

reduced.data$Cluster<- NULL
reduced.data$clusters_ <- clusters_

# KNN Model

library(caret)

min_error <- 1
error <- list()
k_fold <- 5
# Sample the indices 
set.seed(2312122)
indices <- sample(1:nrow(reduced.data), replace = FALSE)
for (j in 1:50){
  for (i in 1:(k_fold-1)) {
    split <- cut(indices, breaks = k, labels = FALSE)
    
    # Performing the KNN model
    knn_model_2 <- knn(train = reduced.data[split != i, -ncol(reduced.data)], 
                       test = reduced.data[split == i, -ncol(reduced.data)], 
                       cl = reduced.data$Class[split != i], 
                       k = j)  # Adjust k as needed
    
    # Predict the validation set
    predicted_knn_2 <- knn_model_2
    
    # Creating the confusion matrix
    conf_matrix_knn_2 <- confusionMatrix(as.factor(predicted_knn_2), 
                                         as.factor(reduced.data$Class[split == i]))
    
    # Calculating the misclassification error
    knn_error2 <- 1 - conf_matrix_knn_2$overall["Accuracy"]
    error <- append(error, knn_error2)
  }
  
  error
  e <- unlist(error)
  mean_e<-mean(e)
  mean_e    # Storing the minimum error to a variable
  if (mean_e < min_error){
    min_error<-mean_e
    v <- j
  }
}
min_error
v
