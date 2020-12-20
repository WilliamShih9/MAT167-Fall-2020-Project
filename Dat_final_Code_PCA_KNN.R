## ------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
dataset = read.csv("hour.csv") # dataset we are using
convert_hour_to_day <- function(hour){
  require(tidyverse)
  day = hour %>%
    group_by(dteday, season, yr, mnth, holiday, weekday, workingday) %>%
    summarize(weathersit = as.integer(round(mean(weathersit))),
              temp = mean(temp),
              atemp = mean(atemp),
              hum = mean(hum),
              windspeed = mean(windspeed),
              casual = sum(casual),
              registered = sum(registered),
              cnt = sum(cnt))
  #day = tibble::rowid_to_column(day, "instant")
  return(day)
}


## ---- echo = FALSE-------------------------------------------------------------------------------------------------------------

day = data.frame(convert_hour_to_day(dataset))
quant_variables = c('temp','atemp','hum','windspeed','casual',
                    'registered')
# function to fit PCA with default scaling is True
pca_func <- function(dat, scaling =TRUE){
 mod =  prcomp(dat,scale. = scaling)
 print(summary(mod))
 plot(cumsum(mod$sdev^2)/sum(mod$sdev^2),type='o',ylab='percentage of cumulative variance explained',xlab = 'number of components')
 return (mod)
}
# fitting to data with all quant variables
PCA = pca_func(day[,quant_variables])


## ---- echo = FALSE-------------------------------------------------------------------------------------------------------------
quant_data = PCA$x[,1:3]
qual_data = c("season", "yr", "mnth", "holiday" ,"weekday", "workingday", "weathersit")
full_data = cbind(quant_data,day[qual_data])
full_data$cnt = day[,"cnt"]#[['cnt']]
# coverting categorical data to factors
for ( j in qual_data){
  full_data[,j] = factor(full_data[,j])#[[j]])
}
set.seed(100)
# splitting data into training and test set where test size is 30%
test_index = sample(1:nrow(full_data),size = round(0.3*nrow(full_data)))
training_data = full_data[-test_index,]
test_data = full_data[test_index,]

# fitting a KNN on the training set to decide the number of neighbors
library(caret)

# analysis with 3 times repeated cross validation
control <- trainControl(method="repeatedcv",repeats = 3) 
# setting a grid of nearest neigbours from 1 to 20 on which the evaluation will be done
grid <- expand.grid(k=1:20)
# fitting knn
knn_training <- train(cnt ~ ., data = training_data, method = "knn", trControl = control,  tuneGrid = grid)

#plotting the knn model's RMSE
plot(knn_training)
print(knn_training)

# prediction using the best model
predictions_test <- predict(knn_training,newdata = test_data )
# plotting the predictions against true value
plot(x=test_data$cnt, y=predictions_test, xlab='true count',ylab='predicted count',main = 'predicted vs actual count',cex=0.6)
# getting the R-squared using predicted values
print(paste('The R-squared is ',cor(predictions_test,test_data$cnt)^2))

