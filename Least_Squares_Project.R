## ----setup, include=FALSE----------------------------------------------


#Sys.setenv(PATH = paste("C:/RBuildtools/4.0/bin", Sys.getenv("PATH"), sep=";"))
#Sys.setenv(BINPREF = "C:/RBuildtools/4.0/mingw_$(WIN)/bin/") 
knitr::opts_chunk$set(echo = FALSE)
library(pracma)
library(tidyverse)
library(scales)
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)

################ List of Functions

# There were two files from this dataset. There is the 'hour.csv' and 'day.csv'.
# This converts the 'hour.csv' to 'day.csv' file (summarizes the data from each hour to each day)
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
  day = tibble::rowid_to_column(day, "instant")
  return(day)
}

# Helper function for plot_by_hour for quantiles
# https://www.tidyverse.org/blog/2020/03/dplyr-1-0-0-summarise/
quibble <- function(x, q = c(0.25, 0.5, 0.75)) {
  tibble(x = quantile(x, q), q = q)
}
 
# Boxplot (lower end is 10th, lower end of box is 25th, middle line is 50th, upper end of box is 75th, upper end is 90th)
# Hour = This is the 'hour.csv' dataset
# Variable = This is the number of users. Either 'registered', 'casual', or 'cnt' (sum of registered or casual).
# Quantiles = Optional. Function needs exactly 5 values 
plot_by_hour <- function(hour, variable, quantiles = c(0.10, 0.25, 0.5, 0.75, 0.90)){
  require(tidyverse)
  require(scales)
  data = hour %>%
    group_by(hr) %>%
    summarise(x = quibble(get(variable), quantiles))
  data = tibble(hr = as.character(data$hr), x = data$x$x, q = data$x$q)
  quant_title = paste(ordinal(quantiles*100), collapse = ", ") 
  quant_title = paste0(quant_title, " percentile ", variable, " of each hour of day")
  plo = ggplot(data, aes(x = reorder(hr, sort(as.numeric(hr))), y = x)) + 
    theme_bw() +
    geom_boxplot(aes(fill = reorder(hr, sort(as.numeric(hr))))) +
    xlab("Hour of Day") +
    ylab(variable) +
    theme(legend.position = "none") +
    ggtitle(quant_title) 
  return(plo)
}

# Day = This is the 'day.csv' dataset
# Variables = These are the variables plotted. The variable plotted should be on the same scale
plot_by_day <- function(day, variables){
  require(tidyverse)
  require(reshape2)
  frame = data.frame(
    day = as.Date(day$dteday)
  )
  for (i in 1:length(variables)){
    frame = cbind(frame, day[variables[i]])
  }
  title = paste(variables, collapse = ", ")
  title = paste0("Date and Value of ", title)
  frame_gg <- melt(frame, id="day")
  plo = ggplot(frame_gg, aes(x = day, y=value, color = variable)) +
    geom_line()+
    xlab("Date") +
    scale_x_date(date_breaks = "3 month", date_labels = "%m-%Y") +
    ylab("Value") +
    theme_bw() +
    ggtitle(title)
  return(plo)
}


# Finds the design matrix. No conversion to dummy variables allowed
# Hour = This is the 'hour.csv' dataset
# x = This is the independent variable(s)
# degrees = The degrees of the polynomials for the independent variable(s): Must be 1 or higher
design_matrix <- function(hour, x, degrees = rep(1, length(x))){
  mat = matrix(, nrow = nrow(hour), ncol = sum(degrees)+1)
  name = c("(Intercept)")
  mat[,1] = 1
  k = 2
  for (i in 1:length(degrees)){
    data = hour[[x[i]]]
    for (j in 1:degrees[i]){
      mat[,k] = data^j
      if (j == 1){
         name = append(name, x[i])
      }
      else{
         name = append(name, paste0(x[i], "^", j))
      }
      k = k + 1
    }
  }
  colnames(mat) = name
  return(mat)
}

# C++ version of the design_matrix function
Rcpp::cppFunction('
NumericMatrix design_matrix_Cpp(List hour, std::vector<std::string> x, IntegerVector degrees){
  int total = 1;
  for (int i = 0; i < degrees.size(); i++){
    total = total + degrees[i];
  }
  CharacterVector name(total);
  name[0] = "(Intercept)";
  NumericVector v = hour[0];
  NumericMatrix mat(v.length(), total);
  std::fill(mat.begin(), mat.begin()+v.length(), 1);
  int k = 1;
  for (int i = 0; i < degrees.size(); i++){
    NumericVector v = hour[x[i]];
    for (int j = 1; j <= degrees[i]; j++){
      mat(_,k) = pow(v, j);
      if (j == 1){
        name[k] = x[i];
      }
      else{
        name[k] = x[i] + "^" + std::to_string(j);
      }
      k++;
    }
    
  }
  colnames(mat) = name;
  return(mat);
}
')

#Finds x for Ax = b using normal equations
normal_equations <- function(A, b){
  x = solve(t(A) %*% A) %*% t(A) %*% b
  return(x)
}

# C++ version for normal_equations
Rcpp::cppFunction('
  arma::vec normal_equations_Cpp(arma::mat A, arma::vec b){
    arma::vec coeffs = inv(trans(A) * A) * trans(A) * b;
    return(coeffs);
  }
', depends='RcppArmadillo')

# C++ version for qr.solve
Rcpp::cppFunction('
  arma::vec qr_solve_Cpp(arma::mat A, arma::vec b){
    arma::mat Q;
    arma::mat R;
    arma::qr_econ(Q, R, A);
    arma::vec x = inv(R)*trans(Q)*b;
    return(x);
  }
', depends = 'RcppArmadillo')

# Finds x for Ax=b using SVD
svd_solve <- function(A, b){
  x = svd(A)
  x = as.numeric(x$v %*% (t(x$u) %*% b / x$d))
  names(x) = colnames(A)
  return(x)
}

# C++ version of svd_solve
Rcpp::cppFunction('
  arma::vec svd_solve_Cpp(arma::mat A, arma::vec b){
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd_econ(U, s, V, A);
    arma::vec x = V * (trans(U)*b/s);
    return(x);
  }
', depends = 'RcppArmadillo')



# Convert condition numbers to LaTeX output
# Digits can be from 0 to 7 inclusive
scientific_to_latex <- function(x, digits = 3){
  x = scientific(as.numeric(x), digits = digits)
  x = paste0("$", substr(x, 1, 1+digits), " \\times 10^{", as.numeric(substr(x, stringr::str_length(x)-2, stringr::str_length(x))), "}$")
  return(x)
}

# Used for curve fitting
# Coefficients = coeffieints for the polynomial term: The first value is power of one, second value is quadratic, third value is cubic, etc..
# Constant = constant term added to each polynomial term
# Range = x values
polynomial <- function(coefficients, constant = 0, range = seq(0, 23, 0.1)){
  A = vector(mode = "numeric", length = length(range))
  for (i in c(1:length(range))){
    A[i] = constant
    for (j in c(1:length(coefficients))){
      A[i] = A[i] + coefficients[j]*range[i]^j
    }
  }
  return(data.frame(range, A))
}

######################## End List of Functions

# 1. Introduction


## ---- echo = FALSE, message = FALSE,  fig.pos = "!h", fig.cap = "Boxplot of registered users for Jan 1, 2011 to December 31, 2011 for each hour of day. The low whisker represents 10th percentile, the box represents the interquartile range, the top whisker represents the 90th percentile. For example, with 731 days in the datset, the low whisker represents the 73rd lowest value.",  fig.width = 8, fig.height = 4----
#dataset = read.csv("hour.csv")
temp <- tempfile()
download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00275/Bike-Sharing-Dataset.zip", temp)
dataset = read.csv(unz(temp, "hour.csv"))

#Figure 1 
plot_by_hour(dataset, "registered")
Figure1 = plot_by_hour(dataset, "registered")


## ---- echo = FALSE, message = FALSE,  fig.pos = "!h", fig.cap = "Boxplot of casual users for Jan 1, 2011 to December 31, 2011 for each hour of day. The low whisker represents 10th percentile, the box represents the interquartile range, the top whisker represents the 90th percentile. For example, with 731 days in the datset, the low whisker represents the 73rd lowest value.",  fig.width = 8, fig.height = 4----

# Figure 2
plot_by_hour(dataset, "casual")
Figure2 = plot_by_hour(dataset, "casual")

## ---- echo = FALSE, message = FALSE,  fig.pos = "!h",  fig.cap = "Number of casual, registered, and sum of casual and registered for all 731 days in the dataset", fig.width = 8, fig.height = 4----
day = convert_hour_to_day(dataset)

# Figure 3
plot_by_day(day, c("casual","registered","cnt"))
Figure3 = plot_by_hour(dataset, "casual")


## ---- echo = FALSE, message = FALSE, fig.pos = "!h", fig.cap = 'Normalized temperature and "feels like" temperature for all 731 days in the datset',  fig.width = 8, fig.height = 4----

# Figure 4
plot_by_day(day, c("temp","atemp"))
Figure4 = plot_by_hour(dataset, "casual")


## ---- message = FALSE, echo = FALSE------------------------------------

#3. Design Matrix 

# This is calculating R-squared of regression models 

Rsq_cnt = sapply(1:15, function(x) {
  summary(lm(dataset$cnt ~ design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,x,1)))))$r.squared
})

Rsq_casual = sapply(1:15, function(x) {
  summary(lm(dataset$casual ~ design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,x,1)))))$r.squared
})

Rsq_registered = sapply(1:15, function(x) {
  summary(lm(dataset$registered ~ design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,x,1)))))$r.squared
})

condition = sapply(1:15, function(x) {
  cond(design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,x,1))))
})

frame = cbind(Rsq_cnt, Rsq_casual, Rsq_registered)

hs = c("$x_{4}$",
  "$x_{4}+x_{4}^{2}$",
  "$x_{4}+x_{4}^{2}+x_{4}^{3}$",
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{4}$",
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{5}$",
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{6}$",  
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{7}$",  
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{8}$",  
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{9}$",  
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{10}$",  
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{11}$",  
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{12}$", 
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{13}$",  
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{14}$",
  "$x_{4}+x_{4}^{2}+\\dots+x_{4}^{15}$"
)
rownames(frame) = hs

frame[,1:3] = round(frame[,1:3], digits = 3)


colnames(frame) = c("$R^{2}_{cnt}$", "$R^{2}_{casual}$", "$R^{2}_{registered}$")

# Table 1
knitr::kable(frame, digits = 3, caption = "R-squared values for linear regression of each of value of $n$, maximum power of polynomial for `hour`, including $x_{1}, x_{2}$, and $x_{9}$", escape = FALSE)


## ---- message = FALSE, echo = FALSE------------------------------------

# Comparing speed of custom design_matrix function versus built in mode.matrix.lm function
Compare = microbenchmark("design_matrix()" = {S = design_matrix(dataset, c("yr","season","hr","temp"), c(1,1,7,1))},
                       "design_matrix_Cpp()" = {S = design_matrix_Cpp(dataset, c("yr","season","hr","temp"),
                                                                  as.integer(c(1,1,7,1)))},
                       "model.matrix.lm()" = 
                         {S = model.matrix(~yr+mnth+hr+I(hr^2)+I(hr^3)+I(hr^4)+I(hr^5)+I(hr^6)+I(hr^7)+temp,
                                                          data = dataset)},
                       times = 100)
# Table 2
knitr::kable(summary(Compare), digits = 3, caption = "Comparison of speed (in milliseconds) of custom design_matrix() functions with built-in R model.matrix.lm() function", escape = FALSE)



## ---- message = FALSE, warning = FALSE, echo = FALSE-------------------

#4. Normal Equations
A = design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,3,1)))
b = dataset$casual

# Relative Error in R
result = sapply(1:15, function(x)  {
  A = design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,x,1)))
  x = tryCatch(expr = 
                 {x = solve(t(A) %*% A) %*% t(A) %*% b
                  x = max((svd_solve(A, b) - x)/svd_solve(A,b))},
               error = function(e){
                 str = paste0('Error, Recripocal $\\kappa(A)$: ', substr(e[[1]], 67, 999))
                 return(str)
               })
  return(x)
})

# Relative Error in C++
result2 = sapply(1:15, function(x)  {
  A = design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,x,1)))
  x = tryCatch(expr = 
                 {x = normal_equations_Cpp(A, b)
                  x = max((svd_solve(A, b) - x)/svd_solve(A,b))},
               error = function(e){
                 str = paste0('Error, Recripocal $\\kappa(A)$: ', substr(e[[1]], 67, 999))
                 return(str)
               })
  return(x)
})


# Convert Results to Latex Format
column1 = scientific_to_latex(condition, 3)
column2 = scientific_to_latex(condition^2, 3)
column3 = scientific_to_latex(result[1:5], 3)
column3b = scientific_to_latex(substr(result[6:15], 32, 999), 3)
column3c = paste0(substr(result[6:15], 1, 31), column3b)
column3d = c(column3, column3c)
column4 = scientific_to_latex(result2, 3)
frame = cbind(column2, column3d, column4)

colnames(frame) = c("$\\kappa(A)^{2}$", "Relative error/error message R", "Relative error C++")

# Table 3
rownames(frame) = hs

knitr::kable(frame, escape = FALSE, caption = "Condition number of each value of $n$ (maximum power of polynomial for `hour`) and relative error of normal equations versus SVD")


## ---- message = FALSE, echo = FALSE------------------------------------
Compare2 = microbenchmark("normal_equations(A, b)" = {S = normal_equations(A, b)},
                       "normal_equations_Cpp(A, b)" = {S = normal_equations_Cpp(A, b)},
                        times = 100)
# Table 4
knitr::kable(summary(Compare2), digits = 3, caption = "Comparison of speed (in milliseconds) of solving $Ax = b$ where $A \\in \\mathbb{R}^{17379 \\times 7}$ using normal equations implemented in R vs Rcpp (C++)")


## ---- message = FALSE, echo = FALSE, message = FALSE-------------------
#5. QR decomposition

Compare3 = microbenchmark("qr.solve(A, b)" = {S = qr.solve(A, b)},
                       "qr_solve_Cpp(A, b)" = {S = qr_solve_Cpp(A, b)},
                        times = 100)
# Table 5
knitr::kable(summary(Compare3), digits = 3, caption = "Comparison of speed (in milliseconds) of solving $Ax = b$ where $A \\in \\mathbb{R}^{17379 \\times 7}$ using QR decomposition implemented in R vs Rcpp (C++)")



## ---- message = FALSE, echo = FALSE------------------------------------
#6. Singular Value Decomposition
Compare4 = microbenchmark("svd_solve(A, b)" = {S = svd_solve(A, b)},
                       "svd_solve_Cpp(A, b)" = {S = svd_solve_Cpp(A, b)},
                        times = 100)
# Table 6
knitr::kable(summary(Compare4), digits = 3, caption = "Comparison of speed (in milliseconds) of solving $Ax = b$ where $A \\in \\mathbb{R}^{17379 \\times 7}$ using SVD implemented in R vs Rcpp (C++)")



## ---- message = FALSE, echo = FALSE------------------------------------
A = design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,7,1)))

Result2 = microbenchmark("normal_equations_Cpp(A, b)" = {S = normal_equations_Cpp(A, dataset$cnt)},
                        "qr.solve(A, b)" = {S = qr.solve(A, dataset$cnt)},
                       "qr_solve_Cpp(A, b)" = {S = qr_solve_Cpp(A, b)},
                       "svd_solve(A, b)" = {S = svd_solve(A, dataset$cnt)},
                       "svd_solve_Cpp(A,b)" = {S = svd_solve_Cpp(A, b)},
                       times = 100)
# Table 7
knitr::kable(summary(Result2), digits = 3, caption = "Comparison of speed (in milliseconds) of solving $Ax = b$ where $A \\in \\mathbb{R}^{17379 \\times 11}$ using QR decomposition or SVD implemented in R vs Rcpp (C++)")



## ---- echo = FALSE, message = FALSE, fig.cap = "The 25th, 50th, and 75th percentile of casual users (731 days) are graphed for each hour of the day. The mean values of season, yr, and temp are entered into the regression solution to give a constant value. Thus, the regression is now only in terms of hour. Then, the mean curve fit of the regression can be compared to the 25th, 50th, and 75th percentiles of casual users.", fig.width = 8, fig.height = 4----
# 7. Final Regression Model
y1 = svd_solve(design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,3,1))), dataset$casual)
y3 = svd_solve(design_matrix_Cpp(dataset, c("yr","season","hr","temp"), as.integer(c(1,1,7,1))), dataset$cnt)

constant_y1 = y1[1] + mean(dataset$yr)*y1[2] + mean(dataset$season)*y1[3] + mean(dataset$temp)*y1[7]

constant_y3 = y3[1] + mean(dataset$yr)*y3[2] + mean(dataset$season)*y3[3] +
mean(dataset$temp)*y3[11]

coeff_y1 = y1[4:6]

coeff_y3 = y3[4:10]

# Fitting curve for casual and cnt
predicted_casual = polynomial(coeff_y1, constant_y1)
predicted_cnt = polynomial(coeff_y3, constant_y3)
predicted_casual = cbind(predicted_casual, "Curve Fit")
predicted_cnt = cbind(predicted_cnt, "Curve Fit")
colnames(predicted_casual) = c("hr","x", "q")
colnames(predicted_cnt) = c("hr", "x", "q")
real_casual = dataset %>%
  group_by(hr) %>%
  summarise(x = quibble(casual, c(0.25,0.5,0.75)))
real_cnt = dataset %>%
  group_by(hr) %>%
  summarise(x = quibble(cnt, c(0.25,0.5,0.75))) 

real_casual$q = paste0(sprintf("%.f", real_casual$x$q*100), "th")
real_cnt$q = paste0(sprintf("%.f", real_cnt$x$q*100), "th")
real_casual$x = real_casual$x$x
real_cnt$x = real_cnt$x$x

casual = rbind(predicted_casual, real_casual)
cnt = rbind(predicted_cnt, real_cnt)

#Figure 5
Figure5 = ggplot(casual, aes(x = hr, y = x, group = q, color = q)) +
  geom_line() +
  xlab("Hour") +
  ylab("Number of casual users") +
  ggtitle("25th, 50th, 75th percentile of casual users compared to mean curve fit (as function of hour)") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2)) +
  theme_bw()
Figure5


## ---- echo = FALSE, message = FALSE, fig.cap = "Same description as for Figure 5, except the variable is now cnt (sum of registered and casual users) rather than casual.", fig.width = 8, fig.height = 4----
#Figure 6
Figure6 = ggplot(cnt, aes(x = hr, y = x, group = q, color = q)) +
  geom_line() +
  xlab("Hour") +
  ylab("Number of casual+registered users") +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2)) +
  ggtitle("25th, 50th, 75th percentile of cnt users compared to mean curve fit (as function of hour)") +
  theme_bw()
Figure6

## ------------------------------------------------------------------------------------------------------------------------------
# 8. Principal Components Analysis
library(tidyverse)
#dataset = read.csv("hour.csv") # dataset we are using
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
# 9. k-Nearest Neighbors
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
# Figure 8
plot(x=test_data$cnt, y=predictions_test, xlab='true count',ylab='predicted count',main = 'predicted vs actual count',cex=0.6)
# getting the R-squared using predicted values
print(paste('The R-squared is ',cor(predictions_test,test_data$cnt)^2))


