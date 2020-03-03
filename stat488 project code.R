###########
#MERGING DATASETS
###########
#commented out because this was already done.

#load training data
# train <- read.csv('train.csv')
# weather.train <- read.csv('weather_train.csv')
# building <- read.csv('building_metadata.csv')

#converting timestamps from factor to actual time stamps
# library(lubridate)
# weather.train$timestamp <-ymd_hms(weather.train$timestamp)
# train$timestamp <- ymd_hms(train$timestamp)

#merging training data
# merge_train <- merge(train, building, by.x = 'building_id', by.y = 'building_id', all.x = TRUE)
# rm(train)
# train <- merge(merge_train, weather.train, by.x = c('site_id', 'timestamp'), by.y = c('site_id', 'timestamp'), all.x = TRUE)
# rm(weather.train)
# rm(merge_train)
# save.image("ashrae_train.RData")
##############################################################################################################
#load packages
library(lubridate)
library(ggplot2)
library(imputeTS)
library(leaps)
library(gam)
library(gbm)

#loading the data
#load('ashrae_train.RData') #would have been used if I had more memory
train <- get(load('merg_train_small.RData'))
rm(merg_train_small)
train$timestamp <- ymd_hms(train$timestamp) #converts from factor to POSIXct to save memory

#Initial EDA
summary(train)

#removed because of too many missing values
train <- subset(train, select = -c(floor_count, year_built))

#imputation
ggplot(aes(x=timestamp, y=dew_temperature), data=train) + geom_point() + geom_smooth()
ggplot(aes(x=timestamp, y=air_temperature), data=train) + geom_point() + geom_smooth()
ggplot(aes(x=timestamp, y=sea_level_pressure), data=train) + geom_point() + geom_smooth()
ggplot(aes(x=timestamp, y=wind_direction), data=train) + geom_point() + geom_smooth()
ggplot(aes(x=timestamp, y=wind_speed), data=train) + geom_point() + geom_smooth()
ggplot(aes(x=timestamp, y=precip_depth_1_hr), data=train) + geom_point() + geom_smooth()
ggplot(aes(x=timestamp, y=cloud_coverage), data=train) + geom_point() + geom_smooth()

train <- na_mean(train, option = 'median')

#encoding primary use
train$primary_use <- as.numeric(train$primary_use)
set.seed(488)
##############################################################################################################

#METER ZERO: ELECTRICITY
meter0.train <- subset(train, meter == 0)

#more EDA
ggplot(aes(x=dew_temperature, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=air_temperature, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=cloud_coverage, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=primary_use, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=square_feet, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=wind_speed, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=wind_direction, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=sea_level_pressure, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=precip_depth_1_hr, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=building_id, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()
ggplot(aes(x=site_id, y=meter_reading), data=meter0.train) + geom_point() + geom_smooth()

#Backward Selection
models <- regsubsets(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + precip_depth_1_hr + sea_level_pressure + wind_direction + wind_speed + cloud_coverage + building_id + site_id, data = meter0.train, method = 'backward', nvmax = 11)
best.features <- summary(models)
data.frame(num_features = which.min(best.features$bic), BIC = min(best.features$bic))
best.features

#split into training and validation sets
meter0.training <- meter0.train[meter0.train$timestamp >= '2016-01-31 18:00:00',]
meter0.validation <- meter0.train[meter0.train$timestamp < '2016-01-31 18:00:00',]

###########
#GAMS
###########
#natural cubic splines
rmlse <- c(rep(NA, 10))
for (d in 1:10) {
  meter0.gam <- gam(meter_reading ~ primary_use + ns(square_feet, df=d) + ns(air_temperature, df=d) + ns(cloud_coverage, df=d) + ns(dew_temperature, df=d) + ns(wind_direction, df=d) + building_id + site_id, data = meter0.training)
  
  meter0.gam.pred <- predict(meter0.gam, newdata = meter0.validation)
  meter0.gam.pred[meter0.gam.pred < 0] = 0
  rmlse[d] <- sqrt((1/nrow(meter0.train)) * sum((log(meter0.gam.pred+1) - log(meter0.validation$meter_reading+1))^2))
}

min(rmlse)
which.min(rmlse)
plot(1:10, rmlse, type = 'b')

#smoothing splines
rmlse <- c(rep(NA, 10))
for (d in 1:10) {
  meter0.gam <- gam(meter_reading ~ primary_use + s(square_feet, df=d) + s(air_temperature, df=d) + s(cloud_coverage, df=d) + s(dew_temperature, df=d) + s(wind_direction, df=d) + building_id + site_id, data = meter0.training)
  
  meter0.gam.pred <- predict(meter0.gam, newdata = meter0.validation)
  meter0.gam.pred[meter0.gam.pred < 0] = 0
  rmlse[d] <- sqrt((1/nrow(meter0.train)) * sum((log(meter0.gam.pred+1) - log(meter0.validation$meter_reading+1))^2))
}

min(rmlse)
which.min(rmlse)
plot(1:10, rmlse, type = 'b')

###########
#GBM
###########
set.seed(488)
#all features
rmlse <- c(rep(NA, 3))
lambdas <- c(0.001, 0.01, 0.1)
i <- 1
for (lambda in lambdas){
  meter0.gbm <- gbm(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + precip_depth_1_hr + sea_level_pressure + wind_direction + wind_speed + cloud_coverage + building_id + site_id, data = meter0.training, shrinkage = lambda, cv.folds = 10)
  meter0.gbm.pred <- predict(meter0.gbm, newdata = meter0.validation)
  meter0.gbm.pred[meter0.gbm.pred < 0] = 0
  rmlse[i] <- sqrt((1/nrow(meter0.train)) * sum((log(meter0.gbm.pred+1) - log(meter0.validation$meter_reading+1))^2))
  i <- i+1
}

min(rmlse)
lambdas[which.min(rmlse)]
plot(lambdas, rmlse, type = 'b')

#subset of features based on BIC
rmlse <- c(rep(NA, 3))
lambdas <- c(0.001, 0.01, 0.1)
i <- 1
for (lambda in lambdas){
  meter0.gbm <- gbm(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + wind_direction + cloud_coverage + building_id + site_id, data = meter0.training, shrinkage = lambda, cv.folds = 10)
  meter0.gbm.pred <- predict(meter0.gbm, newdata = meter0.validation)
  meter0.gbm.pred[meter0.gbm.pred < 0] = 0
  rmlse[i] <- sqrt((1/nrow(meter0.train)) * sum((log(meter0.gbm.pred+1) - log(meter0.validation$meter_reading+1))^2))
  i <- i+1
}

min(rmlse)
lambdas[which.min(rmlse)]
plot(lambdas, rmlse, type = 'b')
##############################################################################################################

#METER ONE: CHILLED WATER
meter1.train <- subset(train, meter == 1)

#more EDA
ggplot(aes(x=dew_temperature, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=air_temperature, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=cloud_coverage, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=primary_use, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=square_feet, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=wind_speed, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=wind_direction, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=sea_level_pressure, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=precip_depth_1_hr, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=building_id, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()
ggplot(aes(x=site_id, y=meter_reading), data=meter1.train) + geom_point() + geom_smooth()

#Backward Selection
models <- regsubsets(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + precip_depth_1_hr + sea_level_pressure + wind_direction + wind_speed + cloud_coverage + building_id + site_id, data = meter1.train, method = 'backward', nvmax = 11)
best.features <- summary(models)
data.frame(num_features = which.min(best.features$bic), BIC = min(best.features$bic))
best.features

#split into training and validation sets
meter1.training <- meter1.train[meter1.train$timestamp >= '2016-01-31 18:00:00',]
meter1.validation <- meter1.train[meter1.train$timestamp < '2016-01-31 18:00:00',]

###########
#GAMS
###########
rmlse <- c(rep(NA, 10))
for (d in 1:10) {
  meter1.gam <- gam(meter_reading ~ primary_use + ns(square_feet, df=d) + ns(dew_temperature, df=d) + ns(sea_level_pressure, df=d) + site_id, data = meter1.training)
  
  meter1.gam.pred <- predict(meter1.gam, newdata = meter1.validation)
  meter1.gam.pred[meter1.gam.pred < 0] = 0
  rmlse[d] <- sqrt((1/nrow(meter1.train)) * sum((log(meter1.gam.pred+1) - log(meter1.validation$meter_reading+1))^2))
}

min(rmlse)
which.min(rmlse)
plot(1:10, rmlse, type = 'b')

#smoothing splines
rmlse <- c(rep(NA, 10))
for (d in 1:10) {
  meter1.gam <- gam(meter_reading ~ primary_use + s(square_feet, df=d) + s(dew_temperature, df=d) + s(sea_level_pressure, df=d) + site_id, data = meter1.training)
  
  meter1.gam.pred <- predict(meter1.gam, newdata = meter1.validation)
  meter1.gam.pred[meter1.gam.pred < 0] = 0
  rmlse[d] <- sqrt((1/nrow(meter1.train)) * sum((log(meter1.gam.pred+1) - log(meter1.validation$meter_reading+1))^2))
}

min(rmlse)
which.min(rmlse)
plot(1:10, rmlse, type = 'b')

###########
#GBM
###########
set.seed(488)
#all features
rmlse <- c(rep(NA, 3))
lambdas <- c(0.001, 0.01, 0.1)
i <- 1
for (lambda in lambdas){
  meter1.gbm <- gbm(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + precip_depth_1_hr + sea_level_pressure + wind_direction + wind_speed + cloud_coverage + building_id + site_id, data = meter1.training, shrinkage = lambda, cv.folds = 10)
  meter1.gbm.pred <- predict(meter1.gbm, newdata = meter1.validation)
  meter1.gbm.pred[meter1.gbm.pred < 0] = 0
  rmlse[i] <- sqrt((1/nrow(meter1.train)) * sum((log(meter1.gbm.pred+1) - log(meter1.validation$meter_reading+1))^2))
  i <- i+1
}

min(rmlse)
lambdas[which.min(rmlse)]
plot(lambdas, rmlse, type = 'b')

#subset of features based on BIC
rmlse <- c(rep(NA, 3))
lambdas <- c(0.001, 0.01, 0.1)
i <- 1
for (lambda in lambdas){
  meter1.gbm <- gbm(meter_reading ~ primary_use + square_feet + dew_temperature + sea_level_pressure + site_id, data = meter1.training, shrinkage = lambda, cv.folds = 10)
  meter1.gbm.pred <- predict(meter1.gbm, newdata = meter1.validation)
  meter1.gbm.pred[meter1.gbm.pred < 0] = 0
  rmlse[i] <- sqrt((1/nrow(meter1.train)) * sum((log(meter1.gbm.pred+1) - log(meter1.validation$meter_reading+1))^2))
  i <- i+1
}

min(rmlse)
lambdas[which.min(rmlse)]
plot(lambdas, rmlse, type = 'b')
##############################################################################################################

#METER TWO: STEAM
meter2.train <- subset(train, meter == 2)

#more EDA
ggplot(aes(x=dew_temperature, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=air_temperature, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=cloud_coverage, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=primary_use, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=square_feet, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=wind_speed, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=wind_direction, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=sea_level_pressure, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=precip_depth_1_hr, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=building_id, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()
ggplot(aes(x=site_id, y=meter_reading), data=meter2.train) + geom_point() + geom_smooth()

#Backward Selection
models <- regsubsets(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + precip_depth_1_hr + sea_level_pressure + wind_direction + wind_speed + cloud_coverage + building_id + site_id, data = meter2.train, method = 'backward', nvmax = 11)
best.features <- summary(models)
data.frame(num_features = which.min(best.features$bic), BIC = min(best.features$bic))
best.features

#split into training and validation sets
meter2.training <- meter2.train[meter2.train$timestamp >= '2016-01-31 18:00:00',]
meter2.validation <- meter2.train[meter2.train$timestamp < '2016-01-31 18:00:00',]

###########
#GAMS
###########
rmlse <- c(rep(NA, 10))
for (d in 1:10) {
  meter2.gam <- gam(meter_reading ~ primary_use + ns(square_feet, df=d) + ns(sea_level_pressure, df=d) + building_id + site_id, data = meter2.training)
  
  meter2.gam.pred <- predict(meter2.gam, newdata = meter2.validation)
  meter2.gam.pred[meter2.gam.pred < 0] = 0
  rmlse[d] <- sqrt((1/nrow(meter2.train)) * sum((log(meter2.gam.pred+1) - log(meter2.validation$meter_reading+1))^2))
}

min(rmlse)
which.min(rmlse)
plot(1:10, rmlse, type = 'b')

#smoothing splines
rmlse <- c(rep(NA, 10))
for (d in 1:10) {
  meter2.gam <- gam(meter_reading ~ primary_use + s(square_feet, df=d) + s(sea_level_pressure, df=d) + building_id + site_id, data = meter2.training)
  
  meter2.gam.pred <- predict(meter2.gam, newdata = meter2.validation)
  meter2.gam.pred[meter2.gam.pred < 0] = 0
  rmlse[d] <- sqrt((1/nrow(meter2.train)) * sum((log(meter2.gam.pred+1) - log(meter2.validation$meter_reading+1))^2))
}

min(rmlse)
which.min(rmlse)
plot(1:10, rmlse, type = 'b')

###########
#GBM
###########
set.seed(488)
#all features
rmlse <- c(rep(NA, 3))
lambdas <- c(0.001, 0.01, 0.1)
i <- 1
for (lambda in lambdas){
  meter2.gbm <- gbm(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + precip_depth_1_hr + sea_level_pressure + wind_direction + wind_speed + cloud_coverage + building_id + site_id, data = meter2.training, shrinkage = lambda, cv.folds = 10)
  meter2.gbm.pred <- predict(meter2.gbm, newdata = meter2.validation)
  meter2.gbm.pred[meter2.gbm.pred < 0] = 0
  rmlse[i] <- sqrt((1/nrow(meter2.train)) * sum((log(meter2.gbm.pred+1) - log(meter2.validation$meter_reading+1))^2))
  i <- i+1
}

min(rmlse)
lambdas[which.min(rmlse)]
plot(lambdas, rmlse, type = 'b')

#subset of features based on BIC
rmlse <- c(rep(NA, 3))
lambdas <- c(0.001, 0.01, 0.1)
i <- 1
for (lambda in lambdas){
  meter2.gbm <- gbm(meter_reading ~ primary_use + square_feet  + sea_level_pressure + building_id + site_id, data = meter2.training, shrinkage = lambda, cv.folds = 10)
  meter2.gbm.pred <- predict(meter2.gbm, newdata = meter2.validation)
  meter2.gbm.pred[meter2.gbm.pred < 0] = 0
  rmlse[i] <- sqrt((1/nrow(meter2.train)) * sum((log(meter2.gbm.pred+1) - log(meter2.validation$meter_reading+1))^2))
  i <- i+1
}

min(rmlse)
lambdas[which.min(rmlse)]
plot(lambdas, rmlse, type = 'b')
##############################################################################################################

#METER THREE: HOT WATER
meter3.train <- subset(train, meter == 3)

#more EDA
ggplot(aes(x=dew_temperature, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=air_temperature, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=cloud_coverage, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=primary_use, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=square_feet, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=wind_speed, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=wind_direction, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=sea_level_pressure, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=precip_depth_1_hr, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=building_id, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()
ggplot(aes(x=site_id, y=meter_reading), data=meter3.train) + geom_point() + geom_smooth()

#Backward Selection
models <- regsubsets(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + precip_depth_1_hr + sea_level_pressure + wind_direction + wind_speed + cloud_coverage + building_id + site_id, data = meter3.train, method = 'backward', nvmax = 11)
best.features <- summary(models)
data.frame(num_features = which.min(best.features$bic), BIC = min(best.features$bic))
best.features

#split into training and validation sets
meter3.training <- meter3.train[meter3.train$timestamp >= '2016-01-31 18:00:00',]
meter3.validation <- meter3.train[meter3.train$timestamp < '2016-01-31 18:00:00',]

###########
#GAMS
###########
rmlse <- c(rep(NA, 10))
for (d in 1:10) {
  meter3.gam <- gam(meter_reading ~ primary_use + ns(square_feet, df=d) + ns(air_temperature, df=d) + ns(dew_temperature, df=d) + building_id + site_id, data = meter3.training)
  
  meter3.gam.pred <- predict(meter3.gam, newdata = meter3.validation)
  meter3.gam.pred[meter3.gam.pred < 0] = 0
  rmlse[d] <- sqrt((1/nrow(meter3.train)) * sum((log(meter3.gam.pred+1) - log(meter3.validation$meter_reading+1))^2))
}

min(rmlse)
which.min(rmlse)
plot(1:10, rmlse, type = 'b')

#smoothing splines
rmlse <- c(rep(NA, 10))
for (d in 1:10) {
  meter3.gam <- gam(meter_reading ~ primary_use + s(square_feet, df=d) + s(air_temperature, df=d) + s(dew_temperature, df=d) + building_id + site_id, data = meter3.training)
  
  meter3.gam.pred <- predict(meter3.gam, newdata = meter3.validation)
  meter3.gam.pred[meter3.gam.pred < 0] = 0
  rmlse[d] <- sqrt((1/nrow(meter3.train)) * sum((log(meter3.gam.pred+1) - log(meter3.validation$meter_reading+1))^2))
}

min(rmlse)
which.min(rmlse)
plot(1:10, rmlse, type = 'b')

###########
#GBM
###########
set.seed(488)
#all features
rmlse <- c(rep(NA, 3))
lambdas <- c(0.001, 0.01, 0.1)
i <- 1
for (lambda in lambdas){
  meter3.gbm <- gbm(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + precip_depth_1_hr + sea_level_pressure + wind_direction + wind_speed + cloud_coverage + building_id + site_id, data = meter3.training, shrinkage = lambda, cv.folds = 10)
  meter3.gbm.pred <- predict(meter3.gbm, newdata = meter3.validation)
  meter3.gbm.pred[meter3.gbm.pred < 0] = 0
  rmlse[i] <- sqrt((1/nrow(meter3.train)) * sum((log(meter3.gbm.pred+1) - log(meter3.validation$meter_reading+1))^2))
  i <- i+1
}

min(rmlse)
lambdas[which.min(rmlse)]
plot(lambdas, rmlse, type = 'b')

#subset of features based on BIC
rmlse <- c(rep(NA, 3))
lambdas <- c(0.001, 0.01, 0.1)
i <- 1
for (lambda in lambdas){
  meter3.gbm <- gbm(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + building_id + site_id, data = meter3.training, shrinkage = lambda, cv.folds = 10)
  meter3.gbm.pred <- predict(meter3.gbm, newdata = meter3.validation)
  meter3.gbm.pred[meter3.gbm.pred < 0] = 0
  rmlse[i] <- sqrt((1/nrow(meter3.train)) * sum((log(meter3.gbm.pred+1) - log(meter3.validation$meter_reading+1))^2))
  i <- i+1
}

min(rmlse)
lambdas[which.min(rmlse)]
plot(lambdas, rmlse, type = 'b')
##############################################################################################################

#PREDICTING ON THE TEST SET
test <- get(load('ashrae_test.RData'))
test <- na_mean(test, option = 'median')
test$primary_use <- as.numeric(test$primary_use)
meter0.test <- subset(test, meter == 0)
meter1.test <- subset(test, meter == 1)
meter2.test <- subset(test, meter == 2)
meter3.test <- subset(test, meter == 3)

#fitting best models
meter0.model <- gbm(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + wind_direction + cloud_coverage + building_id + site_id, data = meter0.training, shrinkage = 0.1, cv.folds = 10)
meter1.model <- gam(meter_reading ~ primary_use + s(square_feet, df=2) + s(dew_temperature, df=2) + s(sea_level_pressure, df=2) + site_id, data = meter1.train)
meter2.model <- gbm(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + precip_depth_1_hr + sea_level_pressure + wind_direction + wind_speed + cloud_coverage + building_id + site_id, data = meter2.train, shrinkage = 0.1, cv.folds = 10)
meter3.model <- gbm(meter_reading ~ primary_use + square_feet + air_temperature + dew_temperature + building_id + site_id, data = meter3.train, shrinkage = 0.1, cv.folds = 10)

#predicting on test sets
pred0 <- predict(meter0.model, newdata = meter0.test)
pred0[pred0 < 0] = 0
pred1 <- predict(meter1.model, newdata = meter1.test)
pred1[pred1 < 0] = 0
pred2 <- predict(meter2.model, newdata = meter2.test)
pred2[pred2 < 0] = 0
pred3 <- predict(meter3.model, newdata = meter3.test)
pred3[pred3 < 0] = 0

#submission file
bho.sub <- data.frame(row_id = c(meter0.test$row_id, meter1.test$row_id, meter2.test$row_id, meter3.test$row_id), meter_reading = c(pred0,pred1,pred2,pred3))
write.csv(bho.sub, file = "bho_submission.csv", row.names = FALSE)