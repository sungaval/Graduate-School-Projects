# Machine Learning final project
# Valentinas Sungaila

# Loading all my packages

library(SASxport)
library(plyr)
library(earth)     # for fitting MARS models
library(caret)
library(ISLR)
library(dplyr)
library(tidyverse)
library(ggplot2)



# Data cleaning 


#set the working directory

setwd("C:/Users/sunga/Desktop/machine learning final project")



# Read in an xpt file

lookup.xport("LLCP2018.xpt")
brfss_Survey_Data_2018  <- read.xport("LLCP2018.xpt")
# how to open .xpt files found here https://www.phusewiki.org/wiki/index.php?title=Open_XPT_File_with_R



# Limit my data a little to be more focused


# let me focus on people who only live in a private residence, (By private residence, we mean someplace like a house or apartment.)
brfss_Survey_Data_2018 = subset(brfss_Survey_Data_2018, PVTRESD3 == 1)
# focus on people currently living in the state they took the survey at
brfss_Survey_Data_2018 = subset(brfss_Survey_Data_2018, CSTATE1 == 1)
# focus on people who are either male or female at birth insterad of people who refused to answer
brfss_Survey_Data_2018 = subset(brfss_Survey_Data_2018, SEX1 <= 2)
# make it so the weight is only in pounds
brfss_Survey_Data_2018 = subset(brfss_Survey_Data_2018, WEIGHT2 <= 0999)
# make it so height is only in feet and inches
brfss_Survey_Data_2018 = subset(brfss_Survey_Data_2018, HEIGHT3 <= 711)
# replace the number of children from 88 for 0 children to 0 for 0 children
brfss_Survey_Data_2018$CHILDREN[brfss_Survey_Data_2018$CHILDREN == 88] = 0



# cleaning up data for variables I do not need because they are duplicates or some variation of the same question

brfss_Survey_Data_2018 <- brfss_Survey_Data_2018[ -c(2:219) ]
brfss_Survey_Data_2018 <- brfss_Survey_Data_2018[ -c(13, 42:43, 47:57) ]
brfss_Survey_Data_2018 <- brfss_Survey_Data_2018[ -c(1,3) ]


# seems to have worked properly but now I have some NA means for some of my variables which means I will need to remove them
# keeping only the columns without na values
# found https://r.789695.n4.nabble.com/how-to-delete-columns-with-NA-values-td1839902.html
# 4th comment

brfss_Survey_Data_2018 = brfss_Survey_Data_2018[,colSums(is.na(brfss_Survey_Data_2018)) == 0]

# make my response which is adults with good or better health into a binary 1,2 variable 

brfss_Survey_Data_2018 = subset(brfss_Survey_Data_2018, X.RFHLTH <= 2)
# 1 means people are in good or better health
# 2 means people are in fair or poor health

table(brfss_Survey_Data_2018$X.RFHLTH)
# looks good
# finished with data cleaning


##########################################################################################
# Since I have a binary response variable if I want to use a binary GLM in my MARS function to check for model fit I need to make my response 
# which is adults with good or better health into a binary 0,1 variable so I can see if there is any effect if I specify a distribution for my response

brfss_Survey_Data_2018_binomial = mutate(brfss_Survey_Data_2018, X.RFHLTH_binary = X.RFHLTH - 1)


# function to move columns around
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

# function creator found here
# https://stackoverflow.com/questions/3369959/moving-columns-within-a-data-frame-without-retyping/18540144#18540144

# The basic options are:
#   
#  first
#  last
#  before
#  after
# Compounded moves are separated by a semicolon.

# example
# moveme(names(df), "g first")
# moveme(names(df), "g first; a last; e before c")


# move my variable of interest which is X.RFHLTH first then X.RFHLTH_binary


brfss_Survey_Data_2018_binomial = brfss_Survey_Data_2018_binomial[moveme(names(brfss_Survey_Data_2018_binomial), 'X.RFHLTH_binary after X.RFHLTH')]

# remove the unnecessary column
brfss_Survey_Data_2018_binomial <- brfss_Survey_Data_2018_binomial[ -1 ]


###-----------------------------------------------------------------------------------------------------------------------------------
# Exploratory analysis
ggplot(data = brfss_Survey_Data_2018, aes(x = as.factor(brfss_Survey_Data_2018$X.INCOMG),
                                          fill = as.factor(brfss_Survey_Data_2018$X.RFHLTH))) +
  geom_bar(stat="count") +
  labs(x = "Income", title = "Income of people with good health vs bad health",
       colour = "Health",
       caption = "(1 = less than 15k, 2 = 15-25k, 3 = 25-35k, 4 = 35-50k, 5 = 50k and above, 9 = refuse to answer )") +
  scale_fill_discrete(name="Health Status",
                      labels=c("Good", "Bad"))



# Histogram of income
hist(brfss_Survey_Data_2018$X.INCOMG,freq = F)




ggplot(data = brfss_Survey_Data_2018, aes(x = as.factor(brfss_Survey_Data_2018$X.SMOKER3),
                                          fill = as.factor(brfss_Survey_Data_2018$X.RFHLTH))) +
  geom_bar(stat="count") +
  labs(x = "Smoking Status", title = "Smokers With Good Health vs Smokers with Poor Health",
       colour = "Health",
       caption = "(1 = smoke every day, 2 = smoke sometimes, 3 = former smoker, 4 = never smokes, 9 = refuse to answer )") +
  scale_fill_discrete(name="Health Status",
                      labels=c("Good", "Bad"))


####-----------------------------------------------------------------------------------------------------------------------------------

# First, split the data into a train set and test set for the normal data.  
# 75% Train
# 25% Test

set.seed(1)

n = length(brfss_Survey_Data_2018$X.RFHLTH) # number of observations

tr = sample(1:n, size = floor(n*.75))
t = (1:n)[-tr]

train_set = data.frame(brfss_Survey_Data_2018[tr, ])
test_set = data.frame(brfss_Survey_Data_2018[t, ])

train_n = length(tr)
test_n = length(t)

###

# my data with a 0,1 response variable  into a train set and test set.  
# 75% Train
# 25% Test


n_binomial = length(brfss_Survey_Data_2018_binomial$X.RFHLTH_binary) # number of observations

tr_binomial = sample(1:n_binomial, size = floor(n_binomial*.75))
t_binomial = (1:n_binomial)[-tr_binomial]

train_set_binomial = data.frame(brfss_Survey_Data_2018_binomial[tr_binomial, ])
test_set_binomial = data.frame(brfss_Survey_Data_2018_binomial[t_binomial, ])

train_n_binomial = length(tr_binomial)
test_n_binomial = length(t_binomial)


###------------------------------------------------------------------------------------------------------------

# running the mars model where the number of knots is decided by R

mars1_train = earth(X.RFHLTH ~ ., data = train_set)
print(mars1_train)
summary(mars1_train)



plot(mars1_train, which = 1)

# The model selection plot that shows us the Generalized cross-validation (GCV) R^2 
#(left-hand y-axis and solid black line) based on the number of terms retained in the model 
#(x-axis) which are constructed from a certain number of original predictors (right-hand y-axis). 
#The vertical dashed lined at 14 tells us the optimal number of terms retained where marginal increases 
#in GCV R^2 are less than 0.001.



# summarize the importance of input variables
evimp(mars1_train)
#Show the variables that are important


## make predictions on training data
predictions <- predict(mars1_train, train_set)

#summarize accuracy of training 
mse <- mean((predictions - train_set$X.RFHLTH)^2)
print(mse)


## make predictions on test data
predictions2 <- predict(mars1_train, test_set)


#summarize accuracy of test
mse2 <- mean((predictions2 - test_set$X.RFHLTH)^2)
print(mse2)


###----------------------------------------------------------------------------


# running the mars model where the number of knots is decided by R but we use 10-fold cross validation

mars2_train = earth(X.RFHLTH ~ ., data = train_set, nfold = 10, ncross = 10, trace = .5, pmethod="cv")
print(mars2_train)
summary(mars2_train)



plot(mars2_train, which = 1)
plotd(mars2_train)
# The model selection plot that shows us the Generalized cross-validation (GCV) R^2 
#(left-hand y-axis and solid black line) based on the number of terms retained in the model 
#(x-axis) which are constructed from a certain number of original predictors (right-hand y-axis). 
#The vertical dashed lined at 14 tells us the optimal number of terms retained where marginal increases 
#in GCV R^2 are less than 0.001.



# summarize the importance of input variables
evimp(mars2_train)
#Show the variables that are important


## make predictions on training data
predictions3 <- predict(mars2_train, train_set)

#summarize accuracy of training 
mse3 <- mean((predictions3 - train_set$X.RFHLTH)^2)
print(mse3)


## make predictions on test data
predictions4 <- predict(mars2_train, test_set)


#summarize accuracy of test
mse4 <- mean((predictions4 - test_set$X.RFHLTH)^2)
print(mse4)

# still get the same thing

# compare MARS models
plot.earth.models(list(mars2_train,mars1_train))

#----------------------------------------------------------------------------------------------------------
# let me specify a binomial family glm for the model since my response is 1 and 2 to see if that makes any difference

# running the mars model where the number of knots is decided by R

mars3_train = earth(X.RFHLTH_binary ~ ., data = train_set_binomial, glm=list(family=binomial))
print(mars3_train)
summary(mars3_train)



plot(mars3_train, which = 1)

# The model selection plot that shows us the Generalized cross-validation (GCV) R^2 
#(left-hand y-axis and solid black line) based on the number of terms retained in the model 
#(x-axis) which are constructed from a certain number of original predictors (right-hand y-axis). 
#The vertical dashed lined at 14 tells us the optimal number of terms retained where marginal increases 
#in GCV R^2 are less than 0.001.



# summarize the importance of input variables
evimp(mars3_train)
#Show the variables that are important


## make predictions on training data
predictions5 <- predict(mars3_train, train_set_binomial)

#summarize accuracy of training 
mse5 <- mean((predictions5 - train_set_binomial$X.RFHLTH_binary)^2)
print(mse5)


## make predictions on test data
predictions6 <- predict(mars3_train, test_set_binomial)


#summarize accuracy of test
mse6 <- mean((predictions6 - test_set_binomial$X.RFHLTH_binary)^2)
print(mse6)


###----------------------------------------------------------------------------


# running the mars model where the number of knots is decided by R but we use 10-fold cross validation
memory.limit(size=51000)
mars4_train = earth(X.RFHLTH_binary ~ ., data = train_set_binomial, nfold = 10, trace = .5, pmethod="cv", glm=list(family=binomial))
print(mars4_train)
summary(mars4_train)



plot(mars4_train, which = 1)
plot(mars4_train)
# The model selection plot that shows us the Generalized cross-validation (GCV) R^2 
#(left-hand y-axis and solid black line) based on the number of terms retained in the model 
#(x-axis) which are constructed from a certain number of original predictors (right-hand y-axis). 
#The vertical dashed lined at 14 tells us the optimal number of terms retained where marginal increases 
#in GCV R^2 are less than 0.001.



# summarize the importance of input variables
evimp(mars4_train)
#Show the variables that are important


## make predictions on training data
predictions7 <- predict(mars4_train, train_set_binomial)

#summarize accuracy of training 
mse7 <- mean((predictions7 - train_set_binomial$X.RFHLTH_binary)^2)
print(mse7)


## make predictions on test data
predictions8 <- predict(mars4_train, test_set_binomial)


#summarize accuracy of test
mse8 <- mean((predictions8 - test_set_binomial$X.RFHLTH_binary)^2)
print(mse8)

# still get the same thing
plot.earth.models(list(mars3_train,mars4_train))

#---------------------------------------------------------------------------------------------



# foud how to  train MARS models with different values of nk, thresh and span (including minspan and endspan)
# nk = Maximum number of model terms before pruning
# thresh = Forward stepping threshold
# span = minspan: Minimum number of observations between knots. endspan: Minimum number of observations before the first and after the final knot


# A simulation can be conducted to show how different values of nk, thresh, minspan and endspan affects the model-training process
# https://blog.zenggyu.com/en/post/2018-06-16/multivariate-adaptive-regression-splines-in-a-nutshell/


results <- crossing(nk = c(5, 10, 20),
                    thresh = c(0, 0.01, 0.1),
                    span = c(1, 5, 10)) %>%
  pmap(function(nk, thresh, span, train_set) {
    fit <- earth(X.RFHLTH ~ ., data = train_set,
                 degree = 1, nprune = NULL,
                 nk = nk, thresh = thresh, minspan = span, endspan = span)
    mutate(train_set, nk = nk, thresh = thresh, span = span,
           y_predicted = predict(fit)[,1], p = length(coef(fit)))
  }, train_set = train_set) %>%
  bind_rows()


# the number of terms included in each model with different parameters (the output is attached below)
results %>%
  select(nk, thresh, span, p) %>%
  distinct() %>%
  mutate(nk = as.factor(nk) %>% fct_relabel(function(x) {paste0("nk=", x)})) %>%
  split(.$nk) %>%
  map(function(x) {
    x %>%
      select(-nk) %>%
      spread(key = span, value = p, sep = "=") %>%
      as.data.frame() %>%
      `rownames<-`(paste0("thresh=", .$thresh))  %>%
      select(-thresh)
  })

3##############################################
# creating segments to show for MARS model

plot(0, 0, col = "white",
     main = "Physical Activity vs Health MARS prediction",
     xlab = "Physical Activity",
     ylab = "Reported Physical Health",
     xlim = c(0, 10), 
     ylim = c(0.8, 1.38))     

segments(x0 = c(0,2), 
         y0 = c(1.36428068,1.134429), 
         x1 = c(2,10), 
         y1 = c(1.134429,1.196296)) 




###################################---------------------------------------------------------------------------------

# lets model the data with a logistic regression and see how it compares.

logit = glm(X.RFHLTH_binary ~ ., data = train_set_binomial, family=binomial)

summary(logit)

## make predictions on training data
predictions9 <- predict(logit, train_set_binomial)

#summarize accuracy of training 
mse9 <- mean((predictions9 - train_set_binomial$X.RFHLTH_binary)^2)
print(mse9)


## make predictions on test data
predictions10 <- predict(logit, test_set_binomial)


#summarize accuracy of test
mse10 <- mean((predictions10 - test_set_binomial$X.RFHLTH_binary)^2)
print(mse10)



