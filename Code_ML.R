#' point of this excercise is to check if it is possible to
#' predict secondary strucutre of protiens using given
#' amino acid sequence enviroment and machine learning

library(dplyr)
library(tidyr)
library(mlr)
library(Peptides)
#library(Biostrings)

## secondary structure database can be downloaded from PDB under this address:
## https://cdn.rcsb.org/etl/kabschSander/ss.txt.gz
## this is "A FASTA-formatted file ("ss.txt.gz")
## with sequences and secondary structure information generated using DSSP"
## more about DSSP: https://en.wikipedia.org/wiki/DSSP_(hydrogen_bond_estimation_algorithm)

# first set path
path <- "~/Path/to/this_directory/"
setwd(path)

# load the database
ssdata <- Biostrings::readAAStringSet("ss.txt", format="fasta") 

# quick look at data
ssdata %>% head()
length(ssdata)
str(ssdata)

# make data frame
names <- names(ssdata)
sequences <- paste(ssdata)
df <- data.frame(names, sequences)
df$sequences <- as.character(df$sequences)

# check data frame
summarizeColumns(df)

# extract info about sequence vs secondary structure (ss)
df <- df %>% separate(names, into = c("names","type"), sep = ":se")

# make seq and ss into columns of df
df_seq <- df %>% filter(type == "quence")
df_ss <- df %>% filter(type ==  "cstr")
df_seq$ss <- df_ss$sequences
df_seq$type <- NULL
df <- df_seq

### DATA DEDUPLICATION
# as some proteins are overrepresented in PDB it is ok to do deduplication at this stage
# find identical sequences and purge them, optionally check if ss are also the same
nrow(df) # all rows
nrow(distinct(df, ss, sequences)) # unique rows
df <- distinct(df, ss, sequences, .keep_all = TRUE)

# higlight gaps by "-"
df$ss <-  gsub(' ', '-', df$ss)

###
#  DATA SUBSETTING
###
# N.B. for reproducible results set seed

fraction_ <- as.integer(nrow(df)*0.1) # we will take just 10% of data for faster computation
df_train <- sample_n(df, fraction_)

# now substract sample from full data
df_no_s <- anti_join(df, df_train)

# then draw another sample and you know that train and test data will not overlap
df_test <- sample_n(df_no_s, 100)


### Functions definitions

#' split_seq this function will take data single row of data frame and split "sequences" column and "ss" columns
#' info single characters and return new data frame with seq and ss columns preserving names information
#' @param df - which row should the function take
#' @param x - from which data frame
split_seq <- function(df, x) {
  a <- df[x,] # take n-th row of data frame
  y <- (strsplit(a$sequences, "")) # split seq
  z <- (strsplit(a$ss, "")) # split ss
  new_df <- data.frame(name = a$names,
                       seq = y[[1]],
                       ss = z[[1]])
  return(new_df)
}

#' add_context this function adds columns with up to 5 aa's before and after given amino acid (aa)
#' @param n - size of margin (number of residues on each site)
#' @param df - data frame
#' @return data frame
add_context <- function(df, n){
  # add lagging
  for (i in seq(1,n)){
    name_X <- paste0("m", i)
    df[[name_X]] <- with(df, as.character(lag(seq, n=i)))
  }
  # add leading
  for (i in seq(1,n)){
    name_X <- paste0("p", i)
    df[[name_X]] <- with(df, as.character(lead(seq, n=i)))
  }
  #dealing with NAs (N.B. as.character is required for this)
  df[is.na(df)]<-"-" 
  return(df)
} # function

### /Functions definitions

###
#  Data processing - this will take some CPU time
###
df_list <- list()
for (i in seq(1,nrow(df_train))){
  # split sequence, then add neighbouring amino acids
  temp_ <- split_seq(x=i, df=df) %>% add_context(., 3)
  # addto list of data frames
  df_list[[length(df_list)+1]] <- temp_
}

# make one data frame from list of df's
train <-  bind_rows(df_list) # error due to different factor levels so they will be coerced to character

# loose names column
train$name <- NULL

# check Data Redundancy - data deduplication
# N.B. you might not want to deduplicate data here if you would like to retain data proportions.
train <- distinct(train)

summarizeColumns(train)

# change columns types to factors
train <- train %>% mutate_all(.funs = function(x){if(is.character(x)) as.factor(x)})


### N.B. we could add another properties here
# for example:
amino_acid <- "a"
as.data.frame(aaComp(amino_acid))


###
#  Machine learning
###

#N.B. cannot check correlations in data as all variables are factors

# make train task
trainTask <- makeClassifTask(data = train, target = "ss")

# Data normalisation - not needed as only factors here - TO DO: comment this out
trainTask <- normalizeFeatures(trainTask,method = "standardize")

### Feature importance estimation
# this is to load java on MacOS and avoid errors Java is needed for FSelector and rJava
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_25.jdk/Contents/Home/jre/lib/server/libjvm.dylib')

#Feature importance
im_feat <- generateFilterValuesData(trainTask, 
              method = c("information.gain", "chi.squared", "kruskal.test")) # if out of memory try just kruskal.test

#plot results
plotFilterValues(im_feat,n.show = 20)

### here we could get rid of not so important features

# choose an algorithm for now we will proceed with:
### Random forest

#get parameters
getParamSet("classif.randomForest")

learner_Forest <- makeLearner("classif.randomForest", predict.type = "response", par.vals = list(ntree = 20, mtry = 3))

# train model
model_Forest <- train(learner_Forest, trainTask)

# but maybe we could check for best parameters
# tune parameters - it is a lot of CPU time

#grid search to find best parameters
rf_param <- makeParamSet(
  makeIntegerParam("ntree",lower = 20, upper = 500),
  makeIntegerParam("mtry", lower = 3, upper = 10),
  makeIntegerParam("nodesize", lower = 10, upper = 50))

# cross validation
set_cv <- makeResampleDesc("CV",iters = 3L)

#do a grid search for best parameters
rf_control <- makeTuneControlGrid()

stune <- tuneParams(learner = learner_Forest, resampling = set_cv, 
                    task = trainTask, par.set = rf_param, control = gscontrol, measures = acc)

# use computed parameters
L_Forest <- setHyperPars(learner_Forest, par.vals = stune$x)

model_Forest <- train(L_Forest, trainTask)
getLearnerModel(model_Forest)

### optionally
# save trained model
saveRDS(model_Forest, "model_Forest.rds")

# and this is how you can load a model
model_Forest <- readRDS("model_Forest.rds")

###
# MAKE PREDICTION
###
test_seq <- df_test[1,] # take first seuqnce
# split sequence, then add neighbouring amino acids
test_seq$name <- NULL
test_seq <- split_seq(x=1, df=test_seq) %>% add_context(., 3)
# change columns types to factors
test_seq <- test_seq %>% mutate_all(.funs = function(x){as.factor(x)})

# In order to get rid of nasty bug and equalise factor levels in train and test data
common <- intersect(names(train), names(test_seq))
for (i in common) {
  if (class(train[[i]]) == "factor") {levels(test_seq[[i]]) <- levels(train[[i]])}
  }

# make prediction
Prediction <- predict(model_Forest, newdata = test_seq)

# return result
Sset <- Biostrings::AAStringSet(c(paste(test_seq$seq, collapse = ''),
                      paste(Prediction$data$response, collapse = '')))
names(Sset) = c("sequence", "prediced ss")
print(Sset)

# compare prediction to real values
P_data <- Prediction$data
P_data <- P_data %>% mutate(match = as.numeric(truth==response))

# calculate Q3 accuracy
Q3 <- sum(P_data$match)/length(P_data$match)
print(paste("Q3 accuracy:",round(Q3*100, digits = 1),"%"))


