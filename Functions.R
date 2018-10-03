library(dplyr)
library(tidyr)
library(Peptides)
library(parallel)
library(pbmcapply)

# detect cores
no_cores <- round(detectCores()*0.8) # do not take all cores admin will be not happy

# make list of amina acids (aa) properties names and suffixes
b <- as.data.frame(aaComp("-"))
rownames(b)
suffixes <- substr(rownames(b), 1, 2)

### TO OD: try to replace for loops with apply families maybe multithread

### Functions definitions

#' split_seq this function will take data single row of data frame and split "sequences" column and "ss" columns
#' info single characters and return new data frame with seq and ss columns preserving names information
#' @param df - which row should the function take
#' @param x - from which data frame
#' @return data frame
#' @export
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
#' @export
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

#' this function will return properties of each amino acid in given column
#' by adding new column with apropriate suffixes to given data frame
#'
#' @param df - data frame
#' @param col - column number
#' @return data frame
#' @export
for_column <- function(col, df){
name_c <- names(df)[col] # - column name
new_names <- paste(name_c, suffixes, sep="_")
dat_f <- data.frame(matrix(NA,1,length(new_names))) # initaite data frame with 9 columns
dat_f <- dat_f[-1,] # dele first row of new data frame

fac <- df[,col] # sequence of column

##old
# for (i in seq(1,length(fac))){
#   aa <- fac[i] # single aa
#   b <- as.data.frame(aaComp(aa)) # get properties using Peptides library
#   c <- b[,-2]
#   dat_f <- rbind(dat_f, c)
#   #cat("\r",paste("processing column:", name_c, round(i*100/length(fac)), "% complete")) # progress
# }

# new with lapply
dat_f <- lapply(seq(1,length(fac)), function(fac=fac){
  aa <- fac[i] # single aa
  b <- as.data.frame(aaComp(aa)) # get properties using Peptides library
  c <- b[,-2]
  return(c)
})
dat_f <- dat_f %>% as.list() %>% do.call(rbind,.) %>% as.data.frame()

colnames(dat_f) <- new_names # add names to data frame
# add column with seq - name_c number 2
dat_f[,name_c] <- df[,name_c]
# add our new data to old data frame
expanded_df <- semi_join(dat_f, df)
return(expanded_df)
} # for_column function 


#' this fucntions will add amino acids properties (as separate columns) for each column in givrn dataframe
#' @param df data frame
#' @return data frame
#' @export
add_aa_properties <- function(df){
  n <- names(df)
  s <- seq(1,length(n))[!(n == "ss")] # take all column indexes but "ss"
  # we assume that names column is already gone at this stage
  df_list <- list() # make list
  
  ##old code slow
  # for (i in s){df_list[[length(df_list)+1]] <- for_column(df, i)}
  ## better - lapply
  #df_list <- lapply(s, for_column, df)
  
  # parallel lapply with progress bar
  df_list <- pbmclapply(s, for_column, df=df, mc.cores = no_cores)
  
  df_result <-  bind_cols(df_list) # combine all columns together
  # change 0, 1 into factors
  df_result <- df_result %>% mutate_all(.funs = function(x){as.factor(x)})
  return(df_result)
}

### /Functions definitions


