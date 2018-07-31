library(dplyr)
library(data.table)
library(stats)
library(tibble)
#setwd("E:/Exploratory Data Analysis/Assignments - 3B/Assignment-1/parliament")
#par <- fread("parliament.csv")

#' @title func
#' @description  cdcvdhc
#' @param df
#' @return res
#' @export split.vectors
split.vectors <- function(df)
{
  num_ind <- which(sapply(df,class) %in% c("integer","numeric","double"))
  num.vect <- colnames(df)[num_ind]
  char_ind <- which(sapply(df,class) %in% c("character","factor"))
  char.vect <- colnames(df)[char_ind]
  res <- list("Numeric_Columns" = num.vect,"Character_Columns" = char.vect)
  return(res)
}

#' @title func
#' @description  cdcvdhc
#' @param df
#' @return res
#' @export numeric_summary

numeric_summary <- function(df)
{
  num_ind <- which(sapply(df,class) %in% c("integer","numeric","double"))
  num_df <- df[,num_ind]
  min <- sapply(num_df,min,na.rm = T)
  max <- sapply(num_df,max,na.rm = T)
  per.25 <- sapply(num_df,quantile,na.rm = T)[2,]
  per.50 <- sapply(num_df,quantile,na.rm = T)[3,]
  per.75 <- sapply(num_df,quantile,na.rm = T)[4,]
  mean <- sapply(num_df,mean,na.rm = T)
  median <- sapply(num_df,median,na.rm = T)
  variance <- sapply(num_df,function(x){options(scipen = 10)
    var(x,na.rm = T)})
  sd <- sapply(num_df,sd,na.rm = T)
  per.na <- round(colSums(is.na(num_df))/nrow(df),3)*100
  count.na <- colSums(is.na(num_df))
  per.out <- sapply(num_df,function(x){
    round(length(boxplot.stats(x)$out)/nrow(num_df),3)*100})
  res <- data.frame(min,max,per.25,per.50,per.75,mean,median,variance,sd,count.na,per.na,per.out)
  res <- as.data.frame(t(res),stringsAsFactors = F)
  res <- rownames_to_column(res,var = "STATISTIC")
  return(res)
}


#' @title func
#' @description  cdcvdhc
#' @param df
#' @return res
#' @export factor_summary

factor_summary <- function(df)
{
  fact.char_ind <- which(sapply(df,class) %in% c("character","factor"))
  fact_df <- df[,fact.char_ind]
  uniq.levels <- sapply(fact_df,function(x){length(levels(as.factor(x)))})
  mode <- sapply(fact_df,function(x){levels(as.factor(x))[which(table(x) == max(table(x)))]})
  mode.freq <- sapply(fact_df,function(x){max(table(x))})
  per.levels <- sapply(fact_df,function(x){
    round(sum(cumsum(sort(table(x)/sum(table(x))*100,decreasing = T)) < 80)/length(levels(as.factor(x)))*100,3)})
  res <- data.frame(uniq.levels,mode,mode.freq,per.levels)
  res <- as.data.frame(t(res),stringsAsFactors = F)
  res <- rownames_to_column(res,var = "STATISTIC")
  return(res)
}
