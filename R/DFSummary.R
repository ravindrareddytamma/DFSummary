library(dplyr)
library(data.table)
library(stats)
library(tibble)
library(tidyr)
library(ggplot2)

#' @title func
#' @description  Returns list of Numeric and Character Columns
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
#' @description  Detects the Numerical Columns in the DataFrame and provides the Summary Statistics
#' @param df
#' @return res
#' @export numeric_summary

numeric_summary <- function(df)
{
  num_ind <- which(sapply(df,class) %in% c("integer","numeric","double"))
  num_df <- subset.data.frame(df,select = names(df)[num_ind])
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
  res <- tibble::rownames_to_column(res,var = "STATISTIC")
  return(res)
}


#' @title func
#' @description  Detects the Character and Factor Columns in the DataFrame and provides the Summary Statistics
#' @param df
#' @return res
#' @export factor_summary

factor_summary <- function(df) 
{
  fact.char_ind <- which(sapply(df, class) %in% c("character","factor"))
  fact_df <- subset(df, select = colnames(df)[fact.char_ind])
  ## Unique Levels
  uniq.levels <- sapply(fact_df, function(x) {
    length(levels(as.factor(x)))
  })
  
  ## Mode of the Data Column
  mode <- sapply(fact_df, function(x) {
    names(table(x)[which.max(table(x))])
  })
  mode.freq <- sapply(fact_df, function(x) {
    table(x)[which.max(table(x))]
  })
  
  ## Percentage of Most repeating Level
  per.levels <- sapply(fact_df, function(x) {
    round(sum(cumsum(sort(table(x)/sum(table(x)) * 100, decreasing = T)) < 
                80)/length(levels(as.factor(x))) * 100, 3)
  })
  res <- data.frame(uniq.levels, mode, mode.freq, per.levels)
  res <- as.data.frame(t(res), stringsAsFactors = F)
  res <- tibble::rownames_to_column(res, var = "STATISTIC")
  return(res)
}

#' @title func
#' @description  Sturge Formula to convert the Numerical column into appropriate Class Intervals. 
#' @param df
#' @return freq.table
#' @export factor_summary


sturge <- function (vect, bins = 0,cat.col = TRUE) 
{
  if (class(vect) == "character") 
    return(transform(table(vect)))
  n <- length(vect)
  if (n == 1) 
    return(table(vect))
  low <- round(min(vect, na.rm = TRUE))
  high <- round(max(vect, na.rm = TRUE))
  k <- round(log2(n))
  if (bins != 0) {
    width <- round((high - low)/(bins - 1))
  }else {
    width <- round((high - low)/k)
  }
  if (width > 0) {
    bins <- seq(low, high + width, width)
  }else {
    stop("Width of Interval is not Correct or NA values in DataSet")
  }
  interval <- cut(vect, bins, dig.lab = 5)
  freq.table <- transform(table(interval))
  ls = list("Freq.Table" = freq.table,"Vect" = interval)
  return(ls)
}

#' @title func
#' @description  Detects the Numerical Columns and performs Pairwise t-test for all Numerical Columns in DataFrame  
#' @param df
#' @return res
#' @export factor_summary


DF.ttest <- function(data,num_cols = c())
{
  require(plyr,quietly = T)
  data <- data.frame(data)
  if(length(num_cols) > 0)
  {
    num_cols <- names(data)[num_cols]
  }else{
    num_cols <- split.vectors(data)$Numeric_Columns
  }
  
  data <- data[,num_cols]
  combos <- combn(ncol(data),2)
  res <- data.frame()
  for(i in 1:ncol(combos))
  {
    result <- t.test(data[,combos[1,i]],data[,combos[2,i]])
    df <- data.frame("Column_1" = colnames(data)[combos[1,i]],
               "Column_1.Mean" = result$estimate[1],
               "Column_2" = colnames(data[combos[2,i]]),
               "Column_2.Mean" = result$estimate[2],
               "t_Value" = as.numeric(sprintf("%.3f", result$statistic)),
               "Deg_of_Freedom"= result$parameter,
               "p_value" = as.numeric(sprintf("%.3f", result$p.value)))
    res <- rbind(res,df)
  }
  res <- res[order(res[,"p_value"],decreasing = T),]
  row.names(res) <- 1:nrow(res)
  return(res)
}

#' @title func
#' @description  Returns the most repeating value in the Vector
#' @param vect
#' @return ind
#' @export mode

mode <- function(vect)
{
  `%>%` <- dplyr::`%>%`
  vect <- vect[!is.na(vect)]
  ind <- which(table(vect) == max(table(vect)))
  ind <- names(ind) %>% as.numeric()
  if(length(ind) > 1){
    if(sum(mean(vect) == ind) > 0)return(mean(vect))
    return(sort(ind)[round(length(ind)/2)])
  }else
  {
    return(ind)
  }
}

#' @title func
#' @description  Gives the Distribution Plots for all Numeric Columns in the DataFrame
#' @param df
#' @return plot
#' @export Numeric.Dist

Numeric.Dist <- function(df)
{
  `%>%` <- dplyr::`%>%`
  num_cols <- split.vectors(df)$Numeric_Columns
  if(length(num_cols)== 0)stop("Error: No Numeric Columns present in Dataset!")
  num_df <- df[,num_cols] %>% as.data.frame()
  names(num_df) <- num_cols
  numeric_df <- num_df %>% tidyr::gather(key,value) 
  mean.df <- data.frame("key" = unique(numeric_df$key),"Mean" = sapply(num_df,mean,na.rm = T),"Label" = rep("Mu",ncol(num_df)))
  median.df <- data.frame("key" = unique(numeric_df$key),"Median" = sapply(num_df,median,na.rm = T),"Label" = rep("Med",ncol(num_df)))
  mode.df <- data.frame("key" = unique(numeric_df$key),"Mode" = sapply(num_df,mode),"Label" = rep("Mod",ncol(num_df)))
  
 numeric_df %>% ggplot2::ggplot(ggplot2::aes(value)) + ggplot2::geom_density(ggplot2::aes(y = ..count..),fill = "steelblue4",size = 1) + 
   ggplot2::geom_vline(data = mean.df,ggplot2::aes(xintercept = Mean),linetype = "dashed",color = "red",size = 0.8) + 
   ggplot2::geom_vline(data = median.df,ggplot2::aes(xintercept = Median),linetype = "dashed",color = "green1",size = 0.8) +
   ggplot2::geom_vline(data = mode.df,ggplot2::aes(xintercept = Mode),linetype = "dashed",color = "yellow2",size = 0.8) + 
   ggplot2::facet_wrap(~key,scales = "free") + ggplot2::geom_text(data = mean.df,aes(x = Mean,label  = Label, y = 0.2),inherit.aes = F,color = "red")+
   ggplot2::geom_text(data = median.df,ggplot2::aes(x = Median,label  = Label, y = 0.4),inherit.aes = F,color = "green1") + 
   ggplot2::geom_text(data = mode.df,ggplot2::aes(x = Mode,label  = Label, y = 0.6),inherit.aes = F,color = "yellow2") + ggplot2::theme_bw()
}

#' @title func
#' @description  Gives the Distribution Plots for all Character Columns in the DataFrame
#' @param df
#' @return plot
#' @export Character.Dist

Character.Dist <- function(df)
{
  `%>%` <- dplyr::`%>%`
  char_cols <- split.vectors(df)$Character_Columns
  if(length(char_cols)== 0)stop("Error: No Character or Factor Columns present in Dataset!")
  char_df <- df[,char_cols] %>% as.data.frame()
  names(char_df) <- char_cols
  character_df <- char_df %>% tidyr::gather(key,value)
  character_df <- character_df %>% dplyr::group_by(key,value) %>% dplyr::summarise("Count" = n())  
  ggplot2::ggplot(character_df,ggplot2::aes(x = stats::reorder(value,Count), y = Count,fill = value)) + ggplot2::geom_bar(stat = "identity") + ggplot2::facet_wrap(~key,scales = "free")+ ggplot2::coord_flip() + ggplot2::guides(fill = F)
  
}

#' @title func
#' @description  Gives the pariwise annova test results for all numerical and categorical columns in the Dataset
#' @param df
#' @return res
#' @export DF.aov_ttest

DF.aov_ttest <- function(df)
{
  `%>%` <- dplyr::`%>%`
  df <- as.data.frame(df)
  n <- split.vectors(df)$Numeric_Columns
  if(length(n) == 0)stop("No Numeric Columns present in the Dataset")
  c <- split.vectors(df)$Character_Columns
  if(length(c) == 0)stop("No Character Columns present in the Dataset")
  num_df <- df[,n] %>% as.data.frame()
  names(num_df) <- n
  char_df <- df[,c] %>% as.data.frame()
  names(char_df) <- c
  char_cols <- c()
  num_cols <- c()
  p.value <- c()
  for(i in 1:ncol(char_df))
  {
    for(j in 1:ncol(num_df))
    {
      if(length(levels(as.factor(char_df[,i]))) <= 2)
      {
        char_cols <- c(char_cols,names(char_df)[i])
        num_cols <- c(num_cols,names(num_df)[j])
        p <- round(t.test(num_df[char_df[,i] == levels(as.factor(char_df[,i]))[1],j],num_df[char_df[,i] == levels(as.factor(char_df[,i]))[2],j])$p.value,4)
        p.value <- c(p.value,p)
      }else{
        char_cols <- c(char_cols,names(char_df)[i])
        num_cols <- c(num_cols,names(num_df)[j])
        r <- summary(aov(num_df[,j] ~ char_df[,i]))
        p <- ifelse(length(r[[1]]) < 5,0,round(r[[1]][5]$`Pr(>F)`[1],4))
        p.value <- c(p.value,p)
      }
    }
  }
  
  res <- data.frame("Character_Column" = char_cols, "Numeric_Column" = num_cols, "P_Value" = p.value)
  res <- res %>% dplyr::arrange(-P_Value)
  return(res)
}


