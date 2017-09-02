##########################################################################################
###### Functions for paper performing comparative analysis of aDNA extraction methods
###### For sourcing into Rstudio
###### Written 5/19-27/2017 by M. Nieves-Colon
##########################################################################################

### Load necessary libraries
library(reshape2)     # for melting data
library(ggplot2)      # for plotting
library(gridExtra)    # for plotting multiple images with ggsave 
library(RColorBrewer) # for color palette
library(tidyr)        # for histogram overlay functions
library(car)          # for Levene's test 


### FUN 1. Get boxplot grouped by extraction method. Arguments are: 
### df = data frame, col = column of interest in string format, y lab text in string format
# Example: get_bp(df, "gc", Average GC content (%)")

get_bp = function(df, col, ylabtext) {
  
  # Subset dataframe by column of interest and melt
    sub = df[, c("Site", "Method", "Tissue", col)]
    df.m = melt(sub)
    print(head(df.m))
  
  # Create ggplot boxplot
    ggplot(df.m, aes(x=Method, y=value, fill=Method)) +
    geom_boxplot() + scale_fill_brewer(palette="Set1") + 
    scale_color_brewer(palette="Set1") +
    labs(x="Extraction Method", y=ylabtext) + 
    theme_bw()
}



### FUN 3. Get boxplot facet wrapped by Tissue. Arguments are: 
### df = data frame, col = column of interest, y lab text in string format
# Example: get_bp_facet(gc.m, "Average GC content (%)")

get_bp_facet = function(df, col, ylabtext) {

  # Subset dataframe by column of interest and melt
  sub = df[, c("Site", "Method", "Tissue", col)]
  df.m = melt(sub)
  print(head(df.m))
    
  # Create ggplot boxplot
  ggplot(df.m, aes(x=Method, y=value, fill=Method)) +
  geom_boxplot() + scale_fill_brewer(palette="Set1") + 
  facet_wrap( ~ Tissue) +
  scale_color_brewer(palette="Set1") +
  labs(x="Extraction Method", y=ylabtext) +
  theme_bw()
}



### FUN 4. Check normality assumptions with histogram, qqplot and sw test.
### Arguments are: x = variable column name of interest
### Example: check_norm(test, test$gc_perc, "GC content", "GC content")

check_norm = function(df, x, text, xlabtext) {

  diff= x[df$Method=="D"]-x[df$Method=="H"]    # Calculate difference of paired values
  #print(x[df$Method=="D"])
  
  # Set EPS device to print qqplot and histogram directly to EPS file. Set A4 size.
  # Will not print to screen.
  setEPS()
  postscript(paste0(text,"_normtest.eps"), width = 8.3, height = 11.7) # A4
  par(mfrow=c(2,1)) # For multipanel plot
  
  hist(diff,main = paste0("Histogram for difference in ",    # Generate histogram
                          text, " (method D - method H)"),
       cex.main = 0.90, col="gray", xlab=xlabtext)
  
  #Generate QQplot of difference between the two methods
  qq.diff = qqnorm(diff, xlab="Standardized normal values", ylab="Observed values",
                   main = paste0("Normal Quantile plot for difference in ",
                                 text, " (method D - method H)"), 
                   cex.main = 0.90) ; qqline(diff, col ="red")
  dev.off() # close EPS device

  shapiro.test(diff)   # Shapiro wilk test for normality print to screen
} 



### FUN 5. Run paired T tests for samples across method
### Arguments are: df = dataframe, dataframe$col, colname = name of variable of interest)
### Example: pairedt(test, test$gc_perc, "GC content")
### Function calculates difference of paired values for t-test, then extracts statistics results and prints
### them. It then makes a new dataframe that appends an existing dataframe in environment (called statdf)
### and appends the stat test results as a new row. This new dataframe is then reassigned with same name
### in global environment. This way the results of each stat test are recursively added to the dataframe.

pairedt = function(df, col, colname) {
  t = t.test(col[df$Method=="D"], col[df$Method=="H"], paired=TRUE)   # test takes difference of paired values
  row = c(colname, "paired_t_test", t$statistic, t$parameter, t$p.value)      
  #print(row)
  print(t)
  extdf <- rbind(statdf, row, stringsAsFactors=FALSE) # fixes "invalid factor level, NA generated" error
  assign("statdf", extdf, .GlobalEnv)           
}



### FUN 6. Run paired wilcoxon tests for samples across method. Use same logic as in above paired t function
### Arguments are: df = dataframe, dataframe$col, colname
### Example: wilcox(test, test$gc_perc, "GC content")

wilcox = function(df, col, colname) {
  w = wilcox.test(col[df$Method=="D"], col[df$Method=="H"], paired=TRUE, alternative = "two.sided")
  row = c(colname, "wilcox_Test", w$statistic, w$parameter, "null", w$p.value)      
  #print(row)
  print(w)
  extdf <- rbind(statdf, row, stringsAsFactors=FALSE) # fixes "invalid factor level, NA generated" error
  assign("statdf", extdf, .GlobalEnv)  
}



#### FUN 9. Per individual barplots 
## Arguments: df = data frame, colname = string column name of interest, mult = factor to 
## multiply data by in case its percent, ylabtext = string for y axis, .. for matching subset column
## Example: get_barplot(data, "perc_end", 100, "Endogenous content (%)", Tissue=="Petrous portion")

get_barplot = function(df, colname, mult, ylabtext,  ...) {
  
  # Subset dataframe by tissue. Must use eval to pass text for subset function
  ssubset <- deparse(substitute(...))
  x = subset(df, eval(parse(text = ssubset)))
  
  # Subset again so just column of interest is present, subset tissue column for sanity check
  sub = x[, c("Sample", "Method", "Tissue", colname)]
  #print(head(sub))

  # Melt dataset and multiply by desired factor (if none then enter 1)
  m = melt(sub)
  m$value <- m$value*mult
  print(head(m))
  
  bg = ggplot(data=m, aes(x=Sample, y=value, fill=Method)) +
    geom_bar(stat="identity", position=position_dodge()) + # no line...
    scale_fill_brewer(palette="Set1", name="Extraction\nMethod") + 
    labs(x="Sample", y=ylabtext) +
    theme_bw()
  
  return(bg)
  
}
  



#### FUN 10. Read length histograms
## Arguments: samplesvec = vector of sample basenames, sufx = filename suffix, 
## tissue = tooth/petrous portion, ylabtext = text for y axis as character.
## Example: get_readhist(samples, ".ds.uniqreads.", "teeth", "read bp")

get_readhist = function(samplesvec, sufx, tissue, ylabtext) {
  
  ##### 1. Import data table with read length counts. Numeric columns must not have 
  #####  characters like n/a. Use for loop to import data for multiple samples
  #####  Import and modify data for all samples at once.
  for (i in samplesvec) {
    print(paste0("Importing data for sample ", i))
    table = read.table(paste0("readlength_hist/",tissue,"/",i, sufx, "hist.txt"), header=F, sep="\t")
    table$sample <- rep(i, nrow(table))
    print(head(table))
    print(paste0("Number of rows in sample ", i, " table: ", nrow(table)))
    cat("\n")
    if(i==samplesvec[1]){mastertable<-table}else{mastertable<-rbind(mastertable, table)}
  }

  ### 2. Split sample name and extraction method into two columns. Since this also removes the sample 
  ### name delimiter fix and re-add to mastertable. 
  temptable = mastertable %>% separate(sample, c("samplename1", "samplename2", "extmethod"), "-")
  temptable$name <- paste(temptable$samplename1, temptable$samplename2, sep="-")
  
  ### 3. Replace original sample name with new name and add in ext method table
  mastertable$sample <- temptable$name
  mastertable$extmethod <-temptable$extmethod
  #print(head(mastertable))

  #### 4. Melt data into long format
  mastertable_long <- melt(mastertable, id=c("sample", "extmethod"))
  #print(head(mastertable_long))
  
  #### 5. Plot overlay histograms with facet wrap
  ggplot(mastertable_long, aes(x=value, fill=extmethod)) +
  facet_wrap( ~ sample, scales = "free") +
  geom_histogram(binwidth=3, alpha=.5, position="identity") +
  scale_fill_brewer(palette="Set1", name="Extraction\nMethod") +
  labs(x=ylabtext, y="Frequency") + theme_bw()
  
}


### References
## 1. http://stackoverflow.com/questions/5142842/export-a-graph-to-eps-file-with-r
## 2. http://stackoverflow.com/questions/14184258/putting-multiple-plots-in-a-a4-sheet-by-using-r-codes
## 3. http://christianlemp.com/blog/2014/02/05/How-I-Manage-Data-Projects-with-RStudio-and-Git.html
## 4. https://stackoverflow.com/questions/28370249/correct-way-to-specifiy-optional-arguments-in-r-functions
## 5. https://stackoverflow.com/questions/7310186/function-in-r-passing-a-dataframe-and-a-column-name
## 6. https://stackoverflow.com/questions/10689055/create-an-empty-data-frame
## 7. https://stackoverflow.com/questions/16819956/invalid-factor-level-na-generated
## 8. https://stackoverflow.com/questions/10904124/global-and-local-variables-in-r
## 9. https://stat.ethz.ch/pipermail/r-help/2004-February/046366.html
## 11. https://stackoverflow.com/questions/11880906/pass-subset-argument-through-a-function-to-subset
## 12. https://cran.r-project.org/web/packages/muStat/muStat.pdf
## 13. http://www.dummies.com/programming/r/how-to-split-strings-in-r/
## 14. http://astrostatistics.psu.edu/su07/R/html/base/html/list.files.html
## 15. https://stackoverflow.com/questions/3969852/update-data-frame-via-function-doesnt-work
## 16. https://stackoverflow.com/questions/18028225/r-list-files-with-multiple-conditions
## 17. https://stackoverflow.com/questions/35618260/removing-legend-ggplot-2-2
## 18. http://ggplot2.tidyverse.org/reference/geom_point.html
## 19. https://stackoverflow.com/questions/11846295/how-to-add-different-lines-for-facets
## 20. https://stackoverflow.com/questions/12188509/cleaning-inf-values-from-an-r-dataframe
