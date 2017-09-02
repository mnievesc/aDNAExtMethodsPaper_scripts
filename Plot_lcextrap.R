###########################################################################################################
###### Script for analyzing DNA preservation dataset - aDNA methods paper
###### Written 6/12/2017 by M. Nieves Colon
###### Requires as input preseq results with library extrapolation estimates per each library per treatment.
###### IIV. Analysis 1.4: Extrapolation curves
##########################################################################################################


#################  1. Import data table and load functions  
##### Import table created in excel. Numeric columns must not have characters like n/a.
data <- read.table("../allsamples_forRanalysis.txt", header=TRUE, sep="\t")
head(data)

#### Import boxplot and normality assumption functions.
source("../AnalysisFunctions.R")


################# Only plotting for PI67 and GB7
################# Using set 1 colors to be consistent with previous graphs. 
################# D = red (1st color) H = blue (2nd color)
set1=brewer.pal(8,"Set1")


####################### 
#####  For PI-67  #####
####################### 

##################  Import and modify data for each sample treatment

##### 1. Import data table. Table created in excel. Numeric columns must not have 
#####  characters like n/a. Use for loop to import data for multiple samples

pi67_lc_d = read.table("preseq_downsampled_perind/PI-67-D.ds.lcextrap.txt", header=TRUE, sep="\t")
pi67_lc_d_long = melt(pi67_lc_d, id="TOTAL_READS")
pi67_lc_h = read.table("preseq_downsampled_perind/PI-67-H.ds.lcextrap.txt",  header=TRUE, sep="\t")
pi67_lc_h_long = melt(pi67_lc_h, id="TOTAL_READS")


pi67_plot = 
  ggplot(pi67_lc_d_long, aes(x=TOTAL_READS, y=subset(pi67_lc_d_long, variable=="EXPECTED_DISTINCT", select=value))) +
  geom_ribbon(data = pi67_lc_d_long, 
              aes(ymin=subset(pi67_lc_d_long, variable=="LOWER_0.95CI", select=value), 
                  ymax=subset(pi67_lc_d_long, variable=="UPPER_0.95CI", select=value)), 
              fill = "grey86", alpha=0.5) + 
  geom_line(data = pi67_lc_d_long, 
            aes(x=TOTAL_READS, y=subset(pi67_lc_d_long, variable=="EXPECTED_DISTINCT", select=value), color=set1[1]),
            size=1)  +
  geom_ribbon(data = pi67_lc_h_long,
              aes(ymin=subset(pi67_lc_h_long, variable=="LOWER_0.95CI", select=value), 
                  ymax=subset(pi67_lc_h_long, variable=="UPPER_0.95CI", select=value)), 
              fill = "grey86", alpha=0.5) +
  geom_line(data = pi67_lc_h_long, 
            aes(x=TOTAL_READS, y=subset(pi67_lc_h_long, variable=="EXPECTED_DISTINCT", select=value), color=set1[2]),
            size=1) +
  labs(y="Expected Distinct Reads", x="Total Sequenced Reads", title = "PI-67") +
  scale_color_brewer(palette="Set1", name="Extraction Method", labels=c("D", "H")) + 
  geom_vline(xintercept = 6876556, linetype = "dotted") +    ylim(0,1250000000) + # for same axes
  theme_bw()




####################### 
#####  For aGB7  #####
####################### 

##################  Import and modify data for each sample treatment

##### 1. Import data table. Table created in excel. Numeric columns must not have 
#####  characters like n/a. Use for loop to import data for multiple samples

gb7_lc_d = read.table("preseq_downsampled_perind/aGB7-D.ds.lcextrap.txt", header=TRUE, sep="\t")
gb7_lc_d_long = melt(gb7_lc_d, id="TOTAL_READS")

gb7_lc_h = read.table("preseq_downsampled_perind/aGB7-H.ds.lcextrap.txt", header=TRUE, sep="\t")
gb7_lc_h_long = melt(gb7_lc_h, id="TOTAL_READS")



#####################  Make line graphs with confidence intervals

gb7_plot = 
  ggplot(gb7_lc_d_long, aes(x=TOTAL_READS, y=subset(gb7_lc_d_long, variable=="EXPECTED_DISTINCT", select=value))) +
  geom_ribbon(data = gb7_lc_d_long, 
              aes(ymin=subset(gb7_lc_d_long, variable=="LOWER_0.95CI", select=value), 
                  ymax=subset(gb7_lc_d_long, variable=="UPPER_0.95CI", select=value)), 
              fill = "grey86", alpha=0.5) + 
  geom_line(data = gb7_lc_d_long, 
            aes(x=TOTAL_READS, y=subset(gb7_lc_d_long, variable=="EXPECTED_DISTINCT", select=value), color=set1[1]),
            size=1)  +
  geom_ribbon(data = gb7_lc_h_long,
              aes(ymin=subset(gb7_lc_h_long, variable=="LOWER_0.95CI", select=value),
                  ymax=subset(gb7_lc_h_long, variable=="UPPER_0.95CI", select=value)),
              fill = "grey86", alpha=0.5) +
  geom_line(data = gb7_lc_h_long,
            aes(x=TOTAL_READS, y=subset(gb7_lc_h_long, variable=="EXPECTED_DISTINCT", select=value), color=set1[2]),
            size=1) +
  labs(y="Expected Distinct Reads", x="Total Sequenced Reads", title = "GB-7") +
  scale_color_brewer(palette="Set1", name="Extraction Method", labels=c("D", "H")) + 
  geom_vline(xintercept = 5192848, linetype = "dotted") +    ylim(0,1250000000) + for same axes
  theme_bw()




################################### 
#####  Save two panel figure  #####
################################### 

pdf("ExtrapCurveMappedQ30.PI67-GB7.DownsampledPerInd.s100k.pdf")
grid.arrange(pi67_plot, gb7_plot)
dev.off()


pdf("ExtrapCurveMappedQ30.PI67-GB7.DownsampledPerInd.s100k.sameaxes.pdf")
grid.arrange(pi67_plot, gb7_plot)
dev.off()


ggsave("pi67.plot.pdf", pi67_plot, device = "pdf")
  
