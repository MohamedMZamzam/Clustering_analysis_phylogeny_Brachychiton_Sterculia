#load libraries#
library (tidyverse)
#Load the excel file with the binary numbers, 0 and 1, that resemble the status of each character state#
G <- readxl::read_xls("D://Al-Azhar University/Dr Ahmed Al Tahir/From Dr Fatima/New folder/Anatomy and Morphology Data Matrix 2003 .xls")
#Generate and tide the data frame#
G <- as.data.frame (G)
M <- G[-1,-1]
name <- G$...1 [2:10]
M <- `rownames<-`(M, name)
#Generate the distance matrix#
n <- dist (M)
#Make H clustering and plot the dendorogram#
clust.res<-hclust(n)


plot(clust.res)

