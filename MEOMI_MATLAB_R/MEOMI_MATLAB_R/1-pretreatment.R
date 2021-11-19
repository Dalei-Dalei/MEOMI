source("D:/manipulation function.R")

##set up parameters
lisan <- 5

## Set the path to read and output files
path_data <- "D:/InSilicoSize10-Ecoli1-heterozygous.tsv"
path_out <- paste("D:/temp-", lisan,".csv", sep = "")

## Read standard network and data
mydata <- read.table(file = path_data,header = T, sep = "\t")
net <- read.table(file = path_standard,header = F, sep = "\t")
##Delete the extra rows
if(colnames(mydata)[1] != "G1")
{
  mydata<-mydata[-1,-1]
} 
temp<-get_grn0(mydata = mydata, bins = lisan)
write.table(temp, file = path_out, append = T, sep = ",", quote = F, row.names = F, col.names = F)
