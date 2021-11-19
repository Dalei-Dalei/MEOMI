source("D:/manipulation function.R")
##set up parameters
lisan <- 5
lamda <- 0.1
## Set the path to read and output files
path_data <- "D:/InSilicoSize10-Ecoli1-heterozygous.tsv"
path_standard <- "D:/DREAM3GoldStandard_InSilicoSize10_Ecoli1.txt"
path_grn <- paste("D:/result_",lisan,"_",lamda,".txt",sep = "")
path_out <- "D:/result.csv" 
## Read standard network and data
mydata <- read.table(file = path_data,header = T, sep = "\t")
net <- read.table(file = path_standard,header = F, sep = "\t")
Gval <-read.table(file = path_grn, header = F, sep = "\t",)
##Delete the extra rows
if(colnames(mydata)[1] != "G1")
{
  mydata<-mydata[-1,-1]
} 
##Processes the standard network to change directed to undirected
for(m in 1: nrow(net)){      
  if(net[m,]$V3==1){
    p <- net[m,]$V1
    q <- net[m,]$V2
    if(net[(net$V1==q&net$V2==p),]$V3 == 0){
      net[(net$V1==q&net$V2==p),]$V3 <- 1
    }
  }
}
Gval[,3] <- Gval[,3]/max(Gval[,3])
r <- canshu(mydata, net, Gval[,c(1:3)])
put_data(path_out, r, data[i]) 
View(r)


  


