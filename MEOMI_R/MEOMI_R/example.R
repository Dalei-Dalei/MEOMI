  source("D:/manipulation function.R")
  source("D:/Conditional interaction information calculation.R")
  source("D:/MEOMI.R")
  ##Set the path to read the file and the path to output the file
  path_data <- "D:/InSilicoSize10-Ecoli1-heterozygous.tsv"
  path_standard <- "D:/DREAM3GoldStandard_InSilicoSize10_Ecoli1.txt"
  Path_out <- "D:/result.csv" 
  ##Read standard network and data
  mydata <- read.table(file = path_data,header = T, sep = "\t")
  net <- read.table(file = path_standard,header = F, sep = "\t")
  ##Delete the extra rows
  if(colnames(mydata)[1] != "G1")
  {
    mydata<-mydata[-1,-1]
  } 
  ##The standard network is processed to change direction to undirection
  for(m in 1: nrow(net)){      
    if(net[m,]$V3==1){
      p <- net[m,]$V1
      q <- net[m,]$V2
      if(net[(net$V1==q&net$V2==p),]$V3 == 0){
        net[(net$V1==q&net$V2==p),]$V3 <- 1
      }
    }
  }
  ##Set up parameters
  lisan <- 5
  lambda <- 0.1
  order <- 4
  Gval<-MEOMI(mydata = mydata, net = net, bins = lisan,lamda = lambda, order = order)
  r <- canshu(mydata, net, Gval[,c(1:3)])
  put_data(Path_out, r, data[i])
  View(r)
  
  
  
  
  
  
 
