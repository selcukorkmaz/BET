rm(list=ls())
setwd("~/Dropbox/GSD/Studies/Web-Tools(Devel)/BAEP_server_Devel/BET/")

### packages (start) #####
library("caret")
library("stringr")
library("XML")
library("tm")
#library("SnowballC")
library("NLP")
library("RCurl")
library("openNLP")
library("shinyIncubator")
library("shinyBS")
library("pROC")
library("FeatureHashing")
library("dplyr")
library("DT")
source("newTextMiner.R")
load("blr.rda")
load("svm.rda")
load("hashSize.rda")
load("TestSetPerformance.rda")
load("oligomericStateList.rda")
load("consistencyList.rda")
### packages (end) #####


data = read.table("Pisa_results.txt", header=T, sep = "\t", stringsAsFactors = F) 
dim(data)  
tail(data)

#data[data$pdbId == "7XIM",]

#data = data[1,]

################## CURRENT PDB RESULTS (start)###############################
listSplit = list()

for(p in 1:dim(data)[1]){

index95 = which(consistencyList$consistencyData95$"PDB ID" %in% as.character(data$pdbId[p]))
if(length(index95) > 0){
result = consistencyList$consistencyData95[index95, c("PDB ID",	"BA Number", "Stoichiometry", "Symmetry","Representative", "Consistency score", "Result")]

#currentResult = result[1:4]

pdbList = result
sequenceCluster = consistencyList$consistencyData95[consistencyList$consistencyData95[,3] %in% pdbList$Representative,]

cols = c("Stoichiometry", "Symmetry")

sequenceCluster$combinedOligomericState <- apply(sequenceCluster[,cols ] , 1 , paste , collapse = "_" )



sequenceCluster2 = unique(sequenceCluster[,c("PDB ID","Representative","Stoichiometry", "Symmetry", "Consistency score", "combinedOligomericState")])




splitData = split(sequenceCluster2,factor(sequenceCluster2$Representative))
signatureDataList = list ()

for(i in 1:length(splitData)){
  
  signatureData = splitData[[i]]
  
  for(j in 1:dim(signatureData)[1]){
    
    signatureData$count[j] = sum(signatureData$combinedOligomericState == signatureData$combinedOligomericState[j])
    
  }
  
  tbl = table(factor(signatureData$`PDB ID`), factor(signatureData$combinedOligomericState))
  
  sequenceCluster3 = unique(signatureData)
  sequenceCluster3$pdbIds = NA
  
  if(length(unique(signatureData$`PDB ID`)) ==1){
    
    sequenceCluster3$pdbIds = signatureData$`PDB ID`
    
  }else{
    
    for(j in 1:dim(tbl)[2]){
      
      
      sequenceCluster3[sequenceCluster3$combinedOligomericState == colnames(tbl)[j],][,"pdbIds"] = paste(names(which(tbl[,j] == 1)), collapse=", ")
      
      
    }
    
  }
  sequenceCluster3
  
  names(sequenceCluster3)[c(2,7:8)] = c("Representative", "# of PDBs", "PDB IDs")
  
  sequenceCluster4 = sequenceCluster3[-6]
  
  
  
  signatureDataList[[i]] = sequenceCluster4
}

lastResult = do.call(rbind.data.frame, signatureDataList)

lastResult2 = unique(lastResult[,-1])

sequenceCluster = lastResult2

sequenceCluster$PDBID = data$pdbId[p]

sequenceClusterLast = sequenceCluster[,c(7,1:6)]
names(sequenceClusterLast) = c("PDB ID", "Representative chain","Stoichiometry", "Symmetry",
                               "Consistency score", "Number of PDB entries", "PDB entries in the cluster")


maxCS = sequenceClusterLast$`Consistency score`[which.max(sequenceClusterLast$`Consistency score`)]

sequenceClusterLast$CS = maxCS

if(maxCS > 0.5 && sum(sequenceClusterLast[,6]) > 3){
  
  for(i in 1: dim(sequenceClusterLast)[1]){
    sequenceClusterLast$Res[i] = if(sequenceClusterLast$CS[i] == sequenceClusterLast$`Consistency score`[i]){1}else{0}
  }
} else{
  
  sequenceClusterLast$Res = 0
}

listSplit[[p]] = sequenceClusterLast[-7]

}
print(paste0(p,"/",dim(data)[1]))
}

sequenceClusterLastRes = do.call(rbind.data.frame, listSplit)

names(sequenceClusterLastRes) = c("PDB ID",	"Representative chain",	"Stoichiometry",	"Symmetry",	"Consistency score",	"Number of PDB entries",	"CS",	"Res")
write.table(sequenceClusterLastRes, "sequenceClusterLastRes_95.txt", quote=F, sep="\t",row.names = F)


