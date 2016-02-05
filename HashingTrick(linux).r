rm(list=ls())
#setwd("/Users/selcukorkmaz/Dropbox/GSD/Studies/Web-Tools(Devel)/PDBTextMiner/")
library("caret")
library("XML")
library("NLP")
library("tm")
library("FeatureHashing")
library("openNLP")
library("numbers")
library("Rmpfr")
set.seed(1234)

positives=read.table("positives.txt",T, sep="\t", quote="")
dim(positives)
head(positives)
negatives=read.table("negatives2.txt", T, sep="\t", quote="")
dim(negatives)
head(negatives)

pos_sample = positives[sample(nrow(positives), 4000),]
neg_sample = negatives[sample(nrow(negatives), 4000),]

data = rbind(pos_sample, neg_sample)
dim(data)
head(data)

sentencesList = as.data.frame(data[1])

wordList = data.frame(matrix(NA,1,2))
names(wordList) = c("words", "label")

for(i in 1:dim(data)[1]){
    review_source <- VectorSource(sentencesList[i,1])
    corpus <- Corpus(review_source)
    corpus <- tm_map(corpus, stripWhitespace) # remove extra white spaces
    corpus <- tm_map(corpus, content_transformer(tolower)) # convert lower cases
    corpus <- tm_map(corpus, removePunctuation) # remove punctuation
    corpus <- tm_map(corpus, removeWords, stopwords("english")) # remove stop words
    #corpus <- tm_map(corpus, removeNumbers) # remove numbers from the text
    #corpus <- tm_map(corpus, stemDocument)
    dtm <- as.data.frame(as.matrix(DocumentTermMatrix(corpus)))
    wordList[i,1] = gsub("[\r\n]", " ", as.String(colnames(dtm)))
    wordList[i,2] = data[i,2]
}
#write.table(wordList, "/Users/selcukorkmaz/Dropbox/GSD/Studies/Web-Tools(Devel)/PDBTextMiner/wordLisr.txt",quote=F, row.names=F, sep="\t")

head(wordList)
dim(wordList)


docs <- VCorpus(VectorSource(wordList[,1]))
dtm <- DocumentTermMatrix(docs)
dim(dtm)


#### Machine learning
set.seed(1234)
splitIndex = createDataPartition(wordList$label, p=0.80, list=F, times = 1) # split data set as 70% training 30% testing data set
trainSplit = wordList[splitIndex,]
dim(trainSplit)
testSplit = wordList[-splitIndex,]
dim(testSplit)

trainData = as.data.frame(trainSplit[,-dim(trainSplit)[2]])
colnames(trainData)="words"
trainClass = as.factor(ifelse(trainSplit[,dim(trainSplit)[2]]==1,"yes","no"))

testData = as.data.frame(testSplit[,-dim(testSplit)[2]])
colnames(testData)="words"
testClass = as.factor(ifelse(testSplit[,dim(testSplit)[2]]==1, "yes","no"))

hashSize = round(nextPrime(dim(dtm)[2])/0.75,0)
save(hashSize, file="hashSize.rda")

train_hashed <- hashed.model.matrix(~ split(words, delim = " ", type = "tf-idf") -1,
                                        data = trainData, hash.size = hashSize, signed.hash = FALSE)

as.integer(which(train_hashed[40, ] != 0))


test_hashed <- hashed.model.matrix(~ split(words, delim = " ", type = "tf-idf") -1,
                                        data = testData, hash.size = hashSize , signed.hash = FALSE)

testSet = list(test_hashed = test_hashed, testClass = testClass)
save(testSet, file="testSet.rda")

as.integer(which(test_hashed[4, ] != 0))


## Parallel Programming Background
#nCores = 2
#library("doMC")
#library("foreach")
#library("digest")

#registerDoMC(nCores)

ctrl = trainControl(method = "repeatedcv", number=10, classProbs = T, summaryFunction=twoClassSummary)

#### 1. boosted logistic regression
lr = train(train_hashed, trainClass, method = "LogitBoost", tuneLength = 10, trControl = ctrl, metric="ROC")
lr

save(lr, file = "blr.rda")


lrpred = predict.train(lr, test_hashed)
lrpred2 = predict.train(lr, as.matrix(test_hashed), type="prob")
lrRes = confusionMatrix(table(lrpred, testClass),positive = "yes")
f1s_lrRes=2*lrRes$byClass[1]*lrRes$byClass[3]/(lrRes$byClass[1]+lrRes$byClass[3])
names(f1s_lrRes)="F1S"
mcc_lr1= mpfr((lrRes$table[1,1]*lrRes$table[2,2]-lrRes$table[1,2]*lrRes$table[2,1]), precBits = 5)
mcc_lr2 = mpfr(((lrRes$table[1,1]+lrRes$table[1,2])), precBits = 5)
mcc_lr3 = mpfr(((lrRes$table[1,1]+lrRes$table[2,1])*(lrRes$table[2,2]+lrRes$table[1,2])*(lrRes$table[2,2]+lrRes$table[2,1])), precBits = 5)            
mcc_lrRes =  as.numeric(mcc_lr1/sqrt(mcc_lr2*mcc_lr3))  
names(mcc_lrRes) = "MCC"
lrAUC = as.numeric(roc(testClass, lrpred2[,"yes"], levels = c("yes","no"), plot=F)$auc)*1
names(lrAUC) = "AUC"
lrRes =  c(lrRes$overall[1:2],lrAUC,lrRes$byClass[c(-5,-7)], f1s_lrRes, mcc_lrRes)
lrRes

#### 2. support vector machines with radial kernel
svm = train(train_hashed, trainClass, method = "svmRadial", tuneLength = 10, trControl = ctrl, metric="ROC")
svm

svmpred = predict.train(svm, test_hashed)
svmRes = confusionMatrix(table(svmpred, testClass),positive = "yes")
svmpred2 = predict.train(svm, as.matrix(test_hashed), type="prob")
f1s_svmRes=2*svmRes$byClass[1]*svmRes$byClass[3]/(svmRes$byClass[1]+svmRes$byClass[3])
names(f1s_svmRes)="F1S"
mcc_svm1= mpfr((svmRes$table[1,1]*svmRes$table[2,2]-svmRes$table[1,2]*svmRes$table[2,1]), precBits = 5)
mcc_svm2 = mpfr(((svmRes$table[1,1]+svmRes$table[1,2])), precBits = 5)
mcc_svm3 = mpfr(((svmRes$table[1,1]+svmRes$table[2,1])*(svmRes$table[2,2]+svmRes$table[1,2])*(svmRes$table[2,2]+svmRes$table[2,1])), precBits = 5)            
mcc_svmRes =  as.numeric(mcc_svm1/sqrt(mcc_svm2*mcc_svm3))            
names(mcc_svmRes) = "MCC"
svmAUC = as.numeric(roc(testClass, svmpred2[,"yes"], levels = c("yes","no"), plot=F)$auc)*1
names(svmAUC) = "AUC"
svmRes =  c(svmRes$overall[1:2],svmAUC,svmRes$byClass[c(-5,-7)], f1s_svmRes, mcc_svmRes)
svmRes



TestSetPerformance = t(rbind(svmRes, lrRes))
colnames(TestSetPerformance) = c("SVM", "BLR")

save(TestSetPerformance, file = "TestSetPerformance.rda")


svmROCdata =  as.data.frame(svmpred2[,2])
names(svmROCdata) = "SVM"
lrROCdata =  as.data.frame(lrpred2[,2])
names(lrROCdata) = "BLR"

class =as.data.frame(testClass)

ROCdata = cbind(svmROCdata,lrROCdata,class)
write.table(ROCdata, "~/Dropbox/GSD/Studies/Web-Tools(Devel)/PDBTextMiner(Doküman)/Diğer/ROCdata.txt",quote=F,row.names=F,sep="\t" )

###########################################################################################################
### 3. random forest
rf = train(as.matrix(train_hashed), trainClass, method = "rf", tuneLength = 10, trControl = ctrl,
           metric="ROC")
rf

rfpred = predict.train(rf, test_hashed)
rfRes = confusionMatrix(table(rfpred, testClass),positive = "yes")
rfpred2 = predict.train(rf, as.matrix(test_hashed), type="prob")
f1s_rfRes=2*rfRes$byClass[1]*rfRes$byClass[3]/(rfRes$byClass[1]+rfRes$byClass[3])
names(f1s_rfRes)="F1S"
mcc_rfRes=(rfRes$table[1,1]*rfRes$table[2,2]-rfRes$table[1,2]*rfRes$table[2,1])/sqrt((rfRes$table[1,1]+rfRes$table[1,2])*(rfRes$table[1,1]+rfRes$table[2,1])*(rfRes$table[2,2]+rfRes$table[1,2])*(rfRes$table[2,2]+rfRes$table[2,1]))
names(mcc_rfRes) = "MCC"
rfAUC = as.numeric(roc(testClass, rfpred2[,"yes"], levels = c("yes","no"), plot=F)$auc)*1
names(rfAUC) = "AUC"
rfRes =  c(rfRes$overall[1:2],rfAUC,rfRes$byClass[c(-5,-7)], f1s_rfRes, mcc_rfRes)
rfRes

###################################################################################################

perf = t(rbind(svmRes,lrRes, rfRes))
colnames(perf) = c("SVM", "BLR", "RF")
perf

save(perf, file = "performance.rda")

svmROC = plot.roc(testClass, svmpred2[,"yes"], col="red")
lrROC = lines.roc(testClass,lrpred2[,"yes"], col="blue")
rfROC = lines.roc(testClass, rfpred2[,"yes"], col="#008600")
legend("bottomright", legend=c("SVM", "LR", "RF"), col=c("red", "blue", "#008600"), lwd=2)

