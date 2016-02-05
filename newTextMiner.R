newTextMiner <- function(PubId, pdbId, source, primaryCitation, oligomerNames, mutantProtein, fi = NULL ){
    
    lastResult = list()
    lastMutantResult = list()
    
    noOligomer = data.frame(matrix(NA,1,2))
    names(noOligomer) = c("pdbId", "PubId")
    
    noFullText = data.frame(matrix(NA,1,2))
    names(noFullText) = c("pdbId", "PubId")
    
    xml_full = 0
    no_oligomericState = 0
    xml_empty = 0
    
    for(i in 1:length(PubId)){
        
        #tryCatch({
            
            #ptm <- proc.time()
            
            result = as.data.frame(matrix(NA,1,6))
            names(result) = c("pdbId", "PubId", "oligomericState", "oligomericSentence","pubType", "primaryCitation")
            
            mutantResult = as.data.frame(matrix(NA,1,4))
            names(mutantResult) = c("pdbId", "PubId", "mutantKeyword", "mutantSentence")
            
            
            if(source[i]=="PDF"){
                
                sanitize <- function(str) {
                    gsub('([#$%&~_\\^\\\\{}\\s\\(\\)])', '\\\\\\1', str, perl = TRUE)
                }
                
                # get a list of all files in the current directory
                dir <- fi
                
                setwd(dir)
                
                files = list.files()
                
                for (i in 1:length(files)){
                    
                    
                    f2 <- sanitize(files[i])
                    system(paste0("pdftotext ", files[i]), wait = TRUE)
                    filetxt <- sub(".pdf", ".txt", f2)
                    
                    paper = sub("([.-])|[[:punct:]]", "\\1", filetxt)
                    
                    text <- readLines(paper, warn=FALSE)
                    
                    # adjust encoding of text - you have to know it
                    Encoding(text) <- "latin1"
                    
                    # Do something with the content - here: get word and character count of all pdfs in the current directory
                    article <- paste(text, collapse="\n") # collapse lines into one long string
                    
                    fullTextTotal = as.String(article)
                    numberOfChar = nchar(fullTextTotal)
                    
                }
                
                
            }
            
            if(source[i]=="PMC"){
                
                doc <- htmlParse(paste("http://www.ncbi.nlm.nih.gov/pmc/articles/",gsub("PMC", "", PubId[i]),"/", sep=""))
                
                sections = xpathSApply(doc, "//div[starts-with(@class, 'tsec sec')]", xmlValue)
                
                #if(length(sections)!=0){
                
                #fullTextTotal = as.String(unlist(sections))
                fullTextTotal = as.String(unlist(sections[-c(length(sections))]))
                
                numberOfChar = nchar(fullTextTotal)
                #}
                
            }
            if(source[i]=="PubMed"){
                
                
                xml.url <- paste("http://www.ebi.ac.uk/europepmc/webservices/rest/search/resulttype=core&query=ext_id:", PubId[i], sep="")
                
                xml_data <- xmlToList(xml.url)
                abstractXML = xml_data$resultList$result$abstractText
                if(length(abstractXML) != 0L){
                    title = xml_data$resultList$result$title
                    fullTextTotal <- as.String(c(paste("[title]", title), abstractXML))
                    numberOfChar = nchar(fullTextTotal)
                    
                }
                
                if(length(abstractXML) == 0L){
                    
                    html <- getURL(paste("http://www.ncbi.nlm.nih.gov/pubmed/?term=",PubId[i], sep=""), followlocation = TRUE)
                    doc = htmlParse(html, asText=TRUE)
                    fullTextTotal <- as.String(xpathSApply(doc, "//p", xmlValue))
                    numberOfChar = nchar(fullTextTotal)
                    
                    
                }
            }
            
            
            if(numberOfChar != 0){ #### sorun var
              
              xml_full =  xml_full+1
              
              sent_token_annotator <- Maxent_Sent_Token_Annotator()
              
              fullTextSentences <- annotate(fullTextTotal, sent_token_annotator)
              ## Extract sentences.
              fullTextSentencesList = fullTextTotal[fullTextSentences]
              
              z=1
              c=1
              for(l in 1:length(fullTextSentencesList)){
                
                review_source <- VectorSource(as.character(fullTextSentencesList)[l])
                corpus <- Corpus(review_source)
                corpus <- tm_map(corpus, stripWhitespace) # remove extra white spaeces
                corpus <- tm_map(corpus, content_transformer(tolower)) # convert lower cases
                corpus <- tm_map(corpus, removePunctuation) # remove punctuation
                corpus <- tm_map(corpus, removeWords, stopwords("english")) # remove stop words
                #corpus <- tm_map(corpus, removeNumbers) # remove numbers from the text
                
                dtm <- as.data.frame(as.matrix(DocumentTermMatrix(corpus)))
                
                abstractWords  = colnames(dtm)
                
                oligomericState = abstractWords[which(abstractWords %in% oligomerNames)]
                mutantProteinInfo = abstractWords[which(abstractWords %in% mutantProtein)]
                
                
                oligomericState = dtm[oligomericState]
                mutantProteinInfo = dtm[mutantProteinInfo]
                
                
                ### mutant
                
                
                if(length(mutantProteinInfo) != 0L) {
                  
                  # oligomericState = oligomericState
                  for(p in 1:length(mutantProteinInfo)){
                    
                    for(j in c:((c-1)+mutantProteinInfo[1,p])){
                      
                      mutantResult[j,1] = paste0("<a href='","http://www.rcsb.org/pdb/explore/explore.do?structureId=",as.character(pdbId[i])," 'target=","'_blank'","'>",as.character(pdbId[i]),"</a>")
                      mutantResult[j,2] = if(source[i] == "CrossRef"){paste0("<a href='","http://dx.doi.org/",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>")}else if(source[i] == "PMC"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pmc/articles/PMC",gsub("PMC","",PubId[i]),"/"," 'target=","'_blank'","'>",PubId[i],"</a>"))}else if(source[i] == "PubMed"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pubmed/?term=",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>"))}else{"PDF"}
                      mutantResult[j,3] = names(mutantProteinInfo[p])
                      mutantResult[j,4] = gsub("\n", " ", as.character(fullTextSentencesList)[l])
                      
                    }
                    c=c+mutantProteinInfo[1,p]
                    
                  }
                  
                  lastMutantResult[[i]] = mutantResult
                }
                
                
                ###oligomeric state
                
                if(length(oligomericState) != 0L) {
                  
                  # oligomericState = oligomericState
                  for(k in 1:length(oligomericState)){
                    
                    
                    for(j in z:((z-1)+oligomericState[1,k])){
                      
                      result[j,1] = paste0("<a href='","http://www.rcsb.org/pdb/explore/explore.do?structureId=",as.character(pdbId[i])," 'target=","'_blank'","'>",as.character(pdbId[i]),"</a>")
                      result[j,2] = if(source[i] == "CrossRef"){paste0("<a href='","http://dx.doi.org/",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>")}else if(source[i] == "PMC"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pmc/articles/PMC",gsub("PMC","",PubId[i]),"/"," 'target=","'_blank'","'>",PubId[i],"</a>"))}else if(source[i] == "PubMed"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pubmed/?term=",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>"))}else{"PDF"}
                      result[j,3] = names(oligomericState[k])
                      result[j,4] = gsub("\n", " ", as.character(fullTextSentencesList)[l])
                      result[j,5] = if(source[i]=="CrossRef"){("CrossRef")}else if(source[i]=="PMC"){("full text")}else if(source[i]=="PubMed"){("abstract")}else if(source[i]=="PDF"){("full text-PDF ")}
                      result[j,6] = primaryCitation[i]
                      
                      paste( sep="")
                      
                    }
                    z=z+oligomericState[1,k]
                    
                  }
                  
                  lastResult[[i]] = result
                }
                
              }
              
              if(is.na(result[,"oligomericState"])){
                
                result[1,1] = paste0("<a href='","http://www.rcsb.org/pdb/explore/explore.do?structureId=",as.character(pdbId[i])," 'target=","'_blank'","'>",as.character(pdbId[i]),"</a>")
                result[1,2] = if(source[i] == "CrossRef"){paste0("<a href='","http://dx.doi.org/",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>")}else if(source[i] == "PMC"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pmc/articles/PMC",gsub("PMC","",PubId[i]),"/"," 'target=","'_blank'","'>",PubId[i],"</a>"))}else if(source[i] == "PubMed"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pubmed/?term=",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>"))}else{"PDF"}
                result[1,3] = "NoResult"
                result[1,4] = "NoResult"
                result[1,5] = if(source[i]=="CrossRef"){("CrossRef")}else if(source[i]=="PMC"){("full text")}else if(source[i]=="PubMed"){("abstract")}else if(source[i]=="PDF"){("full text-PDF ")}
                result[1,6] = primaryCitation[i]
                
                
                
                
                lastResult[[i]] = result[complete.cases(result),]
              }
              
              
              if (is.na(result$oligomericState[1])) {
                
                no_oligomericState = no_oligomericState + 1
                noOligomer[i,1]= paste0("<a href='","http://www.rcsb.org/pdb/explore/explore.do?structureId=",as.character(pdbId[i])," 'target=","'_blank'","'>",as.character(pdbId[i]),"</a>")
                noOligomer[i,2]= if(source[i] == "CrossRef"){paste0("<a href='","http://dx.doi.org/",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>")}else if(source[i] == "PMC"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pmc/articles/PMC",gsub("PMC","",PubId[i]),"/"," 'target=","'_blank'","'>",PubId[i],"</a>"))}else if(source[i] == "PubMed"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pubmed/?term=",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>"))}else{"PDF"}
                
              }
            }
            else{
              
              result[j,1] = paste0("<a href='","http://www.rcsb.org/pdb/explore/explore.do?structureId=",as.character(pdbId[i])," 'target=","'_blank'","'>",as.character(pdbId[i]),"</a>")
              result[j,2] = if(source[i] == "CrossRef"){paste0("<a href='","http://dx.doi.org/",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>")}else if(source[i] == "PMC"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pmc/articles/PMC",gsub("PMC","",PubId[i]),"/"," 'target=","'_blank'","'>",PubId[i],"</a>"))}else if(source[i] == "PubMed"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pubmed/?term=",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>"))}else{"PDF"}
              result[j,3] = "NoResult"
              result[j,4] = "NoResult"
              result[j,5] = if(source[i]=="CrossRef"){("CrossRef")}else if(source[i]=="PMC"){("full text")}else if(source[i]=="PubMed"){("abstract")}else if(source[i]=="PDF"){("full text-PDF ")}
              result[j,6] = primaryCitation[i]
              
              paste( sep="")
              
              z=z+oligomericState[1,k]
              
              lastResult[[i]] = result[complete.cases(result),]
              
              xml_empty = xml_empty + 1
              noFullText[i,1] = paste0("<a href='","http://www.rcsb.org/pdb/explore/explore.do?structureId=",as.character(pdbId[i])," 'target=","'_blank'","'>",as.character(pdbId[i]),"</a>")
              noFullText[i,2]= if(source[i] == "CrossRef"){paste0("<a href='","http://dx.doi.org/",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>")}else if(source[i] == "PMC"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pmc/articles/PMC",gsub("PMC","",PubId[i]),"/"," 'target=","'_blank'","'>",PubId[i],"</a>"))}else if(source[i] == "PubMed"){(paste0("<a href='","http://www.ncbi.nlm.nih.gov/pubmed/?term=",PubId[i]," 'target=","'_blank'","'>",PubId[i],"</a>"))}else{"PDF"}
              
            }
            
            
            
       # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        
        
    }
    
    ### oligomeric state results
    lastResult2 = do.call(rbind.data.frame, lastResult)
    
    lastResult3 = unique(lastResult2[,-3])
    
    if(dim(lastResult2)[1]!=dim(lastResult3)[1]){
      
      
      for(i in 1:dim(lastResult3)[1]){
        
        keyword = lastResult2$oligomericState[which(lastResult3$oligomericSentence[i] == lastResult2$oligomericSentence)]
        
        if(length(keyword > 1)){
          
          keyword = paste(keyword, collapse = ", ")
          
        }
        
        
        lastResult3$oligomericState[i] =  keyword
        
      }
    }
    
    
    ##### mutant results
    lastMutantResult2 = do.call(rbind.data.frame, lastMutantResult)
    
    lastMutantResult3 = unique(lastMutantResult2[,-3])
    
    if(dim(lastMutantResult2)[1]!=dim(lastMutantResult3)[1]){
      
      
      for(i in 1:dim(lastMutantResult3)[1]){
        
        keyword = lastMutantResult2$mutantKeyword[which(lastMutantResult3$mutantSentence[i] == lastMutantResult2$mutantSentence)]
        
        if(length(keyword > 1)){
          
          keyword = paste(keyword, collapse = ", ")
          
        }else{
          
          keyword =  as.character(keyword)
        }
        
        
        lastMutantResult3$mutantKeyword[i] =  keyword
        
      }
    }else{
      
      lastMutantResult3 = lastMutantResult2
    }
    
    
    
    if(dim(lastResult2)[1]!=0){
      nofPubs = unique(lastResult2[,c(1,5)])
      nofPubs2 = data.frame(table(nofPubs[,2]))
      if(!is.na(nofPubs[2,2])){
        if(nofPubs2[1,1] == "full text"){
          nofPMC = nofPubs2[1,2]
          nofPubMed = 0
        }else{
          nofPubMed = nofPubs2[1,2]
          nofPMC = 0
        }
      }
      
      else{
        
        if(nofPubs[1,2] =="full text"){
          nofPMC= nofPubs2[1,2]
          nofPubMed = 0
        }
        else{
          nofPubMed=nofPubs2[1,2]
          nofPMC= 0
        }
      }
      
    }else{
      nofPubMed = 0
      nofPMC= 0
    }
    
    noOligomerResult = noOligomer[complete.cases(noOligomer),]
    
    noFullTextResult = noFullText[complete.cases(noFullText),]
    
    result = list(lastResult2, noOligomerResult, noFullTextResult, nofPMCArticles = nofPMC, nofPubMedArticles = nofPubMed, nofFoundedArticles=xml_full,
                  noOligomericState =  no_oligomericState, noFullText = xml_empty, numberOfuniquePdbId = length(unique(lastResult2$pdbId)), lastResult3, lastMutantResult3)
    
    return(result)
    
}
