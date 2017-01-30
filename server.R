shinyServer(function(input, output, session){
  library("tm")
  library("caret")
  library("DT")
  library("plyr")
  #library("dplyr")
  
  library("data.table")
  library("stringr")
  library("XML")
  #library("NLP")
  #library("SnowballC")
  library("RCurl")
  #detach("package:bitops", unload=TRUE, force =TRUE)
  library("openNLP")
  library("shinyIncubator")
  library("shinyBS")
  #library("pROC")
  library("FeatureHashing")
  library("httr")
  source("newTextMiner.R")
  load("blr.rda")
  load("svm.rda")
  load("hashSize.rda")
  load("TestSetPerformance.rda")
  load("oligomericStateList.rda")
  load("consistencyList.rda")
  load("consistencyScore.rda")
  load("eppic_results.rda")
  load("publication.rda")
  load("pisa_results.rda")
  load("pisa_advanced_results.rda")
  
  
  
  
  ###############################
  ###### Data - Start ###########
  ###############################
  currentData <- reactive({
    
    if(input$startAnalysis == 0)
    {
      return()
    }
    
    isolate({
      if(input$startAnalysis){
        
        withProgress(message = 'Finding publications...',  detail = 'Please wait...', value=0,{
          
          
          PdbId <- reactive({
            pdb = as.character(input$PdbId)
            toupper(pdb)
          })
          
          #dataMentions = read.table("PdbDataMentionUnique.tsv", header=T, fill=T, sep="\t")
          pubData <- publication
          pubData = pubData[!duplicated(pubData$pdbId),]
          rownames(pubData) = pubData$pdbId
          
          #PMCId <- reactive({
          #    as.numeric(input$PMCId)
          #})
          
          data_pmcSingle <- reactive({
            
            #if(!input$dataMentions){
            
            raw = pubData[PdbId(),]
            otherPDB = pubData[pubData[,2] == raw$pubId,]
            
            if(length(otherPDB$pdbId) > 1){
              
              otherPDBResult = otherPDB[otherPDB[,1]!=raw$pdbId,]
              
              otherPDBList = paste(otherPDBResult$pdbId, collapse=", ")
              
              
              raw$otherPDBList = otherPDBList
              
            }else{
              
              raw$otherPDBList = "none"
            }
            
            colnames(raw) = c("PDB ID", "Publication ID", "Source", "Primary Citation", "Related Structures")
            
            #}
            
            #Sys.sleep(2)
            #if(input$dataMentions){
            #  primaryCitation = pubData[PdbId(),]
            #  currentDataention = dataMentions[dataMentions$pdbId==PdbId(),]
            #  raw = rbind(primaryCitation, currentDataention)
            #  raw = unique(raw)
            #
            #  pdbRelatedPub = raw[raw$primary_citation == "yes",]
            #
            #  otherPDB = pubData[pubData$pubId %in% pdbRelatedPub$pubId,]
            #
            #
            #  if(length(otherPDB$pdbId) > 1){
            #
            #    otherPDBResult = otherPDB[otherPDB[,1]!=PdbId(),]
            #
            #    otherPDBList = paste(otherPDBResult$pdbId, collapse=", ")
            #
            #    for(i in 1:length(raw$pdbId)){
            #
            #      if(raw$primary_citation[i] == "yes"){
            #
            #        raw$otherPDBList[i] = otherPDBList
            #      }else{
            #
            #        raw$otherPDBList[i] = "NA"
            #
            #      }
            #
            #    }
            #
            #  }else{
            #
            #    raw$otherPDBList = "none"
            #  }
            #
            #  colnames(raw) = c("PDB ID", "Publication ID", "Source", "Primary Citation", "Related Structures")
            #
            #
            #
            #}
            if(!is.null(return(raw))){
              return(raw)
            }
            
            else if(!is.null(return(raw2))){
              stop("Not available")
            }
            
          })
          
          
          if (input$dataInput==1)  {
            data = data.frame(data_pmcSingle())
          }
          
          else if(input$dataInput==2){  ## Upload data.
            inFile <- input$upload
            
            if (is.null(input$upload))  {return(NULL)}
            
            if (file.info(inFile$datapath)$size <= 10485800){
              data <- read.table(inFile$datapath, header=TRUE, fill=TRUE)
              data= as.data.frame(factor(toupper(data[,1])))
              pubData <- read.table("publication.txt", header=TRUE, sep="\t")
              rownames(pubData) = pubData$pdbId
              primaryCitation = pubData[as.character(data[,1]),]
              otherPDB = pubData[pubData[,2] %in% primaryCitation$pubId,]
              
              
              #if(!input$dataMentions){
              
              if(length(otherPDB$pdbId) > 1){
                
                otherPDBResult = otherPDB[!(otherPDB[,1] %in% primaryCitation$pdbId),]
                
                if(length(otherPDBResult$pdbId)!=0){
                  
                  splitOther = split(otherPDBResult, factor(otherPDBResult$pubId))
                  
                  otherPDBList = list()
                  
                  for(j in 1:length(splitOther)){
                    
                    s = splitOther[[j]]
                    
                    otherPDBsplit = paste(s$pdbId, collapse=", ")
                    
                    s$otherPDBsplit = otherPDBsplit
                    
                    otherPDBList[[j]] = s[1,]
                    
                  }
                  splitOtherPDBResult = as.data.frame(rbindlist( otherPDBList))
                  
                }else{
                  
                  splitOtherPDBResult = otherPDB
                  
                }
                
              }else{
                
                otherPDBResult = otherPDB[!(otherPDB[,1] %in% primaryCitation$pdbId),]
                splitOtherPDBResult = otherPDB
                
                
              }
              
              if(length(otherPDBResult$pdbId)!=0){
                for(i in 1:dim(primaryCitation)[1]){
                  
                  if(primaryCitation[i,"pubId"] %in% splitOtherPDBResult[,"pubId"] ){
                    
                    pubid = as.character(splitOtherPDBResult$pubId[which(primaryCitation[i,"pubId"] == splitOtherPDBResult[,"pubId"])])
                    
                    primaryCitation$otherPDBList[primaryCitation$pubId == pubid] =   splitOtherPDBResult[which(splitOtherPDBResult$pubId == pubid),"otherPDBsplit"]
                    
                    
                  }else{
                    
                    primaryCitation$otherPDBList[i] = "none"
                    
                  }
                  
                }}else{
                  
                  primaryCitation$otherPDBList = "none"
                  
                }
              
              
              data = primaryCitation
              colnames(data) = c("PDB ID", "Publication ID", "Source", "Primary Citation", "Related Structures")
              
              
              #}
              
              #if(input$dataMentions){
              #
              #if(length(otherPDB$pdbId) > 1){
              #
              #  otherPDBResult = otherPDB[!(otherPDB[,1] %in% primaryCitation$pdbId),]
              #  splitOther = split(otherPDBResult, factor(otherPDBResult$pubId))
              #
              #  otherPDBList = list()
              #
              #  for(j in 1:length(splitOther)){
              #
              #    s = splitOther[[j]]
              #
              #    otherPDBsplit = paste(s$pdbId, collapse=", ")
              #
              #    s$otherPDBsplit = otherPDBsplit
              #
              #    otherPDBList[[j]] = s[1,]
              #
              #  }
              #
              # }
              #
              # splitOtherPDBResult = do.call(rbind.data.frame, otherPDBList)
              
              #for(i in 1:dim(primaryCitation)[1]){
              
              # if(primaryCitation[i,"pubId"] %in% splitOtherPDBResult[,"pubId"] ){
              
              
              #primaryCitation$otherPDBList[factor(primaryCitation$pubId) == factor(primaryCitation$pubId[which(primaryCitation[i,"pubId"] == splitOtherPDBResult[,"pubId"])])]
              
              #    pubid = as.character(splitOtherPDBResult$pubId[which(primaryCitation[i,"pubId"] == splitOtherPDBResult[,"pubId"])])
              
              #    primaryCitation$otherPDBList[primaryCitation$pubId == pubid] =   splitOtherPDBResult[which(splitOtherPDBResult$pubId == pubid),"otherPDBsplit"]
              
              
              #   }else{
              
              #     primaryCitation$otherPDBList[i] = "none"
              
              #   }
              
              # }
              
              
              # dataMentions = read.table("PdbDataMentionUnique.tsv", header=T, fill=T, sep="\t")
              # currentDataention = dataMentions[dataMentions$pdbId%in%as.character(data[,1]),]
              # currentDataention$otherPDBList = "none"
              # data = rbind(primaryCitation, currentDataention)
              
              # colnames(data) = c("PDB ID", "Publication ID", "Source", "Primary Citation", "Related Structures")
              
              
              #}
            }
          }
          
          else if (input$dataInput==3){  ## Upload PDF
            
            inPDF <- input$uploadPDF
            
            if (is.null(input$uploadPDF))  {return(NULL)}
            
            else{
              
              #fi <- read.table(inPDF$datapath, header=FALSE, sep="\t", stringsAsFactors = F)
              
              #pdffile = input$uploadPDF
              
              PubId = rep("PDF", 1)
              pdbId = toupper(input$PdbIdPDF)
              source = rep("PDF", 1)
              primaryCitation = rep("yes", 1)
              relatedStructures = rep("none", 1)
              
              #sourcePDF = data.frame(matrix(NA,1,5))
              #sourcePDF[,1] = pdbId
              #sourcePDF[,2] = PubId
              #sourcePDF[,3] = source
              #sourcePDF[,4] = primaryCitation
              #sourcePDF[,5] = relatedStructures
              sourcePDF = data.frame(pdbId,PubId,source,primaryCitation,relatedStructures)
              
              colnames(sourcePDF) = c("PDB ID", "Pub ID", "Source", "Primary citation", "Related Structures")
              
              data = sourcePDF
              
            }
            
            
          }
          
          
          return(data)
          
        }
        )}
    })
    
  })
  
  
  dataPath <- reactive({
    
    input$uploadPDF
    #data2 = read.table(inPDF$datapath, header=FALSE, fill=TRUE, stringsAsFactors = F)
    #assign('data',data2,envir=.GlobalEnv)
    
    
    
    
  })
  
  
  
  
  
  
  ###############################
  ###### Data - End #############
  ###############################
  
  ##### Observe functions #########
  
  #observe({
  
  #    if(input$AdvancedOligomericStatePrediction){
  #        updateCheckboxInput(session, "AdvancedPisaPrediction", label = "PISA", FALSE)
  # updateCheckboxInput(session, "AdvancedEppic", label = "EPPIC", FALSE)
  #updateCheckboxInput(session, "AdvancedOutlierDetection", label = "Sequence cluster", FALSE)
  
  #    }
  #})
  
  
  #observe({
  
  #    if(input$AdvancedPisaPrediction){
  #        updateCheckboxInput(session, "AdvancedOligomericStatePrediction", label = "Text mining", FALSE)
  # updateCheckboxInput(session, "AdvancedEppic", label = "EPPIC", FALSE)
  # updateCheckboxInput(session, "AdvancedOutlierDetection", label = "Sequence cluster", FALSE)
  
  
  #    }
  
  
  #})
  
  
  #observe({
  
  #    if(input$AdvancedEppic){
  #        updateCheckboxInput(session, "AdvancedOligomericStatePrediction", label = "Text mining", FALSE)
  #        updateCheckboxInput(session, "AdvancedPisaPrediction", label = "PISA", FALSE)
  #        updateCheckboxInput(session, "AdvancedOutlierDetection", label = "Sequence cluster", FALSE)
  
  
  #    }
  
  
  #})
  
  
  #observe({
  
  #    if(input$AdvancedOutlierDetection){
  #        updateCheckboxInput(session, "AdvancedOligomericStatePrediction", label = "Text mining", FALSE)
  #        updateCheckboxInput(session, "AdvancedPisaPrediction", label = "PISA", FALSE)
  #updateCheckboxInput(session, "AdvancedEppic", label = "EPPIC", FALSE)
  
  
  #    }
  
  
  #})
  
  
  observe({
    
    if(!input$advancedResults){
      updateCheckboxInput(session, "AdvancedOligomericStatePrediction", label = "Text mining", FALSE)
      updateCheckboxInput(session, "AdvancedPisaPrediction", label = "PISA", FALSE)
      #updateCheckboxInput(session, "AdvancedEppic", label = "EPPIC", FALSE)
      #updateCheckboxInput(session, "AdvancedOutlierDetection", label = "Sequence cluster", FALSE)
      
      
      
    }
    
    
  })
  
  
  
  observe({
    
    if(input$advancedResults){
      updateCheckboxInput(session, "advancedOptions", label = "Options", FALSE)
    }
    
    
  })
  
  
  observe({
    
    if(input$advancedOptions){
      updateCheckboxInput(session, "advancedResults", label = "Results", FALSE)
    }
    
    
  })
  
  #observe({
  #    if(input$selectAll){
  #
  #  updateCheckboxGroupInput(session, inputId = "OligomericStatePrediction",      selected = "OligomericStatePrediction")
  #  updateCheckboxGroupInput(session, inputId = "oliomericPredPlot",      selected = "oliomericPredPlot")
  #  updateCheckboxGroupInput(session, inputId = "outlierDetection",      selected = "outlierDetection")
  #  updateCheckboxGroupInput(session, inputId = "outlierDetectionPlot",      selected = "outlierDetectionPlot")
  #  updateCheckboxGroupInput(session, inputId = "signatureResults",      selected = "signatureResults")
  #  updateCheckboxGroupInput(session, inputId = "mlFilteredSen",      selected = "mlFilteredSen")
  #  updateCheckboxGroupInput(session, inputId = "mlFiltering",      selected = "mlFiltering")
  #  updateCheckboxGroupInput(session, inputId = "noOligomericResults",      selected = "noOligomericResults")
  #  updateCheckboxGroupInput(session, inputId = "noFullTextResults",      selected = "noFullTextResults")
  #  updateCheckboxGroupInput(session, inputId = "publicationInfo",      selected = "publicationInfo")
  #  updateCheckboxGroupInput(session, inputId = "summaryResults",      selected = "summaryResults")
  #  updateCheckboxGroupInput(session, inputId = "mutantResults",      selected = "mutantResults")
  #  updateCheckboxInput(session, "deselectAll", label = "Deselect All", FALSE)
  #  updateCheckboxInput(session, "selectAll",   label = "Select All",   TRUE)
  #
  #    }
  #})
  
  
  
  #observe({
  #    if(input$deselectAll){
  #
  #        updateCheckboxGroupInput(session, inputId = "OligomericStatePrediction",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "oliomericPredPlot",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "outlierDetection",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "outlierDetectionPlot",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "signatureResults",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "mlFilteredSen",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "mlFiltering",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "noOligomericResults",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "noFullTextResults",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "publicationInfo",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "summaryResults",      selected = "")
  #        updateCheckboxGroupInput(session, inputId = "mutantResults",      selected = "")
  #        updateCheckboxInput(session, "deselectAll", label = "Deselect All", TRUE)
  #        updateCheckboxInput(session, "selectAll",   label = "Select All",   FALSE)
  #
  #    }
  #})
  
  
  ##### display data in data upload tab ########
  
  output$RawData <- DT::renderDataTable(
    if(input$publicationInfo){
      datatable(currentData(), escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
        dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf')#, colVis = list(activate = 'click', showAll = TRUE)
        
      ))
      
    })
  
  ## returns input data
  
  
  
  ##### oligomeric keyword list ########
  
  oligomerNames <- reactive({
    
    oligomerNames = c("monomer","monomeric",
                      "dimer","dimeric","homodimer", "heterodimer", "monomerofdimers",
                      "trimer","trimeric", "homotrimer", "heterotrimer",
                      "tetramer","tetrameric", "homotetramer", "heterotetramer", "dimerofdimers",
                      "pentamer","pentameric","homopentamer","heteropentamer",
                      "hexamer","hexameric","homohexamer", "heterohexamer", "trimerofdimers",
                      "heptamer","heptameric", "homoheptamer", "heteroheptamer",
                      "octamer", "octameric","homooctamer", "heterooctamer", "tetramerofdimers",
                      "nonamer", "nonameric", "homononamer", "heterononamer",
                      "decamer", "decameric","homodecamer", "heterodecamer",  "pentamerofdimers",
                      "undecamer","undecameric","homoundecamer", "heteroundecamer",
                      "dodecamer","dodecameric", "homododecamer", "heterododecamer",  "hexamerofdimers",
                      "tridecamer","tridecameric", "homotridecamer", "heterotridecamer",
                      "tetradecamer","tetradecameric", "homotetradecamer", "heterotetradecamer", "heptamerofdimers",
                      "pentadecamer", "pentadecameric","homopentadecamer", "heteropentadecamer",
                      "hexadecamer", "hexadecamerix", "homohexadecamer", "heterohexadecamer", "octamerofdimers",
                      "heptadecamer","heptadecameric", "homoheptadecamer", "heteroheptadecamer",
                      "octadecamer", "octadecameric", "homooctadecamer", "heterooctadecamer", "nonamerofdimers",
                      "nonadecamer","nonadecameric","homononadecamer", "heterononadecamer",
                      "eicosamer","eicosameric", "homoeicosamer", "heteroeicosamer", "decamerofdimers",
                      "undecamerofdimers",
                      "dodecamerofdimers",
                      "tridecamerofdimers",
                      "21mer",
                      "22mer",
                      "23mer",
                      "24mer",
                      "25mer",
                      "26mer",
                      "27mer",
                      "28mer",
                      "29mer",
                      "30mer")
  })
  
  
  ##### experimental evidence list ########
  
  experimentalEvidenceWordList <- reactive({
    
    
    experimentalEvidence = c("size-exclusion chromatography","size exclusion chromatography", "Size-exclusion chromatography","Size exclusion chromatography", "Size-Exclusion Chromatography","Size Exclusion Chromatography", "gel filtration","gel-filtration", "Gel filtration","Gel-filtration", "Gel Filtration","Gel-Filtration","affinity chromatography", "Affinity chromatography", "Affinity Chromatography", "electrophoresis", "Electrophoresis", "gel electrophoresis", "Gel electrophoresis", "Gel Electrophoresis","Native PAGE", "native PAGE","active site titration", "Active site titration", "Active Site Titration","electron microscopy", "Electron microscopy", "Electron Microscopy","mutagenesis", "Mutagenesis", "inhibition kinetics", "Inhibition kinetics", "Inhibition Kinetics", "titration", "Titration", "analytical centrifugation", "Analytical centrifugation", "Analytical Centrifugation","analytical ultracentrifuge", "Analytical ultracentrifuge", "Analytical Ultracentrifuge", "analytical ultracentrifugation", "Analytical ultracentrifugation", "Analytical Ultracentrifugation", "static light scattering", "Static light scattering", "Static Light Scattering", "static-light scattering", "Static-light scattering", "light-scattering","light scattering", "Light scattering","Light-scattering", "Light-Scattering", "dynamic light scattering", "dynamic light scattering", "dynamic-light scattering", "Dynamic light scattering", "Dynamic-light scattering","Dynamic Light Scattering", "Dynamic-Light Scattering", "sedimentation equilibrium", "sedimentation-equilibrium", "Sedimentation equilibrium", "Sedimentation-equilibrium", "Sedimentation-Equilibrium", "Sedimentation Equilibrium", "Small-Angle X-Ray Scattering", "Small-angle X-Ray Scattering", "small-angle X-ray scattering","Small Angle X Ray Scattering", "small-angle X-ray scattering", "small-angle x-ray scattering", "Small Angle X Ray Scattering", "small-angle x-ray scattering", "small angle x ray scattering", "small angle x-ray scattering", "nondenaturing gel","nondenaturing-gel", "Nondenaturing gel", "Nondenaturing-gel","SEC", "AUC","SLS","DLS","CCL","FRET","DSC","SV","NMR","MALS","LS","CL","EM","UDP","SAXS", "elusion profiles", "Elusion profiles", "Elusion Profiles", "elution profile", "Elution profiles", "Elution Profiles", "elution profiles", "Elution profiles", "Elution Profiles")
    
  })
  
  
  mutantProtein <- reactive({
    
    
    mutantProtein = c("mutagenesis", "mutagenic", "interface-disruption", "interface disruption", "disruption", "Mutagenesis", "Mutagenic", "Interface-disruption", "Interface disruption", "Disruption")
    
  })
  
  
  ##### help text in data upload tab ########
  
  alert <- reactive({
    if(!input$startAnalysis) {
      createAlert(session, "help", "popupHelp", title = "Welcome to Quaternary Structure Evaluation Tool (v.0.99)",
                  content = "This web tool is developed for evaluation of quaternary structures of protein structures in the Protein Data Bank. There are four different methods for the evaluation process:
                  <ul>
                  <li>Sequence cluster (SC) result with a consistency score for each PDB ID.</li> <br>
                  <li>Text mining of the primary citations in order to extract oligomeric state and experimental evidence.</li> <br>
                  <li>Stoichiometry and symmetry predictions for quaternary structures using PISA.</li><br>
                  <li>Stoichiometry and symmetry predictions for quaternary structures using EPPIC.</li>
                  </ul><br>
                  
                  In order to use this application:
                  <ul>
                  <li>Enter a <b>PDB ID</b>, <b>upload a file</b> which contains multiple PDB IDs or <b>upload PDF</b> version of a paper.</li> <br>
                  <li>Click <b>Start analysis</b> button to get the consensus result for both oligomeric state and symmetry.</li> <br>
                  <li>Click, <b>SC</b>, <b>Text mining</b>, <b>PISA</b> and <b>EPPIC</b> to see individual results.</li> <br>
                  <li>For advanced results and options click <b>Results</b> and <b>Options</b> checkboxes</li><br>
                  <li>For further information please see the <b>Help</b> tab.</li>
                  </ul>",
                  append = FALSE)
    }
    else {
      closeAlert(session, "popupHelp")
    }
  })
  ##### help text in outlier detection tab ########
  
  alertOD <- reactive({
    if(!input$startAnalysis) {
      createAlert(session, "helpOD", "popupHelpOD", title = "Welcome to Quaternary Structure Evaluation Pipeline!",
                  content = "<ul>
                  <li>First, select a sequence cluser: 95%, 90%, 70% or 40%.</li> <br>
                  <li>Then, click <b>Detect</b> button to detect outlier quaternary structure predictions.</li> <br>
                  <li>For further information please see the <b>Help</b> tab.</li>
                  </ul>",
                  append = FALSE)
    }
    else {
      closeAlert(session, "popupHelpOD")
    }
  })
  
  
  ##### help text in text mining tab ########
  
  alertTM <- reactive({
    
    if(!input$startAnalysis) {
      createAlert(session, "helpTM", "popupHelpTM", title = "Welcome to Quaternary Structure Evaluation Pipeline!",
                  content = "<ul>
                  <li>Select one of the machine learning algorithms from the left panel and click submit to start mining.</li> <br>
                  <li> A summary table will be displayed under the <b>Summary</b> tab.</li> <br>
                  <li> Sentences which are related with the quaternary structure will be displayed under the <b>Sentences</b> tab.</li> <br>
                  <li> To eliminate redundat sentences, these sentences will be classified as quaternary structure related (TRUE) and quaternary structure unrelated (FALSE) under the <b>Machine learning</b> tab. A majority probability will be calculated and oligomeric state prediction will be made for each PDB entry based on this calculated majority probability.</li> <br>
                  <li> The papers without any oligomeric state information will be displayed under the <b>No oligomeric state</b> tab.</li> <br>
                  <li> The papers without full text or abstracts will be displayed under the <b>No full text/abstract result</b> tab.</li> <br>
                  <li>For further information please see the <b>Help</b> tab.</li>
                  </ul>",
                  append = FALSE)
    }
    else {
      closeAlert(session, "popupHelpTM")
    }
  })
  
  
  ##############################################################################################
  ############### PMC Text Miner - Start #######################################################
  ##############################################################################################
  
  TextMinerResult <- eventReactive(input$startAnalysis,{
    
    if(input$OligomericStatePrediction){
      oligomerNames = oligomerNames()
      mutantProtein = mutantProtein()
    
      newTextMiner(currentData()[,2], as.character(currentData()[,1]), as.character(currentData()[,3]), as.character(currentData()[,4]), oligomerNames, mutantProtein, fi =  dataPath())
    }
  })
  
  textMiner <- reactive({
    
    if(input$OligomericStatePrediction){
    
    
    if(input$startAnalysis == 0)
    {
      return()
    }
    
    isolate({
      if(input$startAnalysis)
        withProgress(message = 'Evaluation started...',  detail = 'Please wait...', value=7,{
          
          miningResults = TextMinerResult()
          
          
          if(input$startAnalysis){
            
            fullResult = miningResults[[1]][!(miningResults[[1]][,"oligomericState"] == "Inconclusive"),]
            
            if(dim(fullResult)[1] != 0){
              
              sentencesList = fullResult[4]
              
              wordList = data.frame(matrix(NA,1))
              names(wordList) = c("words")
              
              for(i in 1:dim(sentencesList)[1]){
                review_source <- tm::VectorSource(sentencesList[i,1])
                corpus <- tm::Corpus(review_source)
                corpus <- tm::tm_map(corpus, tm::stripWhitespace) # remove extra white spaeces
                corpus <- tm::tm_map(corpus, tm::content_transformer(tolower)) # convert lower cases
                corpus <- tm::tm_map(corpus, tm::removePunctuation) # remove punctuation
                corpus <- tm::tm_map(corpus, tm::removeWords, stopwords("english")) # remove stop words
                #corpus <- tm_map(corpus, stemDocument)
                dtm <- as.data.frame(as.matrix(tm::DocumentTermMatrix(corpus)))
                wordList[i,1] = gsub("\n", " ", NLP::as.String(colnames(dtm)))
              }
              
              
              test <- hashed.model.matrix(~ split(words, delim = " ", type = "tf-idf"),
                                          data = wordList, hash.size = hashSize, signed.hash = FALSE)
              
              
              if(input$MLalgorithm == "Support vector machines"){
                
                prediction = predict.train(svm, newdat=test)
                
              }
              
              if(input$MLalgorithm == "Boosted logistic regression"){
                
                prediction = predict.train(lr, newdat=test)
                
              }
              
              if(input$MLalgorithm == "Random forest"){
                
                prediction = predict.train(rf, newdat=test)
                
              }
              
              #lrpred = predict.train(lr, test)
              Label = ifelse(prediction=="yes","TRUE","FALSE")
              ml = as.data.frame(cbind(fullResult, Label))
              mlResult = ml[c(1:5)]
              names(mlResult) = c("PDB ID", "PMC ID", "Keyword","Sentence", "Label")
              
              mlResult = ml
              
              #mlResult2 = pmc[[1]]
              
              
              if(!is.null(unlist(mlResult))){
                
                for(i in 1:length(oligomericStateList)){
                  mlResult$oligomericState = gsub(paste(oligomericStateList[[i]], collapse="|"), oligomericStateList[[i]][length(oligomericStateList[[i]])], mlResult$oligomericState)
                }  
                
                
                mlResult[,"newOligomericState"] = mlResult$oligomericState
                
                experimentalList = list()
                experimentalEvidenceLast = list()
                
                
                if(length(unique(mlResult$pdbId)) > 1){
                  
                  sentenceSplit = split(mlResult, mlResult$pdbId)
                  
                }else{
                  
                  sentenceSplit = list(mlResult$oligomericSentence)
                }
                
                
                
                for(senSplt in 1:length(sentenceSplit)){
                  
                  if(length(unique(mlResult$pdbId)) > 1){
                    sentence = removeNumbers(stripWhitespace(sentenceSplit[[senSplt]]$oligomericSentence))
                  }else{
                    
                    sentence = removeNumbers(stripWhitespace(sentenceSplit[[senSplt]]))
                  }
                  
                  
                  experimentalEvidenceWordList = list()
                  
                  
                  for(sen in 1:length(sentence)){
                    
                    
                    experimentalWords = as.data.frame(str_locate_all(pattern =paste(experimentalEvidenceWordList(), collapse="\\b|\\b"),sentence[sen]))
                    
                    if(!is.na(experimentalWords[1,][1])){
                      
                      sameSentenceList=list()
                      
                      for(exp in 1:dim(experimentalWords)[1]){
                        
                        sameSentenceList[[exp]]=substr(sentence[sen], experimentalWords[exp,][1], experimentalWords[exp,][2])
                        
                        
                      }
                      
                      experimentalEvidenceWordList[[sen]] = paste(unlist(sameSentenceList), collapse = ", ")
                      
                      
                    }else{
                      
                      experimentalEvidenceWordList[[sen]] = "not found"
                    }
                  }
                  
                  expr = as.data.frame(unlist(experimentalEvidenceWordList))
                  
                  experimentalEvidenceLast[[senSplt]] = cbind(as.data.frame(sentenceSplit[[senSplt]]),expr)
                  
                }
                
                
                experimentalEvidenceLast2 = as.data.frame(rbindlist(experimentalEvidenceLast))
                
                # mlResult$experimentalEvidence = experimentalEvidenceLast2[,"unlist(experimentalEvidenceWordList)"]
                
                # names(mlResult)[8] = "experimentalEvidence"
                
                
                
                
                if(length(unique(mlResult$pdbId)) == 1){
                  
                  mlResult$experimentalEvidence = experimentalEvidenceLast2[,"unlist(experimentalEvidenceWordList)"]
                  
                }
                if(length(unique(mlResult$pdbId)) > 1){
                  
                  mlResult = experimentalEvidenceLast2
                }
                
                names(mlResult)[dim(mlResult)[2]] = "experimentalEvidence"
                
                
                
                
                
                
                MLlist = list()
                
                MLfilter = split(mlResult, mlResult$pdbId)
                
                for(p in 1:length(MLfilter)){
                  
                  
                  if("TRUE"%in%MLfilter[[p]]$Label){
                    MLlist[[p]] = MLfilter[[p]][MLfilter[[p]]$Label=="TRUE",]
                    MLlist[[p]]$MLfilter="TRUE"
                    
                  }else{
                    MLlist[[p]] = MLfilter[[p]]
                    MLlist[[p]]$MLfilter="FALSE"
                    
                  }
                  
                }
                
                
                mlResult2 = as.data.frame(rbindlist(MLlist))
                
                evidenceSplitList = list()
                
                evidenceSplit = split(mlResult2, mlResult2$pdbId)
                
                for(evidence in 1:length(evidenceSplit)){
                  
                  
                  if("TRUE" %in% grepl(paste(experimentalEvidenceWordList(), collapse="|"), evidenceSplit[[evidence]]$experimentalEvidence)){
                    evidenceSplitList[[evidence]] = evidenceSplit[[evidence]][evidenceSplit[[evidence]]$experimentalEvidence!="not found",]
                    
                  }else{
                    
                    evidenceSplitList[[evidence]] = evidenceSplit[[evidence]]
                  }
                }
                
                
                mlResult2 = as.data.frame(rbindlist(evidenceSplitList))
                
                pdbID = split(mlResult2, mlResult2$pdbId)
                
                
                resultOligomericState = list()
                
                for(z in 1:length(pdbID)) {
                  
                  pubSubset = pdbID[z]
                  nameOligomericState = factor(pubSubset[[1]]$newOligomericState)
                  levelsOligomericState = levels(nameOligomericState)
                  tbl = table(nameOligomericState)
                  rowNamesOligomericState = rownames(tbl)
                  probabilityOligomericState = tbl/sum(tbl)
                  for(i in 1:dim(pubSubset[[1]])[1]){
                    
                    for(j in 1:length(levelsOligomericState)){
                      
                      if(nameOligomericState[i] == rowNamesOligomericState[j]) {
                        
                        pubSubset[[1]]$probability[i] = probabilityOligomericState[j]
                      }
                    }
                  }
                  resultOligomericState[z] = pubSubset
                }
                
                ### find unique probabilites of each pdb entry ###
                
                data = as.data.frame(rbindlist(resultOligomericState))
                
                
                if(length(unique(data$pdbId))>1){
                  
                  splitData = split(data, data$pdbId)
                  
                  resultUniqueValues = list()
                  
                  
                  for(k in 1:length(splitData)) {
                    
                    subset = splitData[k][[1]]
                    
                    resultUniqueValues[k][[1]] = unique(subset)
                    
                  }
                  
                  resultUniqueValues = as.data.frame(rbindlist(resultUniqueValues))
                  resultUniqueValues = resultUniqueValues[,c("PubId", "pdbId", "newOligomericState", "probability", "MLfilter", "experimentalEvidence", "pubType", "primaryCitation")]
                  
                }else{
                  
                  splitData = data
                  
                  
                  resultUniqueValues = unique(splitData[,c("PubId", "pdbId", "newOligomericState", "probability", "MLfilter", "experimentalEvidence", "pubType", "primaryCitation")])
                  
                  
                }
                
                ### find only one probability of each pdb entry ###
                
                
                uniqueData = unique(resultUniqueValues)
                
                
                splitUniqueData = split(uniqueData, uniqueData$pdbId)
                
                resultUniqueForEach = list()
                
                for(k in 1:length(splitUniqueData)) {
                  a = splitUniqueData[[k]]
                  resultUniqueForEach[[k]] = a[which(a$probability == max(a$probability)),]
                  
                  
                }
                
                
                resultUniqueProbabilityList= as.data.frame(rbindlist(resultUniqueForEach))
                
                
                uniqueData2 = unique(resultUniqueProbabilityList)
                
                
                #            splitUniqueData2 = split(uniqueData2, uniqueData2$pdbId)
                
                #            resultUniqueForEach2 = list()
                
                #            for(t in 1:length(splitUniqueData2)) {
                #                reliabilityResult = splitUniqueData2[[t]]
                #                #publicationType = if(reliabilityResult$pubType[1]== "full text")(1)else(0)
                #                #pCitation = if(reliabilityResult$primaryCitation[1]== "yes")(2)else(0)
                #                expEvidence = if(reliabilityResult$experimentalEvidence[1]== "not found")(0)else(1)
                #                mlFiltering = if(reliabilityResult$MLfilter[1]== "filtered")(1)else(0)
                #
                #                reliability = sum(expEvidence,mlFiltering)
                #
                #                if(reliability == 0){predictionReliability = "low"}
                #                if(reliability == 1){predictionReliability = "moderate"}
                #                if(reliability == 2){predictionReliability = "high"}
                #                #if(reliability>=5){predictionReliability = "very high"}
                #                reliabilityResult$predictionReliability = predictionReliability
                #
                #                resultUniqueForEach2[[t]] = reliabilityResult
                #
                #
                #            }
                
                resultUniqueProbabilityList= uniqueData2
                
                
                resultUniqueProbabilityList[,4] = format(round(resultUniqueProbabilityList[,4],2), nsmall = 2)
                
                
                oligomericSentence = mlResult[,c("pdbId","PubId", "newOligomericState","experimentalEvidence", "oligomericSentence")]
                names(oligomericSentence) = c("PDB ID", "Publication ID", "Keyword", "Experimental evidence", "Sentence")
                
              }
              
              
              
              
              resultUniqueProbabilityList2 = unique(resultUniqueProbabilityList[,-6])
              
              if(dim(resultUniqueProbabilityList)[1]!=dim(resultUniqueProbabilityList2)[1]){
                
                
                for(i in 1:dim(resultUniqueProbabilityList2)[1]){
                  
                  keyword = as.character(resultUniqueProbabilityList$experimentalEvidence[which(resultUniqueProbabilityList2$pdbId[i] == resultUniqueProbabilityList$pdbId)])
                  
                  if(length(keyword) > 1){
                    
                    keyword = paste(keyword, collapse = ", ")
                    
                  }
                  
                  resultUniqueProbabilityList2$experimentalEvidence[i] =  keyword
                  
                }
                
                
                names(resultUniqueProbabilityList2) = c("Publication ID", "PDB ID","Oligomeric state", "Probability", "ML Filtering", "Publication type", "Primary citation",
                                                        "Experimental evidence")
                
                resultUniqueProbabilityList2 = resultUniqueProbabilityList2[,c("PDB ID", "Oligomeric state", "Experimental evidence", "ML Filtering", "Probability", "Publication ID", "Publication type", "Primary citation")]
                
                
              }else{
                
                resultUniqueProbabilityList2 = resultUniqueProbabilityList
                
                names(resultUniqueProbabilityList2) = c("Publication ID", "PDB ID","Oligomeric state", "Probability", "ML Filtering",  "Experimental evidence", "Publication type", "Primary citation")
                
                resultUniqueProbabilityList2 = resultUniqueProbabilityList2[,c("PDB ID", "Oligomeric state", "Experimental evidence", "ML Filtering", "Probability", "Publication ID", "Publication type", "Primary citation")]
                
                
                
              }
              
              resultUniqueProbabilityList = resultUniqueProbabilityList2
              
              
              mutantTable = miningResults[[11]]
              
              mutantTable = mutantTable[,c("pdbId","PubId","mutantKeyword","mutantSentence")]
              names(mutantTable) = c("PDB ID", "Publication ID", "Keyword", "Sentence")
              
              
              #   for(i in 1:dim(resultUniqueProbabilityList)[1]){
              #
              #    if(resultUniqueProbabilityList$`PDB ID`[i] %in% mutantTable$`PDB ID`){
              #
              #        resultUniqueProbabilityList$probableMutant[i] = "yes"
              #
              #    }else{
              #
              #        resultUniqueProbabilityList$probableMutant[i] = "no"
              #
              #    }
              #
              #    }
              
              
              for(i in 1:dim(resultUniqueProbabilityList)[1]){
                
                
                if(substring(resultUniqueProbabilityList$`PDB ID`[i], 88, 91) %in% currentData()[,1]){
                  
                  resultUniqueProbabilityList$relatedStructures[i] = currentData()[which(substring(resultUniqueProbabilityList$`PDB ID`[i], 88, 91) == currentData()[,1]),5]
                  
                }
                
              }
              
              
              names(resultUniqueProbabilityList)[9] = c("Related structures")
              
              noResult = miningResults[[1]][(miningResults[[1]][,"oligomericState"] == "Inconclusive"),]
              
              
              
              
              if(dim(noResult)[1] != 0){
                
                
                noResult$experimentalEvidence = "Inconclusive"
                noResult$MLfiltering = "Inconclusive"
                noResult$Probability = "Inconclusive"
                #noResult$PredictionReliability = "Inconclusive"
                #noResult$ProbableMutant = "Inconclusive"
                noResult$RelatedStructures = "Inconclusive"
                
                names(noResult) = c("PDB ID", "Publication ID", "Oligomeric state", "oligomericSentence", "Publication type",
                                    "Primary citation", "Experimental evidence", "ML Filtering", "Probability",
                                    "Related structures")
                
                noResultLast = noResult[,c("PDB ID","Oligomeric state","Experimental evidence","ML Filtering","Probability",
                                           "Publication ID","Publication type","Primary citation",
                                           "Related structures")]
                
                
                resultCombined = rbind(resultUniqueProbabilityList, noResultLast)
                
              }else{
                
                resultCombined = resultUniqueProbabilityList
              }
              
              
              mlResult = mlResult[,c("pdbId", "PubId","oligomericState","experimentalEvidence","oligomericSentence", "Label")]
              names(mlResult) = c("PDB ID", "Publication ID", "Keyword", "Experimental evidence","Sentence", "Label")
              
              mlFiltered = mlResult[mlResult$Label == "TRUE",]
              
              summaryList = list(miningResults$nofFoundedArticles,miningResults$nofPMCArticles,miningResults$nofPubMedArticles,miningResults$noOligomericState, miningResults$noFullText, miningResults$numberOfuniquePdbId)
              
              summaryResult=data.frame(do.call(rbind, summaryList))
              
              summaryResult = cbind(NA,summaryResult)
              summaryResult[1,1] = "Total number of articles"
              summaryResult[2,1] = "Number of PMC articles"
              summaryResult[3,1] = "Number of PubMed articles"
              summaryResult[4,1] = "Number of articles with no oligomeric state"
              summaryResult[5,1] = "Number of articles with no full text/abstract"
              summaryResult[6,1] = "Number of unique PDB IDs"
              summaryResult[7,1] = "Machine learning algorithm"
              summaryResult[7,2] = input$MLalgorithm
              names(summaryResult) = c("Source", "Count")
              
              
              
            }else{
              
              
              if(dim(miningResults[[1]])[1] !=0){
                noResult = miningResults[[1]][(miningResults[[1]][,"oligomericState"] == "Inconclusive"),]
              }else{
                
                noResult = data.frame(matrix(NA,1,6))
                names(noResult) = c("PDB ID", "Publication ID", "Oligomeric state", "oligomericSentence", "Publication type",
                                    "Primary citation")
                noResult$`PDB ID` = currentData()[i,1]
                noResult$`Publication ID` = "Inconclusive"
                noResult$`Oligomeric state` = "Inconclusive"
                noResult$oligomericSentence = "Inconclusive"
                noResult$`Publication type` = "Inconclusive"
                noResult$`Primary citation` = "Inconclusive"
                
              }
              
              noResult$experimentalEvidence = "Inconclusive"
              noResult$MLfiltering = "Inconclusive"
              noResult$Probability = "Inconclusive"
              #noResult$PredictionReliability = "Inconclusive"
              #noResult$ProbableMutant = "Inconclusive"
              noResult$RelatedStructures = "Inconclusive"
              
              names(noResult) = c("PDB ID", "Publication ID", "Oligomeric state", "oligomericSentence", "Publication type",
                                  "Primary citation", "Experimental evidence", "ML Filtering", "Probability",
                                  "Related structures")
              
              noResultLast = noResult[,c("PDB ID","Oligomeric state","Experimental evidence","ML Filtering","Probability",
                                         "Publication ID","Publication type","Primary citation",
                                         "Related structures")]
              
              
              resultCombined = noResultLast
              
              
              
              
              
              
              mlResult = NULL
              #names(mlResult) = c("PDB ID", "Publication ID", "Keyword", "Experimental evidence","Sentence", "Label")
              
              mlFiltered = NULL
              
              summaryList = list(miningResults$nofFoundedArticles,miningResults$nofPMCArticles,miningResults$nofPubMedArticles,miningResults$noOligomericState, miningResults$noFullText, miningResults$numberOfuniquePdbId)
              
              summaryResult=data.frame(do.call(rbind, summaryList))
              
              summaryResult = cbind(NA,summaryResult)
              summaryResult[1,1] = "Total number of articles"
              summaryResult[2,1] = "Number of PMC articles"
              summaryResult[3,1] = "Number of PubMed articles"
              summaryResult[4,1] = "Number of articles with no oligomeric state"
              summaryResult[5,1] = "Number of articles with no full text/abstract"
              summaryResult[6,1] = "Number of unique PDB IDs"
              summaryResult[7,1] = "Machine learning algorithm"
              summaryResult[7,2] = input$MLalgorithm
              names(summaryResult) = c("Source", "Count")
              
              
              oligomericSentence = NULL
              mutantTable = NULL
              
              
              
              
            } 
            
            noOligomericResult = miningResults[[2]]
            noFullText = miningResults[[3]]
            names(noOligomericResult) = names(noFullText) = c("PDB ID", "Publication ID")
            
            result = list(resultUniqueProbabilityList = resultCombined, oligomericSentence = oligomericSentence, noOligomericResult = noOligomericResult, noFullText = noFullText, summaryResult = summaryResult, allSentences = mlResult, mlFiltered = mlFiltered, mutantInfo = mutantTable)
            result
            
            
          }
          
        })})
    
  }
    
    })
  
  
  ############################################################################################
  ############### Text mining  - Start ######################################################
  ############################################################################################
  output$TextMining <- DT::renderDataTable({
    
    if(input$OligomericStatePrediction){
    alert()
    
    if (!input$OligomericStatePrediction){
      return()
      
    }
    
    else{
      
      datatable(textMiner()[[1]],escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
        dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf'),   deferRender = TRUE
        
      ))
      
      
    }
    
    }
  })
  ############################################################################################
  ############### Text mining - End ####################################################
  ############################################################################################
  
  
  
  ############################################################################################
  ############### Sentences - Start ####################################################
  ############################################################################################
  output$sentences<- DT::renderDataTable({
    
    if(input$OligomericStatePrediction){
    
      validate(textMiner()[[2]])
    
    }
  })
  ############################################################################################
  ############### Sentences - End ####################################################
  ############################################################################################
  
  ############################################################################################
  ############### Machine learning - Start ####################################################
  ############################################################################################
  output$machineLearning<- DT::renderDataTable({
    
    if(input$OligomericStatePrediction){
        
        if (!input$mlFiltering){
          return()
          
        }
        
        if(input$mlFiltering){
          
          data = textMiner()[[6]]
          datatable(data, escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
            dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf')#, colVis = list(activate = 'click', showAll = TRUE)
            
          ))%>% formatStyle(
            "Sentence", "Label",
            backgroundColor = styleEqual(c("FALSE", "TRUE"), c('#F8766D', '#00BA38')))
          
        }
    }
  })
  ############################################################################################
  ############### Machine learning - End ####################################################
  ############################################################################################
  
  ############################################################################################
  ############### ML filtered sentences - Start ####################################################
  ############################################################################################
  output$mlFilteredSentences<- DT::renderDataTable({
    
    if(input$OligomericStatePrediction){
      
      if (!input$mlFilteredSen){
        return()
        
      }
      
      if(input$mlFilteredSen){
        datatable(textMiner()[[7]],escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
          dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf'),   deferRender = TRUE
          
        ))
      }
    }
  })
  ############################################################################################
  ############### ML filtered sentences - End ####################################################
  ############################################################################################
  
  ############################################################################################
  ############### Mutant Info - Start ####################################################
  ############################################################################################
  
  output$mutantResults<- DT::renderDataTable({
    
    if(input$OligomericStatePrediction){
        if (!input$mutantResults){
          return()
          
        }
        
        if(input$mutantResults){
          
          
          datatable(textMiner()[[8]],escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
            dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf'),   deferRender = TRUE
            
          ))
        }
    }
  })
  
  
  ############################################################################################
  ############### Mutant Info - End ####################################################
  ############################################################################################
  
  
  
  
  ############################################################################################
  ############### No oligomeric result learning - Start ######################################
  ############################################################################################
  
  
  output$noOligomericResult<- DT::renderDataTable({
    
    if(input$OligomericStatePrediction){
        
        if (!input$noOligomericResults){
          return()
          
        }
        
        
        if(!length(textMiner()[[3]])==0 && input$noOligomericResults){
          
          
          
          datatable(textMiner()[[3]],escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
            dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf'),   deferRender = TRUE
            
          ))
          
        }
        
    }
  })
  
  
  ############################################################################################
  ############### No oligomeric result learning - Start ######################################
  ############################################################################################
  
  
  output$noFullText<- DT::renderDataTable({
    
    if(input$OligomericStatePrediction){
    
        if (!input$noFullTextResults){
          return()
          
        }
        
        
        if(!length(textMiner()[[4]])==0 && input$noFullTextResults){
          
          
          
          datatable(textMiner()[[4]], escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
            dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf'),   deferRender = TRUE
            
          ))
          
        }
    
    }  
    
  })
  
  ############################################################################################
  ############### No oligomeric result learning - End ######################################
  ############################################################################################
  
  ############################################################################################
  ############### Summary - Start ######################################
  ############################################################################################
  
  output$summary<- DT::renderDataTable({
    
    
    if(input$OligomericStatePrediction){
        if (!input$summaryResults){
          return()
          
        }
        
        
        if(!length(textMiner()[[5]])==0 && input$summaryResults){
          
          
          
          
          datatable(textMiner()[[5]],escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
            dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf'),   deferRender = TRUE
            
          ))
          
        }
    }
  },escape=FALSE, rownames=FALSE)
  
  ############################################################################################
  ############### Summary - End ######################################
  ############################################################################################
  
  ############################################################################################
  ############### Help - Model performances  - Start ######################################
  ############################################################################################
  output$testPerformance<- renderPrint({
    
    
    TestSetPerformance
    
    
    
  })
  
  
  #output$testPerformanceROC<- renderPlot({
  
  #svmROC = plot.roc(testSet$testClass, svmPerformance$svmpred2[,"yes"], col="red")
  #lrROC = lines.roc(testSet$testClass,lrPerformance$lrpred2[,"yes"], col="blue")
  #rfROC = lines.roc(testSet$testClass, rfPerformance$rfpred2[,"yes"], col="#008600")
  #legend("bottomright", legend=c("SVM", "BLR", "RF"), col=c("red", "blue", "#008600"), lwd=2)
  #},  height = 500, width = 500)
  
  
  
  
  
  output$baKeywords<- renderPrint({
    
    
    keywords = oligomerNames()
    print(keywords)
    
    
  })
  
  output$eeKeywords<- renderPrint({
    
    
    keywords = experimentalEvidenceWordList()
    print(keywords)
    
    
  })
  
  ############################################################################################
  ############### Help - Model performances  - End ######################################
  ############################################################################################
  
  
  outlierDetect <- reactive({
    
    if(input$seqCluster == "95"){
      
      # 95% sequence cluster
      index95 = which(consistencyData95$"PDB ID" %in% as.character(currentData()[,1]))
      if(length(index95) != 0){
        result = consistencyData95[index95, c("PDB ID",	"BA Number", "Stoichiometry", "Symmetry","Representative", "Consistency score", "Result")]
        
        result
        
      }
      
      
    }
    
    else if(input$seqCluster == "90"){
      
      # 90% sequence cluster
      index90 = which(consistencyData90$"PDB ID" %in% as.character(currentData()[,1]))
      result = consistencyData90[index90, c("PDB ID",	"BA Number", "Stoichiometry", "Symmetry","Representative", "Consistency score", "Result")]
      result
    }
    
    else if(input$seqCluster == "70"){
      
      # 70% sequence cluster
      index70 = which(consistencyData70$"PDB ID" %in% as.character(currentData()[,1]))
      result = consistencyData70[index70, c("PDB ID",	"BA Number", "Stoichiometry", "Symmetry","Representative", "Consistency score", "Result")]
      
      result
    }
    
    else if(input$seqCluster == "40"){
      
      # 40% sequence cluster
      index40 = which(consistencyData40$"PDB ID" %in% as.character(currentData()[,1]))
      result = consistencyData40[index40, c("PDB ID",	"BA Number", "Stoichiometry", "Symmetry","Representative", "Consistency score", "Result")]
      
      result
      
    }
    
  })
  
  
  outlierDetect2 <- reactive({
    
    
    result = outlierDetect()
    result$Result = as.character(result$Result)
    
    lastResult2 = sequenceCluster()
    
    result$nofPDB = NA
    
    resUnique = unique(lastResult2$Representative)
    
    for(i in 1:length(resUnique)){
      
      ind = lastResult2$Representative %in% resUnique[i]
      
      a = lastResult2[ind,]
      
      ind2 = result$Representative %in% resUnique[i]
      
      
      result[ind2,"nofPDB"] = sum(a$`# of PDBs`)
      
    }
    
    for(i in 1:dim(result)[1]){
      
      if(result$nofPDB[i] < 3){
        
        result$`Consistency score`[i] = "not enough chains"
        result$Result[i] = "no result"
      }
    }
    
    result = result[,c("PDB ID", "BA Number", "Stoichiometry", "Symmetry","Representative", "Consistency score", "Result")]
    
    result
    
  })
  
  
  ############################################################################################
  ############### Outlier detection  - Start ######################################
  ############################################################################################
  
  #     output$consistency<- DT::renderDataTable({
  #
  #
  #
  #         if(!input$signatureResults)
  #         {
  #             return()
  #         }
  #
  #
  #         #isolate({
  #         if(input$signatureResults && input$startAnalysis)
  #         #withProgress(message = 'Detecting outliers...',  detail = 'Please wait...', value=7,{
  #
  #             datatable(outlierDetect2()[-5],escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
  #             dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf'),   deferRender = TRUE)#, lengthMenu = c(5,50,34000))
  #
  # )
  
  #})
  
  ############################################################################################
  ############### Outlier detection  - End ######################################
  ############################################################################################
  outlierProcess <- reactive({
    
    result = outlierDetect()
    cols = c("Stoichiometry", "Symmetry")
    result$oligomericState <- apply(result[,cols ] , 1 , paste , collapse = "_" )
    names(result)[c(1,6)]=c("pdbId","prob")
    unique(result[,c("pdbId","prob","oligomericState")])
    
    
    
  })
  
  ############################################################################################
  ############### Outlier detection plot- Start ######################################
  ############################################################################################
  
  
  output$section1 <- renderText({
    if (input$startAnalysis && input$OligomericStatePrediction){
      'Text mining results'
    }
  })
  
  #output$section2 <- renderText({
  #    if (input$startAnalysis && input$oliomericPredPlot && !input$OligomericStatePrediction)
  #        {'Please click "Text mining" method on the left.'}
  #        else if (input$startAnalysis && input$oliomericPredPlot ){
  #            'Distribution plot for text mining results'
  #        }
  #})
  
  #output$section3 <- renderText({
  #    if (input$startAnalysis && input$mlFilteredSen){
  #        'Sentences filtered by machine learning algorithm	'
  #    }
  #})
  
  #output$section4 <- renderText({
  # if (input$startAnalysis && input$signatureResults){
  #
  #     if ( input$seqCluster == "95"){
  #         'Sequence cluster results (95% sequence cluster)'
  #     }
  #
  #     else if ( input$seqCluster == "90"){
  #         'Sequence cluster results (90% sequence cluster)'
  #     }
  #
  #     else if ( input$seqCluster == "70"){
  #         'Sequence cluster results (70% sequence cluster)'
  #     }
  #
  #     else if ( input$seqCluster == "40"){
  #         'Sequence cluster results (40% sequence cluster)'
  #     }
  # }
  #})
  
  #output$section5 <- renderText({
  #    if (input$startAnalysis && input$outlierDetectionPlot){
  #
  #        if ( input$seqCluster == "95"){
  #            'Distribution plot for sequence cluster (95% sequence cluster)'
  #        }
  #
  #        else if ( input$seqCluster == "90"){
  #            'Distribution plot for sequence cluster (90% sequence cluster)'
  #        }
  #
  #        else if ( input$seqCluster == "70"){
  #            'Distribution plot for sequence cluster (70% sequence cluster)'
  #        }
  #
  #        else if ( input$seqCluster == "40"){
  #            'Distribution plot for sequence cluster (40% sequence cluster)'
  #        }
  #    }
  #})
  
  
  output$section6 <- renderText({
    if (input$startAnalysis && input$mlFiltering){
      'All extracted sentences from the paper(s)'
    }
  })
  
  output$section7 <- renderText({
    if (input$startAnalysis && input$eppicPrediction){
      'EPPIC prediction'
    }
  })
  
  #output$section8 <- renderText({
  #    if (input$startAnalysis && input$noFullTextResults){
  #        'Papers with no full text/abstract'
  #    }
  #})
  
  #output$section9 <- renderText({
  #    if (input$startAnalysis && input$summaryResults){
  #        'Summary'
  #    }
  #})
  
  
  #output$section10 <- renderText({
  #    if (input$startAnalysis && input$publicationInfo){
  #        'Publication information'
  #    }
  #
  #})
  
  output$section11 <- renderText({
    if (input$startAnalysis && input$signatureResults){
      
      if ( input$seqCluster == "95"){
        '95% Sequence cluster results'
      }
      
      else if ( input$seqCluster == "90"){
        '90% Sequence cluster results'
      }
      
      else if ( input$seqCluster == "70"){
        '70% Sequence cluster results'
      }
      
      else if ( input$seqCluster == "40"){
        '40% Sequence cluster results'
      }
    }
    
  })
  
  
  
  output$section13 <- renderText({
    if (input$startAnalysis && input$pisaPrediction){
      'PISA prediction'
    }
  })
  
  output$section14 <- renderText({
    if (input$startAnalysis && input$Display){
      '3D view'
    }
  })
  
  
  output$section15 <- renderText({
    
    if (input$startAnalysis && input$advancedPisa){
      'PISA advanced results'
    }
  })
  
  output$section16 <- renderText({
    
    if (input$startAnalysis){
      'Oligomeric state'
    }
  })
  
  output$section17 <- renderText({
    if (input$startAnalysis){
      'Symmetry'
    }
    
  })
  
  #heightSize <- reactive({
  #
  #    if(input$startAnalysis == 0)
  #    {
  #        100
  #    }
  #    if(input$startAnalysis){
  #    resultHeight = outlierProcess()
  #    90*((length(unique(resultHeight$pdbId))^2)/(length(unique(resultHeight$pdbId))^1.4))
  #    }
  #})
  
  
  #heightSizeNULL <- reactive({
  #
  #    if(input$startAnalysis == 0)
  #    {
  #        100
  #    }
  #})
  
  #react<- reactiveValues(h = 200)
  
  
  #output$consistencyPlot<- renderPlot({
  #
  #    if(!input$outlierDetectionPlot)
  #    {
  #        return()
  #    }
  #        isolate({
  #        if(input$outlierDetectionPlot){
  #        withProgress(message = 'Creating plot...',  detail = 'Please wait...', value=7,{
  #
  #            result = outlierProcess()
  #
  #            result <- result %>%
  #            group_by(pdbId, oligomericState) %>%
  #            summarise(prob = sum(as.numeric(prob))) %>%
  #            mutate(csum = cumsum(as.numeric(prob)))
  #
  #            ggplot(result, aes(x = pdbId, y = as.numeric(prob), fill = oligomericState, order=prob)) +
  #            geom_bar(stat = "identity", fill = "#01DF3A")  + xlab("PDB ID") + ylab("Consistency score") +
  #            coord_flip() + theme_bw() +
  #            facet_grid(. ~ oligomericState, scales = "free_y", shrink = FALSE)+
  #            guides(fill=FALSE) + ylim(c(0,1)) + geom_hline(yintercept = 0.5, colour = "black")
  #            #scale_fill_manual(labels=c("A2_C1","A2_C2","AB_C1","AB_C2", "A8_D4"), values=c("#04B431","#04B431","#04B431","#04B431","red"))
  #
  #            }
  #            )
  #        }})
  #
  #},  height = function(){
  #
  #    if(input$startAnalysis && input$outlierDetectionPlot){
  #
  #    isolate({
  #        resultHeight = outlierProcess()
  #        90*((length(unique(resultHeight$pdbId))^2)/(length(unique(resultHeight$pdbId))^1.4))
  #
  #        }
  #    )}
  #
  #    else
  #    {
  #        "auto"
  #    }
  #}
  #)
  
  ############################################################################################
  ############### Outlier detection plot- End ######################################
  ############################################################################################
  
  
  
  ############################################################################################
  ############### Oligomeric state prediction plot- Start ######################################
  ############################################################################################
  
  
  #output$ospPlot<- renderPlot({
  #
  #    if(!input$oliomericPredPlot)
  #    {
  #        return()
  #    }
  #
  #
  #if(input$OligomericStatePrediction){
  #
  #            result = textMiner()[[1]]
  #               
  #               #c("PDB ID", "Oligomeric state", "Experimental evidence", "ML Filtering", "Probability", "Publication ID", "Publication type", "Primary citation")]
  #
  #            result = unique(result[,c("PDB ID","Oligomeric state","Probability")])
  #
  #            names(result) = c("pdbId","oligomericState", "prob")
  #           result$pdbId = substr(result$pdbId, 88, 91)
  #
  #            result <- result %>%
  #            group_by(pdbId, oligomericState) %>%
  #            summarise(prob = sum(as.numeric(prob))) %>%
  #            mutate(csum = cumsum(as.numeric(prob)))
  #
  #            ggplot(result, aes(x = pdbId, y = as.numeric(prob), fill = oligomericState, order=prob)) +
  #            geom_bar(stat = "identity", fill = "#01DF3A")  + xlab("PDB ID") + ylab("Probability") +
  #            coord_flip() + theme_bw() +
  #            facet_grid(. ~ oligomericState, scales = "free_y", shrink = FALSE)+
  #            guides(fill=FALSE) + ylim(c(0,1))
  #}
  #
  #}, height = function(){
  #
  #    if(input$startAnalysis && input$oliomericPredPlot){
  #
  #        isolate({
  #
  #                result = unique(textMiner()[[1]])
  #                90*((length(unique(result$'PDB ID'))^2)/(length(unique(result$'PDB ID'))^1.4))
  #
  #            })
  #    }
  #
  #    else
  #    {
  #         "auto"
  #
  #    }
  #}
  #)
  
  
  
  ############################################################################################
  ############### Oligomeric state prediction plot- End ######################################
  ############################################################################################
  
  
  
  ############################################################################################
  ############### Sequence cluster - Start ######################################
  ############################################################################################
  
  sequenceCluster <- reactive({
    
    pdbList = outlierDetect()
    
    if(input$seqCluster == "95"){
      
      sequenceCluster = consistencyData95[consistencyData95[,3] %in% pdbList$Representative,]
      
    }
    
    else if(input$seqCluster == "90"){
      
      sequenceCluster = consistencyData90[consistencyData90[,3] %in% pdbList$Representative,]
      
    }
    
    else if(input$seqCluster == "70"){
      
      sequenceCluster = consistencyData70[consistencyData70[,3] %in% pdbList$Representative,]
      
    }
    
    
    else if(input$seqCluster == "40"){
      
      sequenceCluster = consistencyData40[consistencyData40[,3] %in% pdbList$Representative,]
      
    }
    
    
    cols = c("Stoichiometry", "Symmetry")
    
    sequenceCluster$combinedOligomericState <- apply(sequenceCluster[,cols ] , 1 , paste , collapse = "_" )
    
    
    
    sequenceCluster2 = unique(sequenceCluster[,c("PDB ID","Representative","Stoichiometry", "Symmetry", "Consistency score", "combinedOligomericState")])
    
    
    
    
    splitData = split(sequenceCluster2,factor(sequenceCluster2$Representative))
    signatureDataList = list ()
    
    for(i in 1:length(splitData)){
      
      signatureData = splitData[[i]]
      
      for(j in 1:dim(signatureData)[1]){
        
        signatureData$count[j] = sum(sequenceCluster$combinedOligomericState == signatureData$combinedOligomericState[j])
        
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
    
    lastResult = as.data.frame(rbindlist(signatureDataList))
    
    lastResult2 = unique(lastResult[,-1])
    
    sequenceCluster = lastResult2
    
    sequenceCluster$PDBID = currentData()[,1]
    
    sequenceClusterLast = sequenceCluster[,c(7,1:6)]
    names(sequenceClusterLast) = c("PDB ID", "Representative chain","Stoichiometry", "Symmetry",
                                   "Consistency score", "Number of Quaternary Structures", "PDB entries in the cluster")
    
    
    maxCS = sequenceClusterLast$`Consistency score`[which.max(sequenceClusterLast$`Consistency score`)]
    
    sequenceClusterLast$CS = maxCS
    
    if(maxCS > 0.5 && sum(sequenceClusterLast[,6]) > 2){
      
      for(i in 1: dim(sequenceClusterLast)[1]){
        sequenceClusterLast$Res[i] = if(sequenceClusterLast$CS[i] == sequenceClusterLast$`Consistency score`[i]){1}else{0}
      }
    } else{
      
      sequenceClusterLast$Res = 0
    }
    
    resultSeq =  sequenceClusterLast[-8]
    
    resultSeq[,5] = round(resultSeq[,5], 2)
    
    resultSeq
  })
  
  ############################################################################################
  ############### Sequence cluster - End ######################################
  ############################################################################################
  
  
  ############################################################################################
  ############### Signature results - Start ######################################
  ############################################################################################
  
  output$signature<- DT::renderDataTable({
    
    if(!input$signatureResults){
      
      return()
      
    }
    
    if(input$signatureResults && input$startAnalysis){
      
      df = as.data.frame(sequenceCluster())
      
      datatable(df, escape=FALSE,  rownames=FALSE, class = 'cell-border hover stripe', extensions = c('Buttons'), options = list(dom = 'T<"clear">lfrtip',buttons = c('copy', 'excel', 'pdf'),  columnDefs = list(list(targets = 7, visible = FALSE)))
                
                
      )%>% formatStyle(
        'Res',
        target = 'row',
        backgroundColor = styleEqual(c(1), c('#00BA38'))
      )
      
    }
    
  })
  
  
  
  ############################################################################################
  ############### Signature results - End ######################################
  ############################################################################################
  
  
  
  ############################################################################################
  ############### PISA - Start ######################################################
  ############################################################################################
  selectedPisa <- reactive({
    
    pdbId = as.character(currentData()[,1])
    
    pisa = pisa_results
    pisa = unique(pisa)
    
    pisa$pdbId = toupper(pisa$pdbId)
    rownames(pisa) = pisa$pdbId
    names(pisa) = c("PDB ID", "Assessment", "Structure",  "Stoichiometry", "Symmetry")
    
    selectedPisa = pisa[as.character(pdbId),]
    
    if(!is.na(selectedPisa[1,1])){
      
      selectedPisa
    }else{
      
      selectedPisa$`PDB ID` = pdbId
      selectedPisa$Structure = "Inconclusive"
      selectedPisa$Stoichiometry = "Inconclusive"
      selectedPisa$Symmetry = "Inconclusive"
      selectedPisa
    }
    
  })
  
  
  output$pisa <- DT::renderDataTable({
    
    if (!input$startAnalysis){
      return()
      
    }
    
    if (input$startAnalysis && input$pisaPrediction){
      
      datatable(selectedPisa(), escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
        dom = 'T<"clear">lfrtip', buttons = c('copy', 'excel', 'pdf'), deferRender = TRUE
        
      ))
      
      
    }
    
    
  })
  
  
  ############################################################################################
  ############### PISA - End ######################################################
  ############################################################################################
  
  ############################################################################################
  ############### EPPIC - Start ######################################################
  ############################################################################################
  selectedEppic <- reactive({
    
    pdbId = currentData()[,1]
    
    eppic = eppic_results
    #eppic = unique(eppic)
    
    eppic$pdbId = toupper(eppic$pdbId)
    rownames(eppic) = eppic$pdbId
    names(eppic) = c("PDB ID", "Structure",  "Stoichiometry", "Symmetry")
    
    selectedEppic = eppic[as.character(pdbId),]
    
    
    if(!is.na(selectedEppic$`PDB ID`)){
      
      selectedEppic
    }else{
      
      selectedEppic$`PDB ID` = pdbId
      selectedEppic$Structure = "Inconclusive"
      selectedEppic$Stoichiometry = "Inconclusive"
      selectedEppic$Symmetry = "Inconclusive"
      selectedEppic
    }
    
  })
  
  
  output$eppic <- DT::renderDataTable({
    
    if (!input$startAnalysis){
      return()
      
    }
    
    if (input$startAnalysis && input$eppicPrediction){
      
      datatable(selectedEppic(), escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
        dom = 'T<"clear">lfrtip', buttons = c('copy', 'excel', 'pdf'), deferRender = TRUE
        
      ))
      
      
    }
    
    
  })
  
  
  ############################################################################################
  ############### PISA - End ######################################################
  ############################################################################################
  
  
  
  ############################################################################################
  ############### 3D view with jmol - Start ######################################
  ############################################################################################
  
  
  threeDim <- reactive({
    
    pdbID = currentData()[,1]
    
    mol <- HTML(paste('<script src="/js/3Dmol-min.js"></script><div style="height:400px;width:400px;position: absolute;"class="viewer_3Dmoljs" data-href="jmol/pdbFiles/',pdbID,'_pisa.pdb" data-backgroundcolor="0xffffff"  data-style="cartoon:color=spectrum" ></div>', sep=''))
    
    
    
    
    return(mol)
    
    
  })
  
  
  
  
  output$jMolRes<- renderPrint({
    
    if(input$startAnalysis && input$pisaPrediction && input$Display ){
      
      
      return(threeDim())
      
    }
    
    
    #HTML('<script src="http://proteinformatics.charite.de/ngl/js/build/ngl.embedded.min.js"></script><script>if(!Detector.webgl)Detector.addGetWebGLMessage();document.addEventListener("DOMContentLoaded", function(){NGL.init(function(){varstage = new NGL.Stage( "viewport");stage.loadFile("rcsb://1crn");});});</script><div id="viewport" style="width:400px; height:300px;"></div>')
    
  })
  
  ############################################################################################
  ############### 3D view with jmol - End ######################################
  ############################################################################################
  
  
  ############################################################################################
  ############### PISA advanced results - Start ######################################
  ############################################################################################
  
  
  selectedPisaAdvanced <- reactive({
    
    pdbId = currentData()[,1]
    
    pisaAdvanced = pisa_advanced_results
    #pisaAdvanced = unique(pisaAdvanced)
    
    #pisaAdvanced$pdbId = toupper(pisaAdvanced$pdbId)
    rownames(pisaAdvanced) = pisaAdvanced$pdbId
    names(pisaAdvanced) = c("PDB ID", "Assessment", "Structure", "Formula", "Composition", "ASA", "BSA", "Dissociation energy","Entropy", "Dissociation area", "Internal energy", "MM size")
    selectedPisaAdvanced = pisaAdvanced[as.character(pdbId),]
    
    selectedPisaAdvanced
    
  })
  
  output$pisaAdvancedResults <- DT::renderDataTable({
    
    if (!input$startAnalysis){
      return()
      
    }
    
    if (input$startAnalysis && input$pisaPrediction && input$advancedPisa){
      
      datatable(selectedPisaAdvanced(), escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
        dom = 'T<"clear">lfrtip', buttons = c('copy', 'excel', 'pdf'), deferRender = TRUE
        
      ))
      
      
    }
    
    
  })
  ############################################################################################
  ############### 3D view with jmol - End ######################################
  ############################################################################################
  
  
  
  ############################################################################################
  ############### Sequence cluster2 - Start ######################################
  ############################################################################################
  
  sequenceCluster2 <- reactive({
    
    pdbList = outlierDetect()
    
    if(input$seqCluster == "95"){
      
      sequenceCluster = consistencyData95[consistencyData95[,3] %in% pdbList$Representative,]
      
    }
    
    else if(input$seqCluster == "90"){
      
      sequenceCluster = consistencyData90[consistencyData90[,3] %in% pdbList$Representative,]
      
    }
    
    else if(input$seqCluster == "70"){
      
      sequenceCluster = consistencyData70[consistencyData70[,3] %in% pdbList$Representative,]
      
    }
    
    
    else if(input$seqCluster == "40"){
      
      sequenceCluster = consistencyData40[consistencyData40[,3] %in% pdbList$Representative,]
      
    }
    
    
    cols = c("Stoichiometry", "Symmetry")
    
    sequenceCluster$combinedOligomericState <- apply(sequenceCluster[,cols ] , 1 , paste , collapse = "_" )
    
    
    
    sequenceCluster2 = unique(sequenceCluster[,c("PDB ID","Representative","Stoichiometry", "Symmetry", "Consistency score", "combinedOligomericState")])
    
    
    
    
    splitData = split(sequenceCluster2,factor(sequenceCluster2$Representative))
    signatureDataList = list ()
    
    for(i in 1:length(splitData)){
      
      signatureData = splitData[[i]]
      
      for(j in 1:dim(signatureData)[1]){
        
        signatureData$count[j] = sum(sequenceCluster$combinedOligomericState == signatureData$combinedOligomericState[j])
        
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
    
    lastResult = as.data.frame(rbindlist(signatureDataList))
    lastResult$PDBID = NA
    
    if(length(unique(lastResult$Representative)) > 1){
      
      splitReps = split(lastResult, factor(lastResult$Representative))
      
      RepsList = list()
      
      for(j in 1:length(splitReps)){
        
        splitReps[[j]]$PDBID = toupper(splitReps[[j]]$`PDB ID`[which(splitReps[[j]]$`PDB ID` %in% toupper(currentData()[,1]))])[1]
        
        #RepsList[[j]] =splitReps[[j]][which.max(splitReps[[j]][,4]),]
        RepsList[[j]] =splitReps[[j]]
        
      }
      
      seqCluster2 = as.data.frame(rbindlist(RepsList))
      
    }else{
      
      lastResult$PDBID = toupper(lastResult$`PDB ID`[which(lastResult$`PDB ID` %in% toupper(currentData()[,1]))])[1]
      
      seqCluster2 =lastResult
    }
    
    seqRes = unique(seqCluster2[-1])
    
    seqRes2 = seqRes[,c(7,1:6)]
    
    names(seqRes2) = c("PDB ID", "Representative chain", "Stoichiometry", "Symmetry", "Consistency score", "Number of Quaternary Structure", "PDB IDs")
    seqRes2
    
  })
  
  ############################################################################################
  ############### Sequence cluster2 - End ######################################
  ############################################################################################
  
  currentPDB <- reactive({
    
    current = consistencyScore
    c = current[current$pdbId%in%currentData()[,1],]
    c[,c("pdbId", "bioAssemblyNr95", "stoich95", "symmetry95")]
    
    
  })
  
  
  ############################################################################################
  ############### Combined results- Start ######################################
  ############################################################################################
  
  
  combined <- reactive({
    if (input$startAnalysis){
      
      
      current = currentPDB()
      current2 = current[,1:3]
      names(current2) = c("PDB ID", "BA Number", "PDB")
      
      
      if(length(unique(sequenceCluster2()[, "PDB ID"])) > 1){
        
        byPDB = list()
        
        splitByPDB = split(sequenceCluster2(), sequenceCluster2()[,1])
        
        for(i in 1:length(splitByPDB)){
          
          byPDB[[i]] = splitByPDB[[i]][which.max(splitByPDB[[i]][,5]),]
          
          
        }
        
        seqCluster2 = as.data.frame(rbindlist(byPDB))
        
        seqCluster3 = seqCluster2[,c(1,3,5,6)]
        
        seqCluster3$Stoichiometry = as.character(seqCluster3$Stoichiometry)
        
        for (i in 1:dim(seqCluster3)[1]){
          
          if(seqCluster3[i,3] <= 0.5 || seqCluster3[i,4] < 3){
            
            seqCluster3[i,2] = "Inconclusive"
            
          }
          
        }
        
        seqCluster3 = seqCluster3[,1:2]
        
        names(seqCluster3) = c("PDB ID", "Sequence Cluster")
        
      }else{
        
        seqCluster2 = sequenceCluster2()
        seqCluster3 = seqCluster2[,c(1,3,5,6)]
        
        
        if(max(seqCluster3[,3]) <= 0.5 || sum(seqCluster3[,4]) < 3){
          
          seqCluster3[,2] = "Inconclusive"
          
        }
        
        seqCluster4 = seqCluster3[which.max(seqCluster3[,3]),]
        
        
        seqCluster3 = seqCluster4[,1:2]
        
        names(seqCluster3) = c("PDB ID", "Sequence Cluster")
      }
      pisaRes = selectedPisa()
      pisaRes2 = pisaRes[,c(1,4)]
      names(pisaRes2) = c("PDB ID","PISA")
      
      eppicRes = selectedEppic()
      eppicRes2 = eppicRes[,c(1,3)]
      names(eppicRes2) = c("PDB ID","EPPIC")
      
      if(input$OligomericStatePrediction){
      
      tm = textMiner()[[1]]
      tm2 = tm[,c(1:2,4,5)]
      
      if(("TRUE" %in% as.character(tm2[,3])) && (max(tm2[,4]) > 0.51)){
        
        tm2[,1] = substring(tm2[,1], 88, 91)
        tm2 =tm2[,1:2]
      }else{
        
        tm2[,1] = substring(tm2[,1], 88, 91)
        tm2[,2] = "Inconclusive"
        
        tm2 =tm2[,1:2]
        
        
      }
      
      names(tm2) = c("PDB ID", "Text Mining")
      
      }
      
      #seqCluster2 = seqCluster[which.max(seqCluster[,4]),]
      
      
      if(dim(current2)[1] > 1){
        
        
        #merge1 = merge(current2,seqCluster3,by = 'PDB ID')
        #merge2 = merge(pisaRes2,tm2,by = 'PDB ID')
        
        #consensus = merge(merge1, merge2, by = 'PDB ID')
        if(input$FourMethodConsensus){
        consensus = join_all(list(current2, seqCluster3, pisaRes2, eppicRes2, tm2), by = 'PDB ID', type = 'left')
        }else{consensus = join_all(list(current2, seqCluster3, pisaRes2, eppicRes2), by = 'PDB ID', type = 'left')}
        
      }else{
        
        if(input$FourMethodConsensus){
         consensus = cbind(current[,1:3], seqCluster3[,2], pisaRes[,4], eppicRes2[,2], tm2[,2])
        }else{consensus = cbind(current[,1:3], seqCluster3[,2], pisaRes[,4], eppicRes2[,2])}
        
      }
      
      if(input$FourMethodConsensus){
        
        names(consensus) = c("pdbId", "BaNumber", "CurrentStoichiometry", "SequenceClusterStoichiometry", "PISAStoichiometry", "EPPICStoichiometry", "TextMiningStoichiometry")
        
      }else{
        names(consensus) = c("pdbId", "BaNumber", "CurrentStoichiometry", "SequenceClusterStoichiometry", "PISAStoichiometry", "EPPICStoichiometry")
      }
      
      if(consensus$CurrentStoichiometry != "Inconclusive"){
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A"] ="monomer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A2"] ="dimer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "AB"] ="dimer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A3"] ="trimer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A2B"] ="trimer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "ABC"] ="trimer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A4"] ="tetramer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A2B2"] ="tetramer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "ABCD"] ="tetramer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A5"] ="pentamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A6"] ="hexamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A3B3"] ="hexamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A2B2C2"] ="hexamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A7"] ="heptamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A8"] ="octamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A4B4"] ="octamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A2B2C2D2"] ="octamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A9"] ="nonamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A10"] ="decamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A5B5"] ="decamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A11"] ="undecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A12"] ="dodecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A6B6"] ="dodecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A13"] ="tridecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A14"] ="tetradecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A15"] ="pentadecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A16"] ="hexadecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A17"] ="heptadecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A18"] ="octadecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A19"] ="nonadecamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A20"] ="eicosamer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A21"] ="21-mer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A22"] ="22-mer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A23"] ="23-mer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A12B12"] ="24-mer"
        consensus$CurrentOligomericState[consensus$CurrentStoichiometry == "A24"] ="24-mer"
      }else{
        
        consensus$CurrentOligomericState = "Inconclusive"
      }
      
      if(consensus$SequenceClusterStoichiometry != "Inconclusive"){
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A"] ="monomer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A2"] ="dimer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "AB"] ="dimer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A3"] ="trimer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A2B"] ="trimer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "ABC"] ="trimer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A4"] ="tetramer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A2B2"] ="tetramer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "ABCD"] ="tetramer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A5"] ="pentamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A6"] ="hexamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A3B3"] ="hexamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A2B2C2"] ="hexamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A7"] ="heptamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A8"] ="octamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A4B4"] ="octamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A2B2C2D2"] ="octamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A9"] ="nonamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A10"] ="decamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A5B5"] ="decamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A11"] ="undecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A12"] ="dodecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A6B6"] ="dodecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A13"] ="tridecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A14"] ="tetradecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A15"] ="pentadecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A16"] ="hexadecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A17"] ="heptadecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A18"] ="octadecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A19"] ="nonadecamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A20"] ="eicosamer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A21"] ="21-mer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A22"] ="22-mer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A23"] ="23-mer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A12B12"] ="24-mer"
        consensus$SequenceClusterOligomericState[consensus$SequenceClusterStoichiometry == "A24"] ="24-mer"
        
      }else{
        
        consensus$SequenceClusterOligomericState = "Inconclusive"
      }
      
      if(consensus$PISAStoichiometry != "Inconclusive"){
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A"] ="monomer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A2"] ="dimer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "AB"] ="dimer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A3"] ="trimer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A2B"] ="trimer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "ABC"] ="trimer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A4"] ="tetramer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A2B2"] ="tetramer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "ABCD"] ="tetramer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A5"] ="pentamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A6"] ="hexamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A3B3"] ="hexamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A2B2C2"] ="hexamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A7"] ="heptamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A8"] ="octamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A4B4"] ="octamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A2B2C2D2"] ="octamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A9"] ="nonamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A10"] ="decamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A5B5"] ="decamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A11"] ="undecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A12"] ="dodecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A6B6"] ="dodecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A13"] ="tridecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A14"] ="tetradecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A15"] ="pentadecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A16"] ="hexadecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A17"] ="heptadecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A18"] ="octadecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A19"] ="nonadecamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A20"] ="eicosamer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A21"] ="21-mer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A22"] ="22-mer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A23"] ="23-mer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A12B12"] ="24-mer"
        consensus$PISAOligomericState[consensus$PISAStoichiometry == "A24"] ="24-mer"
        
      }else{
        
        consensus$PISAOligomericState = "Inconclusive"
      }
      
      
      if(consensus$EPPICStoichiometry != "Inconclusive"){
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A"] ="monomer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A2"] ="dimer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "AB"] ="dimer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A3"] ="trimer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A2B"] ="trimer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "ABC"] ="trimer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A4"] ="tetramer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A2B2"] ="tetramer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "ABCD"] ="tetramer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A5"] ="pentamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A6"] ="hexamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A3B3"] ="hexamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A2B2C2"] ="hexamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A7"] ="heptamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A8"] ="octamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A4B4"] ="octamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A2B2C2D2"] ="octamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A9"] ="nonamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A10"] ="decamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A5B5"] ="decamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A11"] ="undecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A12"] ="dodecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A6B6"] ="dodecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A13"] ="tridecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A14"] ="tetradecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A15"] ="pentadecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A16"] ="hexadecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A17"] ="heptadecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A18"] ="octadecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A19"] ="nonadecamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A20"] ="eicosamer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A21"] ="21-mer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A22"] ="22-mer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A23"] ="23-mer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A12B12"] ="24-mer"
        consensus$EPPICOligomericState[consensus$EPPICStoichiometry == "A24"] ="24-mer"
        
      }else{
        
        consensus$EPPICOligomericState = "Inconclusive"
      }
      
      
      
      if(input$FourMethodConsensus){ 
        consensusLast = consensus[,c("pdbId", "BaNumber", "CurrentOligomericState", "SequenceClusterOligomericState", "PISAOligomericState", "EPPICOligomericState", "TextMiningStoichiometry")]
      }else{
        consensusLast = consensus[,c("pdbId", "BaNumber", "CurrentOligomericState", "SequenceClusterOligomericState", "PISAOligomericState", "EPPICOligomericState")]
      }
      for(i in 1:dim(consensusLast)[1]){
        if(input$FourMethodConsensus){ 
          vote = c(as.character(consensusLast$SequenceClusterOligomericState[i]), as.character(consensusLast$PISAOligomericState[i]), as.character(consensusLast$EPPICOligomericState[i]), as.character(consensusLast$TextMiningStoichiometry[i]))
        }else{
        vote = c(as.character(consensusLast$SequenceClusterOligomericState[i]), as.character(consensusLast$PISAOligomericState[i]), as.character(consensusLast$EPPICOligomericState[i]))
        }
        v = as.data.frame(table(vote))
        v2 = v[v$Freq>2,]
        
        consensusLast$ConsensusResult[i] = if(max(table(vote)) <3){"Inconclusive"} else if(length(v2$vote)>0){names(which.max(table(vote)))}else if(min(v$Freq)>2){names(which.max(table(vote)))}else{"Inconclusive"}
      }
      
      
      for(i in 1:dim(consensusLast)[1]){
        consensusLast$ResultC[i] = if(as.character(consensusLast$CurrentOligomericState)[i] == as.character(consensusLast$ConsensusResult)[i]){1}else if(as.character(consensusLast$ConsensusResult)[i] == "Inconclusive"){2}else{0}
        
        consensusLast$ResultSC[i] = if(as.character(consensusLast$CurrentOligomericState)[i] == as.character(consensusLast$SequenceClusterOligomericState)[i]){1}else if(as.character(consensusLast$SequenceClusterOligomericState)[i] == "Inconclusive"){2}else{0}
        
        consensusLast$ResultP[i] = if(as.character(consensusLast$CurrentOligomericState)[i] == as.character(consensusLast$PISAOligomericState)[i]){1}else if(as.character(consensusLast$PISAOligomericState)[i] == "Inconclusive"){2}else{0}
        
        consensusLast$ResultE[i] = if(as.character(consensusLast$CurrentOligomericState)[i] == as.character(consensusLast$EPPICOligomericState)[i]){1}else if(as.character(consensusLast$EPPICOligomericState)[i] == "Inconclusive"){2}else{0}
        
        if(input$FourMethodConsensus){ 
         consensusLast$ResultTM[i] = if(as.character(consensusLast$CurrentOligomericState)[i] == as.character(consensusLast$TextMiningStoichiometry)[i]){1}else if(as.character(consensusLast$TextMiningStoichiometry)[i] == "Inconclusive"){2}else{0}
        }
      }
      
      if(input$FourMethodConsensus){
       names(consensusLast) = c("PDB ID", "BA Number", "PDB", "Sequence Cluster",
                                "PISA", "EPPIC", "Consensus", "ResultC", "ResultSC", "ResultP", "ResultE", "ResultTM")
      }else{
       names(consensusLast) = c("PDB ID", "BA Number", "PDB", "Sequence Cluster",
                                "PISA", "EPPIC", "Consensus", "ResultC", "ResultSC", "ResultP", "ResultE")
      }
      consensusLast
      
    }
  })
  
  
  ###############################################################################################
  
  combinedSymmetry <- reactive({
    if (input$startAnalysis){
      
      
      current = currentPDB()
      current2 = current[,c(1,2,4)]
      names(current2) = c("PDB ID", "BA Number", "PDB")
      
      
      if(length(unique(sequenceCluster2()[, "PDB ID"])) > 1){
        
        byPDB = list()
        
        splitByPDB = split(sequenceCluster2(), sequenceCluster2()[,1])
        
        for(i in 1:length(splitByPDB)){
          
          byPDB[[i]] = splitByPDB[[i]][which.max(splitByPDB[[i]][,5]),]
          
          
        }
        
        seqCluster2 = as.data.frame(rbindlist(byPDB))
        
        seqCluster3 = seqCluster2[,c(1,4,5,6)]
        
        seqCluster3$Symmetry = as.character(seqCluster3$Symmetry)
        
        for (i in 1:dim(seqCluster3)[1]){
          
          if(seqCluster3[i,3] <= 0.5 || seqCluster3[i,4] < 3){
            
            seqCluster3[i,2] = "Inconclusive"
            
          }
          
        }
        
        seqCluster3 = seqCluster3[,1:2]
        
        names(seqCluster3) = c("PDB ID", "Sequence Cluster")
        
      }else{
        
        seqCluster2 = sequenceCluster2()
        seqCluster3 = seqCluster2[,c(1,4,5,6)]
        
        
        if(max(seqCluster3[,3]) <= 0.5 || sum(seqCluster3[,4]) < 3){
          
          seqCluster3[,2] = "Inconclusive"
          
        }
        
        seqCluster4 = seqCluster3[which.max(seqCluster3[,3]),]
        seqCluster3 = seqCluster4[,1:2]
        names(seqCluster3) = c("PDB ID", "Sequence Cluster")
      }
      
      pisaRes = selectedPisa()
      pisaRes2 = pisaRes[,c(1,5)]
      names(pisaRes2) = c("PDB ID","PISA")
      
      eppicRes = selectedEppic()
      eppicRes2 = eppicRes[,c(1,4)]
      names(eppicRes2) = c("PDB ID","EPPIC")
      
      
      if(dim(current2)[1] > 1){
        
        consensusLast = join_all(list(current2,seqCluster3,pisaRes2,eppicRes2), by = 'PDB ID', type = 'left')
        
      }else{
        
        consensusLast = cbind(current[,c(1,2,4)], seqCluster3[,2], pisaRes[,5], eppicRes2[2])
        
      }
      
      names(consensusLast) = c("pdbId", "BaNumber", "CurrentSymmetry", "SequenceClusterSymmetry", "PISASymmetry", "EPPICSymmetry")
      
      for(i in 1:dim(consensusLast)[1]){
        vote = c(as.character(consensusLast$SequenceClusterSymmetry[i]), as.character(consensusLast$PISASymmetry[i]),as.character(consensusLast$EPPICSymmetry[i]))
        
        consensusLast$ConsensusResult[i] = if(max(table(vote)) == 1){"Inconclusive"}else{names(which.max(table(vote)))}
      }
      
      for(i in 1:dim(consensusLast)[1]){
        consensusLast$ResultC[i] = if(as.character(consensusLast$CurrentSymmetry)[i] == as.character(consensusLast$ConsensusResult)[i]){1}else if(as.character(consensusLast$ConsensusResult)[i] == "Inconclusive"){2}else{0}
        
        consensusLast$ResultSC[i] = if(as.character(consensusLast$CurrentSymmetry)[i] == as.character(consensusLast$SequenceClusterSymmetry)[i]){1}else if(as.character(consensusLast$SequenceClusterSymmetry)[i] == "Inconclusive"){2}else{0}
        
        consensusLast$ResultP[i] = if(as.character(consensusLast$CurrentSymmetry)[i] == as.character(consensusLast$PISASymmetry)[i]){1}else if(as.character(consensusLast$PISASymmetry)[i] == "Inconclusive"){2}else{0}
        
        consensusLast$ResultE[i] = if(as.character(consensusLast$CurrentSymmetry)[i] == as.character(consensusLast$EPPICSymmetry)[i]){1}else if(as.character(consensusLast$EPPICSymmetry)[i] == "Inconclusive"){2}else{0}
        
      }
      
      names(consensusLast) = c("PDB ID", "BA Number", "PDB", "Sequence Cluster",
                               "PISA", "EPPIC", "Consensus", "ResultC", "ResultSC", "ResultP", "ResultE")
      
      consensusLast
      
    }
    
  })
  
  
  #############################################################################################
  
  
  
  output$combinedResults <- DT::renderDataTable({
    
    
    df = as.data.frame(combined())
    
    #datatable(currentData(), ,escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
    #dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf')#, colVis = list(activate = 'click', showAll = TRUE)
    
    
    datatable(df, escape=FALSE,  rownames=FALSE, class = 'cell-border hover stripe', extensions = c('Buttons'), options = list(dom = 'T<"clear">lfrtip',buttons = c('copy', 'excel', 'pdf'),  columnDefs = list(list(targets = c(7:10), visible = FALSE)))
              
    )%>% formatStyle(
      columns = c("Consensus", "Sequence Cluster", "PISA", "EPPIC"), valueColumns = c("ResultC","ResultSC", "ResultP", "ResultE"),
      backgroundColor = styleEqual(c(0, 1, 2), c('#F8766D', '#00BA38', '#619CFF'))
    )
    
    
  })
  # columns = c("Consensus", "Sequence Cluster", "PISA", "EPPIC", "Text Mining"), valueColumns = c("ResultC","ResultSC", "ResultP", "ResultE", "ResultTM"),
  
  output$combinedSymmetryResults <- DT::renderDataTable({
    
    df = as.data.frame(combinedSymmetry())
    
    #datatable(currentData(), ,escape=FALSE, rownames=FALSE,  class = 'cell-border hover stripe', extensions = c('Buttons', 'Responsive'), options = list(
    #dom = 'T<"clear">lfrtip',  buttons = c('copy', 'excel', 'pdf')#, colVis = list(activate = 'click', showAll = TRUE)
    
    
    datatable(df, escape=FALSE,  rownames=FALSE, class = 'cell-border hover stripe', extensions = c('Buttons'), options = list(dom = 'T<"clear">lfrtip',buttons = c('copy', 'excel', 'pdf'),  columnDefs = list(list(targets = c(7:10), visible = FALSE)))
              
    )%>% formatStyle(
      columns = c("Consensus", "Sequence Cluster", "PISA", "EPPIC"), valueColumns = c("ResultC","ResultSC", "ResultP", "ResultE"),
      backgroundColor = styleEqual(c(0, 1, 2), c('#F8766D', '#00BA38', '#619CFF'))
    )
    
    
    
    
  })
  
  ############################################################################################
  ############### Combined results - End ######################################
  ############################################################################################
  #output$path <- renderDataTable({
  #    
  #    if(!input$startAnalysis){
  #
  #        invisible()
  #
  #    }
  #
  #else if(is.na(currentData()[,1])){
  #    
  #    Error = "Invalid PDB ID"
  #    
  #    as.data.frame(Error)
  #    
  #}else{
  #    
  #    invisible()
  
  #}
  
  
  
  #},rownames=FALSE, options = list(iDisplayLength = 1,bSearchable = FALSE
  #,bFilter=FALSE,bPaginate=FALSE,bAutoWidth=TRUE
  #,bInfo=0,bSort=0
  #, "sDom" = "rt"
  #))
  })