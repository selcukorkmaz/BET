library("shinythemes")
library("shinyBS")

shinyUI(
fluidPage(
theme = "css/mytheme.css",

#titlePanel("Biological Assembly Prediction Tool"),


sidebarLayout(
	sidebarPanel(width=3,


conditionalPanel("input.tabs=='Evaluation results'",
            h5(HTML('<p><b>Enter a PDB ID or upload a file</b></p>'),

#br(),

#checkboxInput(inputId = "pmc", label = "PMC Full Text Mining", value = FALSE),
#conditionalPanel(condition = "input.pmc",

radioButtons(inputId = "dataInput", "", list("PDB ID"=1, "Upload file"=2, "Upload PDF"=3), selected=1),
#bsTooltip(id = "input.dataInput=='1'", title="Single PDB entry", placement = "bottom", trigger = "focus",options = NULL),
#bsTooltip(id = "dataInput", title="Use the first option for the single PDB entry. Use the second option for the multiple PDB entries. Use the third option for mining a PDF paper.", placement = "bottom", trigger = "hover", options = NULL),


        conditionalPanel(condition="input.dataInput=='1'",
br(),


 textInput(inputId = "PdbId", label = "Enter a PDB ID", value = "1Z77"),
bsTooltip(id = "PdbId", title="Enter a 4 characters long ID", placement = "bottom", trigger = "hover", options = NULL)


#        selectizeInput("MLalgorithm", "Select an ML algorithm", choices = c("Support vector machines", "Boosted logistic regression", "Random forest"), multiple = FALSE, selected = "Boosted logistic regression")
        ),

        conditionalPanel(condition="input.dataInput=='2'",
# HTML('<br>'),
            fileInput("upload", "", multiple = FALSE, accept = "text/plain"),
bsTooltip(id = "upload", title="Upload a .txt file which includes multiple PDB entries. First row must be header.", placement = "bottom", trigger = "hover", options = NULL)
# HTML('<p><b>NOTE: First row must be header.</b></p>')
        ),

conditionalPanel(condition="input.dataInput=='3'",
    fileInput("uploadPDF", "Upload a paper", multiple = FALSE, accept = "pdf"),
    bsTooltip(id = "uploadPDF", title="Upload a PDF paper for a PDB entry.", placement = "bottom", trigger = "hover", options = NULL),
    textInput(inputId = "PdbIdPDF", label = "Enter a PDB ID"),
    bsTooltip(id = "PdbIdPDF", title="Enter a PDB entry for the uploaded PDF paper.", placement = "bottom", trigger = "hover", options = NULL)

),

br(),


actionButton(inputId = "startAnalysis", label = "Evaluate", icon = NULL),
br(),
br(),

#h6("Check the box to include data mentions"),
#checkboxInput("dataMentions", HTML('<p><font size="2">Include data mentions</font></p>'), FALSE),

h5("Quaternary structure evaluation methods"),


fluidRow(column(3, checkboxInput(inputId = "signatureResults", label = "SC", value = FALSE)),
column(3, checkboxInput(inputId = "OligomericStatePrediction", label = "TM", value = FALSE)),
column(3, checkboxInput(inputId = "pisaPrediction", label = "PISA", value = FALSE)),
column(3, checkboxInput("eppicPrediction", "EPPIC", FALSE))
),

bsTooltip(id = "signatureResults", title="Sequence cluster", placement = "bottom", trigger = "hover", options = NULL),

bsTooltip(id = "OligomericStatePrediction", title="Text mining", placement = "bottom", trigger = "hover", options = NULL)







),

h5("Advanced results and options"),


fluidRow(column(6, checkboxInput(inputId = "advancedResults", label = "Results", value = FALSE)),
column(3),
column(6, checkboxInput(inputId = "advancedOptions", label = "Options", value = FALSE))
),


#checkboxInput(inputId = "advancedResults", label = "Advanced results", value = FALSE),

    conditionalPanel(condition = "input.advancedResults",
        checkboxInput(inputId = "mlFiltering", label = "Sentences", value = FALSE),
        checkboxInput(inputId = "advancedPisa", label = "Pisa advanced results", value = FALSE),
        checkboxInput(inputId = "Display", label = "Display 3D structure", value = FALSE)

),

#conditionalPanel(condition = "input.AdvancedOutlierDetection",

#    checkboxInput(inputId = "outlierDetectionPlot", label = "Sequence cluster plot", value = FALSE)

#),


#conditionalPanel(condition = "input.AdvancedOligomericStatePrediction",


#fluidRow(column(6, checkboxInput("selectAll", "Select All", FALSE)),
#column(1),
#column(6, checkboxInput("deselectAll", "Deselect All", FALSE))),

#h5("Oligomeric state"),
#checkboxInput(inputId = "oliomericPredPlot", label = "Oligomeric state plot", value = FALSE),
#h5("Sentences"),
#checkboxInput(inputId = "mlFilteredSen", label = "ML filtered sentences", value = FALSE),
#checkboxInput(inputId = "mlFiltering", label = "All sentences", value = FALSE),
#h5("Other"),
#checkboxInput(inputId = "mutantResults", label = "Mutant protein information", value = FALSE),
#checkboxInput(inputId = "noOligomericResults", label = "Papers with no oligomeric results", value = FALSE),
#checkboxInput(inputId = "noFullTextResults", label = "Papers with no full text/abstract results", value = FALSE),
#checkboxInput(inputId = "publicationInfo", label = "Publication information", value = FALSE),
#checkboxInput(inputId = "summaryResults", label = "Summary", value = FALSE)
#),




#conditionalPanel(condition = "input.AdvancedPisaPrediction",
#checkboxInput(inputId = "advancedPisa", label = "Pisa advanced results", value = FALSE),

#checkboxInput(inputId = "Display", label = "Display 3D structure", value = FALSE)


#),


conditionalPanel(condition = "input.AdvancedEppic"
),



#checkboxInput(inputId = "advancedResults", label = "Advanced results", value = FALSE),


conditionalPanel(condition = "input.advancedOptions",


selectizeInput("MLalgorithm", "Select an ML algorithm", choices = c("Support vector machines", "Boosted logistic regression"), multiple = FALSE, selected = "Boosted logistic regression"),


selectizeInput("seqCluster", "Select a sequence cluster", choices = c("95% sequence cluster" = "95","90% sequence cluster" = "90", "70% sequence cluster" = "70", "40% sequence cluster" = "40"), multiple = FALSE, selected = "70")

),


br()

#checkboxInput(inputId = "advancedOptions", label = "Advanced options", value = FALSE),

),

conditionalPanel("input.tabs=='Help'"

)

    ),


mainPanel(

navbarPage("Quaternary Structure Evaluation Tool", id="tabs", inverse = TRUE, collapsible = TRUE, fluid = TRUE, position = "fixed-top", #class("navbar navbar-inverse"),

        tabPanel("Evaluation results",


# tabsetPanel(
                tabPanel(


h4(textOutput(outputId = "section1")),
            bsAlert("help"),




h4(textOutput(outputId = "section16")),
#        DT::dataTableOutput('path'),

       DT::dataTableOutput('combinedResults'),
#br(),
h4(textOutput(outputId = "section17")),

       DT::dataTableOutput('combinedSymmetryResults'),





#h4(textOutput(outputId = "section4")),
#tags$head(tags$style("#consistency  {white-space: nowrap;  }")),
       DT::dataTableOutput('consistency'),

#h4(textOutput(outputId = "section5")),
#                           plotOutput('consistencyPlot', height = "auto"),

h4(textOutput(outputId = "section11")),
                            DT::dataTableOutput('signature'),


h4(textOutput(outputId = "section1")),
#tags$head(tags$style("#TextMining  {white-space: nowrap;  }")),

                           DT::dataTableOutput('TextMining'),

#h4(textOutput(outputId = "section2")),
#                            plotOutput('ospPlot', height = "auto"),


h4(textOutput(outputId = "section13")),
                            DT::dataTableOutput('pisa'),


h4(textOutput(outputId = "section7")),
                            DT::dataTableOutput('eppic'),


#h4(textOutput(outputId = "section3")),
#                            DT::dataTableOutput('mlFilteredSentences'),

h4(textOutput(outputId = "section6")),
                            DT::dataTableOutput('machineLearning'),

#h4(textOutput(outputId = "section12")),
#                            DT::dataTableOutput('mutantResults'),

#h4(textOutput(outputId = "section7")),
#                            DT::dataTableOutput('noOligomericResult'),

#h4(textOutput(outputId = "section8")),
#                            DT::dataTableOutput('noFullText'),

#h4(textOutput(outputId = "section10")),
#                            DT::dataTableOutput('RawData'),

#h4(textOutput(outputId = "section9")),
#                            DT::dataTableOutput('summary'),



h4(textOutput(outputId = "section15")),
                            DT::dataTableOutput('pisaAdvancedResults'),

h4(textOutput(outputId = "section14")),
                            uiOutput('jMolRes')

            )

#        )
),

        tabPanel("Help",

            tabsetPanel(
                tabPanel('About',
                h4("Introduction", id = "intro"),
tags$p(align = "justify", a("The Protein Data Bank", href = "http://www.rcsb.org/pdb/home/home.do", target="_blank"), " (PDB) provides detailed information about the three-dimensional 
       (3D) structures of biological macromolecules, including proteins and nucleic acids. The PDB was founded in 1971 with only 7 structures and it contains more than 
       118,000 structures as of May 2016. The majority of the 3D structures of quaternary structures in the PDB are determined by X-ray crystallography. However, it is not 
       possible to distinguish biologically relevant contacts from crystal contacts in a crystal lattice by using crystallography alone. Therefore, further experiments, 
       such as gel filtration, size exclusion chromatography, analytical ultracentrifugation, etc., are needed to assign the correct oligomeric state of a quaternary structure. 
       Moreover, the correct quaternary structure information may be inferred by comparison with similar proteins through homology modeling and be provided by the structure 
       depositor as metadata. Furthermore, it can be predicted using analytical methods, such as", a("PISA", href = "http://www.ebi.ac.uk/pdbe/pisa/", target="_blank"), "(Proteins, Interfaces, Structures and Assemblies) 
       (Krissinel and Henrick, 2007) and", a("EPPIC", href = "http://www.eppic-web.org/ewui/", target="_blank"), "(Evolutionary Protein-Protein Interface Classifier) (Duarte et al. 2012)"),

tags$p(align = "justify", "The PDB grows by more than 10,000 structures annually thanks to the researchers all around the world. 
       However, because of the incomplete or unclear experimental data, as well as the errors in the data upload 
       processes, the quaternary structure annotations in the PDB are not always correct and reliable (Capitani et al. 2015). 
       In spite the great efforts to reduce the number of erroneous structures, it is reported that there are 
       significant number of incorrect quaternary structures in the PDB. The error rate is 14% according to 
       Levy (2007), while more recently Baskaran et al. (2014) reported a lower bound of error rate as 7%. 
       Therefore, an extensive study is needed to detect the incorrect structures throughout the archive and to 
       assign the most likely quaternary structures for the possibly incorrect structures."),

tags$p(align = "justify", "Here we developed a web-based tool to detect incorrect biological assembly predictions throughout the PDB repository and to assign the most
       probable biological assemblies for the detected incorrect structures by using four different approaches:"),


tags$ul(
  tags$li("Determine a representative quaternary structure for a given certain sequence identity threshold using sequence cluster approach and consistency score calculation."),
  tags$li("Extract correct oligomeric state information alongside with the experimental evidence from a primary paper for a given PDB entry using a text mining approach."),
  tags$li('Rebuild the quaternary structure using PISA software and predict stoichiometry and symmetry of a protein structure using BioJava. '),
  tags$li('Predict stoichiometry and symmetry of a protein structure using EPPIC software.')
),

tags$p(align = "justify", "Finally, we aggregated results from these four different methods to provide a consensus result in order to predict the most likely quaternary structure of each biological macromolecule in the PDB."),

tags$p("All source code is available at" , a("Github.", href = "https://github.com/selcukorkmaz/BET", target="_blank"), "Please see", a("Manual", href="#manual")," for more detailed information."),
br(),
h4("General workflow of the tool"),
br(),
HTML('<center><img src="screenShots/workFlow.jpg" width = "50%"></center>'),
br(),
br(),

h4("Workflow of the sequence cluster (SC) approach"),
br(),
tags$img(src = "screenShots/sequenceCluster.jpg", width = "100%"),
br(),
br(),

h4("Workflow of the text mining (TM) approach"),
br(),
tags$img(src = "screenShots/MLprocedure.jpg", width = "100%"),
br(),
br(),

h4("Workflow of the PISA approach"),
br(),
tags$img(src = "screenShots/pisaProcedure.jpg", width = "100%"),
br(),
br(),

h4("Consensus approach"),
br(),
tags$img(src = "screenShots/consensus.jpg", width = "100%")


),
tabPanel('Manual', id = "manual",



h4("Usage of the tool"),

tags$ol(
    tags$li(align = "justify", "The tool has a very simple user interface, which only requires PDB entries as an input. Users can either enter a single PDB entry into the box in the tool or 
            upload a .txt file, which includes multiple PDB IDs, using the option in the tool for multiple PDB entries. Furthermore, users can upload a PDF version of a paper
            as a third option for the text mining method. After entering or uploading the proper input to the tool, one can click", tags$b("Evaluate"),"button to start evaluation process for PDB entries."),

br(),

    tags$img(src = "screenShots/homePage.jpg", width = "100%"),


br(),
br(),


    tags$li(align = "justify", "Initially, the tool gives results for both oligomeric state and symmetry information. The first table contains oligomeric state results from four different 
evaluation methods, including sequence clustering, text mining, PISA and EPPIC, and the second table includes results from sequence clustering, PISA and EPPIC. 
Finally, the results are aggregated and majority rule is applied to obtain a consensus result. Each result from different methods is compared to current result in the PDB. 
If the two results do not match, red color is used to indicate the inconsistency between the results. On the other hand, green color is used to show the consistency between 
evaluation methods and the current PDB results, and blue color is used to demonstrate the inconclusive results."),
br(),

tags$img(src = "screenShots/resultPage.jpg", width = "100%"),

br(),
br(),

    tags$li(align = "justify", "The evaluation methods can be selected individually to see the more detailed results for sequence clustering, text mining, PISA and EPPIC. 
            Sequence cluster table gives stoichiometry and symmetry information in the cluster alongside with the consistency score, the number of quaternary structures in 
            the cluster and the PDB entries in the cluster. If the maximum consistency score in the cluster is greater than 0.5 and there are enough biological
            assemblies in the cluster, a green background color will be used to indicate the row that has the highest consistency score and the representative biological 
            assembly for the cluster.  Text mining results table contains the oligomeric state prediction keyword, experimental evidence information, ML (machine-learning) 
            filtering result as TRUE or FALSE, majority probability, publication ID as PMC or PubMed, publication type as full-text or abstract and primary citation 
            indication as yes or no. PISA predictions table gives PDB ID, assessment as stable, grey or unstable, structure as BA (biological assembly) or AU 
            (asymmetric unit) and information that define oligomeric state: stoichiometry and symmetry. Likewise, EPPIC prediction table includes stoichiometry and 
            symmetry information for a particular PDB entry."),
br(),
    tags$img(src = "screenShots/detailedResultPage.jpg", width = "100%"),
br(),
br(),

tags$li(align = "justify", "One can display all extracted sentences from the primary publication alongside with the keywords, experimental evidence information and ML filtering label. 
        The sentences classified as positive by the ML algorithm will be displayed with a green background color, while negative sentences will be displayed with a red 
        background color."),
br(),
tags$img(src = "screenShots/sentences.jpg", width = "100%"),
br(),
br(),

tags$li("PISA advanced results contain a table, which include accessible surface area, buried surface area, dissociation energy, entropy, dissociation area and internal 
        energy, and a 3D visualization tool for the protein structure."),
br(),
tags$img(src = "screenShots/pisa.jpg", width = "100%")

)
                ),




                tabPanel('Model performances',
                    h4("Data set"),
                    HTML('<ul>
                            <li>The dataset used in this study is collected from PDB primary citations, which contained 5500 positive (BA relevant) sentences  and 5500 negative (BA unrelevant) senteces. </li>
                            <li>The dataset is splitted as 80% training set and 20% testing set.</li>
                        </ul>'),
                br(),

                    h4("Model building"),
                    HTML('<ul>
                            <li>We used three machine learning algorithms in this tool: Support Vector Machines (SVM), Boosted Logistic Regression (BLR), Random Forest (RF). </li>
                            <li>We made a grid search and used 10-fold cross-validation to select optimal tuning parameters in the training set. </li>
                            <li>We tested model performance on a new test set. </li>
                        </ul>'),
                br(),

                    h4("Test set performaces"),


                    verbatimTextOutput("testPerformance"),

                    br(),

                    h4("ROC curves"),

                    tags$img(src = "screenShots/ROCplot.pdf", width = "100%")


                    ),



            tabPanel('Keywords',
                     
                h5("These keywords will be used to search and extract quaternary structure related sentences throughout the paper(s)."),
                    verbatimTextOutput("baKeywords"),
                
                h5("These keywords will be used to search and extract experimental evidence."),
                  verbatimTextOutput("eeKeywords")
        
            ),

            tabPanel('News',HTML('<p><b>This page will be available soon...</b></p>')),
            tabPanel('Authors',
                     br(),
                     HTML('<p><a href="http://yunus.hacettepe.edu.tr/~selcuk.korkmaz/" target="_blank"> <b>Selcuk Korkmaz, PhD</b></a><p>'),
                     HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
                     HTML('<p><a href="mailto:selcukorkmaz@gmail.com" target="_blank">selcukorkmaz@gmail.com</a><p>'),
                     
                     HTML('<p><a href="https://www.linkedin.com/in/peterrose" target="_blank"> <b>Peter Rose, PhD</b></a><p>'),
                     HTML('<p>Site Head <a href="http://www.rcsb.org/pdb/home/home.do" target="_blank">RCSB Protein Data Bank West</a> UCSD<p>'),
                     HTML('<p><a href="mailto:peter.rose@rcsb.org" target="_blank">peter.rose@rcsb.org</a><p>'),
                     
                     HTML('<p> <b>Jose Duarte, PhD</b><p>'),
                     HTML('<p>Senior Scientist <a href="http://www.rcsb.org/pdb/home/home.do" target="_blank">RCSB Protein Data Bank West</a> UCSD<p>'),
                     HTML('<p><a href="mailto:jose.duarte@rcsb.org" target="_blank">jose.duarte@rcsb.org</a><p>'),
                     
                     HTML('<p><a href="http://www.spice-3d.org"> <b>Andreas Prlić, PhD</b></a><p>'),
                     HTML('<p>Technical and Scientific Team Lead <a href="http://www.rcsb.org/pdb/home/home.do" target="_blank">RCSB Protein Data Bank West</a> UCSD<p>'),
                     HTML('<p><a href="mailto:andreas.prlic@rcsb.org" target="_blank">andreas.prlic@rcsb.org</a><p>'),
                     
                     HTML('<p><a href="http://yunus.hacettepe.edu.tr/~dincer.goksuluk" target="_blank"> <b>Dincer Goksuluk</b></a><p>'),
                     HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
                     HTML('<p><a href="mailto:dincer.goksuluk@gmail.com" target="_blank">dincer.goksuluk@gmail.com</a><p>'),
     
                     HTML('<p><a href="http://gokmenzararsiz.simplesite.com" target="_blank"> <b>Gokmen Zararsiz, PhD</b></a><p>'),
                     HTML('<p>Erciyes University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
                     HTML('<p><a href="mailto:gokmenzararsiz@hotmail.com" target="_blank">gokmenzararsiz@hotmail.com</a><p>')  
                     
                     
            ),
            tabPanel('Citation',HTML('<p><b>This page will be available soon...</b></p>'))
        ))),

        tags$head(
            tags$style("body {padding-top: 53px};")

        ),


            HTML('<a href="#" class="back-to-top">Back to Top</a>'),

HTML('<script src="js/jquery.js"></script>'),

tags$style(type="text/css",
".shiny-output-error { visibility: hidden; }",
".shiny-output-error:before { visibility: hidden; }"
)




#HTML('<script src="JSmol.min.js"></script>')



#HTML('<link rel="stylesheet" href="css/foundation.css">')



#HTML('<script>jmolInitialize(".")</script>')





#HTML('<applet name="Test_Protein" code="JmolApplet" archive="jmol/JmolAppletSigned.jar" width="95%" height="80%" mayscript="true" > <param name="progressbar" value="true" /#><param name="progresscolor" value="blue" /><param name="boxmessage" value="loading ..." /><param name="emulate" value="chime" /><param name="bgcolor" value="#000000" /><param name="load" #value="4a5a_pisa.pdb" /><param name="script" value=";reset;select all ;zoom 110;wireframe 70;move 0 -360 0 0 0 0 0 0 3 60;" /></applet>')


    )
)))










