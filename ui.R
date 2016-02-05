library("shinythemes")
library("shinyBS")

shinyUI(
fluidPage(
theme = "css/mytheme.css",

#titlePanel("Biological Assembly Prediction Tool"),


sidebarLayout(
	sidebarPanel(width=3,


conditionalPanel("input.tabs=='Analysis'",
            h5(HTML('<p><b>Enter a PDB ID or upload a file</b></p>'),

#br(),

#checkboxInput(inputId = "pmc", label = "PMC Full Text Mining", value = FALSE),
#conditionalPanel(condition = "input.pmc",

radioButtons("dataInput", "", list("PDB ID"=1, "Upload file"=2, "Upload PDF"=3), selected=1),



        conditionalPanel(condition="input.dataInput=='1'",
br(),


 textInput(inputId = "PdbId", label = "Enter a PDB ID", value = "2PMF")


#        selectizeInput("MLalgorithm", "Select an ML algorithm", choices = c("Support vector machines", "Boosted logistic regression", "Random forest"), multiple = FALSE, selected = "Boosted logistic regression")
        ),

        conditionalPanel(condition="input.dataInput=='2'",
# HTML('<br>'),
            fileInput("upload", "", multiple = FALSE),
            HTML('<p><b>NOTE: First row must be header.</b></p>')
        ),

conditionalPanel(condition="input.dataInput=='3'",
    fileInput("uploadPDF", "", multiple = FALSE)
),

br(),


actionButton(inputId = "startAnalysis", label = "Start analysis", icon = NULL),
br(),
br(),

#h6("Check the box to include data mentions"),
#checkboxInput("dataMentions", HTML('<p><font size="2">Include data mentions</font></p>'), FALSE),

h5("Biological assembly evaluation methods"),


fluidRow(column(4, checkboxInput(inputId = "signatureResults", label = "Sequence cluster", value = FALSE)),
column(4, checkboxInput(inputId = "OligomericStatePrediction", label = "Text mining", value = FALSE)),
column(3, checkboxInput(inputId = "pisaPrediction", label = "PISA", value = FALSE))
#column(3, checkboxInput("eppic", "EPPIC", FALSE))
)





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


selectizeInput("MLalgorithm", "Select an ML algorithm", choices = c("Support vector machines", "Boosted logistic regression"), multiple = FALSE, selected = "Support vector machines"),


selectizeInput("seqCluster", "Select a sequence cluster", choices = c("95% sequence cluster" = "95","90% sequence cluster" = "90", "70% sequence cluster" = "70", "40% sequence cluster" = "40"), multiple = FALSE, selected = "70")

),


br()

#checkboxInput(inputId = "advancedOptions", label = "Advanced options", value = FALSE),

),

conditionalPanel("input.tabs=='Help'"

)

    ),


mainPanel(

navbarPage("Biological Assembly Evaluation Tool v.0.2.2", id="tabs", inverse = TRUE, collapsible = TRUE, fluid = TRUE, position = "fixed-top", #class("navbar navbar-inverse"),

        tabPanel("Analysis",


# tabsetPanel(
                tabPanel(


h4(textOutput(outputId = "section1")),
            bsAlert("help"),

h4(textOutput(outputId = "section16")),
       dataTableOutput('combinedResults'),


#h4(textOutput(outputId = "section4")),
#tags$head(tags$style("#consistency  {white-space: nowrap;  }")),
#        dataTableOutput('consistency'),

#h4(textOutput(outputId = "section5")),
#                           plotOutput('consistencyPlot', height = "auto"),

h4(textOutput(outputId = "section11")),
                            dataTableOutput('signature'),


h4(textOutput(outputId = "section1")),
#tags$head(tags$style("#TextMining  {white-space: nowrap;  }")),

                           dataTableOutput('TextMining'),

#h4(textOutput(outputId = "section2")),
#                            plotOutput('ospPlot', height = "auto"),


h4(textOutput(outputId = "section13")),
                            dataTableOutput('pisa'),


#h4(textOutput(outputId = "section3")),
#                            dataTableOutput('mlFilteredSentences'),

h4(textOutput(outputId = "section6")),
                            dataTableOutput('machineLearning'),

#h4(textOutput(outputId = "section12")),
#                            dataTableOutput('mutantResults'),

#h4(textOutput(outputId = "section7")),
#                            dataTableOutput('noOligomericResult'),

#h4(textOutput(outputId = "section8")),
#                            dataTableOutput('noFullText'),

#h4(textOutput(outputId = "section10")),
#                            dataTableOutput('RawData'),

#h4(textOutput(outputId = "section9")),
#                            dataTableOutput('summary'),



h4(textOutput(outputId = "section15")),
                            dataTableOutput('pisaAdvancedResults'),

h4(textOutput(outputId = "section14")),
                            uiOutput('jMolRes')

            )

#        )
),

        tabPanel("Help",

            tabsetPanel(
                tabPanel('About',
                h4("Introduction"),
tags$p("", a("Protein Data Bank", href = "http://www.rcsb.org/pdb/home/home.do", target="_blank"), " (PDB) is a repository, which consists of detailed information about the three-dimensional (3D) structures of macromolecules, such as proteins and nucleic acids. As of August 2015, there are more than 110,000 biological macromolecular structures (i.e. protein, nucleic acid) in the PDB. These structure predictions are either from author deposition, using experimental methods, such as X-ray crystallography, nuclear magnetic resonance (NMR), electron microscopy etc., or PISA (Proteins, Interfaces, Structures and Assemblies) predictions, which is a computational method based on chemical thermodynamics."),

tags$p("Despite the enormous efforts, it is estimated that nearly 15% of the biological assembly (BA) structures in the PDB archive is incorrect. The problem arises beacause of the multiple predictions (either author or PISA) of the same PDB entry. For example, there are 13 BA predictions for PDB ID:",a("1BEN", href = "http://www.rcsb.org/pdb/explore/explore.do?structureId=1ben", target="_blank"),"(see the screenshot below)."),
br(),
tags$img(src = "screenShots/1ben.jpg", width = "100%"),
br(),
br(),
tags$p("If we look at the detailed BA information, we can see that there are some inconsistencies in stoichiometry and symmetry between different BA predictions. For example, 1BEN is a hetero-dimer (AB) and has a cyclic (C1) symmetry according to the BA1 prediction, whereas it is a dodecamer (A6B6) and has a dihedral (D3) symmetry according to the BA3 prediction."),
br(),

tags$img(src = "screenShots/mp.pdf", width = "100%"),

br(),
br(),

h4("Inconsistency in Stoichiometry and Symmetry"),

tags$p("To overcome this problem and to identify incorrect BA predictions (outliers) in the PDB, we can use stoichiometry and symmetry information for each PDB entry and try to detect inconsistencies between multiple BA predictions. For this purpose, we make use of the homologous protein chain sequences. Since a great deal of protein chains in the PDB are similar in terms of sequence similarity,  we can cluster these chains and select only one representative from each group of similar chains. First, we find sequence clusters for 95%, 90%, 70% and 40% sequence similarity. Then, we calculate consistency probability for each PDB entry in each sequence cluster using following formula:"),

tags$p(align="center", "P(Consistencyi) = P(Stoichiometryi)xP(Symmetryi|Stoichiometryi)"),

tags$p("Finally, if consistency probability < 0.5 for any PDB entry, then we declare it as possible outlier."),

h5("Example"),

tags$p("95% Sequence cluster and representative chain is 1W28.A"),

tags$ul(
tags$li("A3 stoichiometry and C3 symmetry: 13 instances"),
tags$li("A stoichiometry and C symmetry: 5 instances"),
tags$li('P(Consistency = A3_C3) = P(A3)xP(C3|A3) = (13/18)x(13/13) = 0.722'),
tags$li('P(Consistency = A1_C1) = P(A1)xP(C1|A1) = (5/18)x(5/5) = 0.278'),
tags$li("Since consistency probability for A1 stoichiometry and C1 symmetry is less than our treshold, 0.5, we can conclude that these biological assembly predictions could be incorrect and further investigations might be neccessary to reveal the correct biological asseblies. In our example, we detect PDB entries",  a("1E5I", href = "http://www.rcsb.org/pdb/explore/explore.do?structureId=1e5i", target="_blank"),",", a("1W28", href = "http://www.rcsb.org/pdb/explore/explore.do?structureId=1w28", target="_blank"), ",", a("1W2A", href = "http://www.rcsb.org/pdb/explore/explore.do?structureId=1w2a", target="_blank"),",", a("1W2N", href = "http://www.rcsb.org/pdb/explore/explore.do?structureId=1w2n", target="_blank"), "and", a("1W2O", href = "http://www.rcsb.org/pdb/explore/explore.do?structureId=1w2o", target="_blank"), "as possible outliers."),
br(),
tags$img(src = "screenShots/outlierDetection.jpg", width = "100%"),
br(),
br()
),

tags$p("In order to predict the most probable BA for these kind of inconsistent structures throughout the PDB archive, we developed a text mining approach. This tool mines the primary citations of the protein structures in order to extract correct BA information and predict oligomeric state for BA structrues in the PDB (see", a("Usage", href="#usage")," for more detailed information).")



),
tabPanel('Usage', id = "usage",

h4("Workflow of the tool"),

tags$img(src = "screenShots/workFlow.pdf", width = "100%"),

h4("Usage of the tool"),

tags$ol(
    tags$li("First, upload a data set using", tags$b("Data upload"), "tab. Data set requires PDB ID and PMC/PubMed ID. Please note that, even though the PMC ID has a 'PMC' prefix, the ID must be without that prefix (i.e. 1616964 instead of PMC1616964, see example dataset). There are three options in the", tags$b("Data upload"), "tab. Users can load example data set to test the tool or upload their own dataset using either", tags$b("Single PDB entry"), "or", tags$b("Upload a file"), "options."),

br(),
    tags$img(src = "screenShots/ss2.jpg", width = "100%"),


br(),
br(),


    tags$li("After uploading the dataset properly, click", tags$b("Text mining"), "tab. Select one of three machine learning algorithms from the side-panel on the left: Support Vector Machines, Boosted Logistic Regression and Random Forest. Finally, click submit to start text mining."),
br(),

tags$img(src = "screenShots/ss3.jpg", width = "100%"),

br(),
br(),

    tags$li("After text mining process is done, the results will be appeared under six subtabs in this section.", tags$b("Summary"), "displays a summary table of the results. Users can overview the results quickly, such as number of articles, number of articles founded as HTML, number of articles without oligomeric state, number of articles without full text or abstract and number of unique PDB ID."),
br(),
    tags$img(src = "screenShots/summary.jpg", width = "100%"),
br(),
br(),

tags$li("Users can see the extracted sentences from the primary citations under", tags$b("Sentences"),"sub-tab. These sentences extracted from papers based on a keyword list (see full list in the", tags$b("Keywords"), "sub-tab under the", tags$b("Help"),"tab). The first column shows the PDB IDs. These IDs are also links for the corresponding RCSB web-site. Second column shows the hit keyword as oligomeric state and last column shows the corresponding sentence."),
br(),
tags$img(src = "screenShots/sentences.jpg", width = "100%"),
br(),
br(),

tags$li("Since we use a set of biological assembly related keyword list, all extracted sentences will be related with the biological assembly information somehow. However, some of these sentences are redundant. For example, some sentences will be related with the assymetric unit rather than the biological assembly and some sentences will be about some other protein structures and not the one we interested in. Therefore, we apply machine learning algorithms to eliminate these redundant sentences. Users can select one of three machine learning algorithms, including Support Vector Machines, Boosted Logistic Regression and Random Forest, at the begining of the text mining analysis (see the model performances of the machine learning algorithms in the", tags$b("Model performance "), "tab under the",  tags$b("Help"),"tab). In this section, there will be an extra column,", tags$b("Label"),", which classify these sentecences as biological assembly related (TRUE) and biological assembly unrelated (FALSE)."),
br(),
tags$img(src = "screenShots/ml.jpg", width = "100%"),
br(),
br(),

tags$li("Finally, click", tags$b("Oligomeric state prediction"), "sub-tab to see the oligomeric state predictions for each PDB ID. The probability will be calculated based on the majority of the remaining sentences after machine learning filtering. PDB ID and PMC ID columns are also links for the corresponding protein structre in the RCSB website and the papers in the PMC or PubMed databases."),
br(),
tags$img(src = "screenShots/ss4.jpg", width = "100%"),
br(),
br(),

tags$li("Since we use a keyword list, there maybe some papers which have a full text or abstact but not the keyword that we are using for search and extract sentences. In this case users can use", tags$b("No oligomeric state "), "sub-tab and easily reach that papers using the referring links."),
br(),
tags$img(src = "screenShots/nor.jpg", width = "100%"),
br(),
br(),
tags$li("Since some papers are old enough, there may not be available full texts or even abstracts. These papers will be appeared in the", tags$b("No full text/abstract"), "sub-tab with their corresponding PDB ID and PMC or PubMed links."),

tags$img(src = "screenShots/nftar.jpg", width = "100%")






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



            tabPanel('Keywords',h5("These keywords will be used to search and extract biological assembly related sentences throughout the paper(s)."),

            verbatimTextOutput("baKeywords")

            ),

            tabPanel('News',HTML('<p><b>This page will be available soon...</b></p>')),
            tabPanel('Authors',HTML('<p><b>This page will be available soon...</b></p>')),
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










