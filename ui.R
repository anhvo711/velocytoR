#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(velocyto.R)
library(pagoda2)
library(shinycssloaders)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- tagList(
  dashboardPage(skin = "purple",
                dashboardHeader(title = "velocytoR"),
                dashboardSidebar(
                  sidebarMenu(
                    id = "tabs",
                    menuItem("User Guide", tabName = "introTab", icon = icon("info-circle")),
                    menuItem("Data Input", tabName = "dataInputTab", icon = icon("upload")),
                    menuItem("Normalize & Filter", tabName = "normalizeFilterTab", icon = icon("filter")),
                    menuItem("Pagoda2 Processing", tabName = "pagodaProcessTab", icon = icon("th")),
                    menuItem("Velocity Estimation", tabName = "velocityEstimationTab", icon = icon("chart-line")),
                    menuItem("Gene Visualization", tabName = "visualizationGene", icon = icon("chart-line")),
                    menuItem("Cell Trajectory Modeling", tabName = "cellTrajectoryTab", icon = icon("chart-line"))
                  )
                ),
                dashboardBody(
                  shinyjs::useShinyjs(),
                  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
                  tabItems(
                    tabItem(tabName = "introTab",
                            p("Velocyto is a library for the analysis of RNA velocity."),
                            p("This application is a web interface to the velocyto.py ", 
                              a("analysis pipeline.", href="http://velocyto.org/velocyto.py/tutorial/analysis.html#analysis")),
                            br(),
                            tags$h2("1. Data Input"),
                            p("To start, upload a .loom file to the webpage and wait for summarized plot to appear. 
                              Here, you can download the plot or continue on to the next step of the pipeline."),
                            p("Note: Only .loom files are accepted. If you need to generate a .loom file, please contact the GCBC."),
                            br(),
                            tags$h2("2. Normalize & Filter"),
                            p("The data can be filtered by excluding genes. There are two methods: regex expression and gene selection. Click \"Filter\" to filter the data 
                              and show the adjusted plot."),
                            p("If no filtering is needed, leave the fields blank and clikc \"Filter\" to show the adjusted plot after splicing from the loom file."),
                            br(),
                            tags$h2("3. Pagoda2 Processing"),
                            p("The processing step uses a Pagoda2 object and consists of 3 steps: PCA reduction, creation of a KNN graph, and creation of the embedding plot.
                              Fill out the fields accordingly and click \"Find Embedding\" to show a plot of the cell clusters. Preferences for this plot can be changed under \"Plot Embedding\""),
                            br(),
                            tags$h2("4. Velocity Estimation"),
                            p("First, choose options for calculating the velocities. Then click \"Calculate Velocities\" to display the velocities plot. This plot
                            can be downloaded inv arious formats."),
                            br(),
                            tags$h2("5. Gene Visualization"),
                            p("From the list of genes taken from the .loom file, choose a gene to view its various plots. The number of KNN cells to be used to in plotting can be changed as needed."),
                            br(),
                            tags$h2("6. Cell Trajectory Modeling"),
                            p("This tab can be used to model central trajectories by directed diffusion on embedding. The main parameters are set up by sigma (which limits the range of how far 
                              a cell can jump in terms of distance) and n (how many nearest neighbors are being considered for jumps). The results are sensitive to these parameters, 
                              as we don’t have a good way of assessing how much the directional velocity component should compare with random Brownian motion of a cell with the manifold."),
                            p("Warning: This step can take some time.")
                    ),
                    tabItem(tabName = "dataInputTab",
                            fileInput("inFile", "Choose a Loom File", FALSE, accept = c(".loom")),
                            fluidRow(
                              column(8,
                                     withSpinner(plotOutput("colSumsHistPlot", height = "50vh"), type = 3, color.background = "#ecf0f5")
                              ),
                              column(4,
                                     conditionalPanel(condition = "output.fileUploaded",
                                                      h4("Plot Download Options"),
                                                      numericInput("histogramDownloadHeight", "Plot height (in inches):", value = 15),
                                                      numericInput("histogramDownloadWidth", "Plot width (in inches):", value = 15),
                                                      radioButtons("histogramDownloadAs", "Download File Type:", choices = list("PDF" = ".pdf", "SVG" = ".svg", "PNG" = ".png", "JPEG" = ".jpeg"), selected = ".pdf"),
                                                      downloadButton("downloadHistogram", "Download Plot"),
                                                      br(),
                                                      actionButton("goToNormalize", "Next Step", class = "next", icon = icon("arrow-circle-right"))
                                     )
                              )
                            )
                            
                    ),
                    tabItem(tabName = "normalizeFilterTab",
                            textInput("regexFilter", "Remove Genes that Match This Regex Expression", placeholder = "Eg. ^MT- for genes that start with 'MT-'"),
                            span("To use multiple regex expressions, separate each expression by |. (Eg. ^MT|^Rb)"),
                            br(),
                            br(),
                            selectizeInput("filterSpecGenes", label="Select Genes to Remove", choices= NULL, multiple=TRUE),
                            actionButton("filterButton", "Filter Data"),
                            tags$hr(),
                            conditionalPanel(condition = "output.filterButtonClicked",
                                             fluidRow(
                                               br(),
                                               column(8,
                                                      withSpinner(plotOutput("adjustedVariancePlot", height = "50vh"), type = 3, color.background = "#ecf0f5")
                                               ),
                                               column(4,
                                                      h4("Plot Download Options"),
                                                      numericInput("varianceDownloadHeight", "Plot height (in inches):", value = 15),
                                                      numericInput("varianceDownloadWidth", "Plot width (in inches):", value = 15),
                                                      radioButtons("varianceDownloadAs", "Download File Type:", choices = list("PDF" = ".pdf", "SVG" = ".svg", "PNG" = ".png", "JPEG" = ".jpeg"), selected = ".pdf"),
                                                      downloadButton("downloadVariance", "Download Plot"),
                                                      br(),
                                                      actionButton("goToProcessing", "Next Step", class = "next", icon = icon("arrow-circle-right"))
                                               )
                                             )
                            )
                    ),
                    tabItem(tabName = "pagodaProcessTab",
                            fluidRow(
                              column(4,
                                     h3("Calculate PCA Reduction"),
                                     numericInput("numPC", "Number of PCs to Use", value = 100, max = 500, min = 20),
                                     numericInput("numODGenes", "Number of Top Overdispersed Genes to be Used", value = 3000, max = 30000, min = 1000)
                              ),
                              column(4,
                                     h3("Make KNN Graph"),
                                     numericInput("numK", "Number of Nearest Neighbor", value = 30, max = 100, min = 10),
                                     selectInput("centerBoolean", "Center Data before Making Graph", choices = c("Yes"="TRUE", "No"="FALSE"), selected = "Yes")
                              ),
                              column(4,
                                     h3("Plot Embedding"),
                                     selectInput("showLegend", "Show Graph Legend", choices = c("Yes"="TRUE", "No"="FALSE"), selected = "Yes"),
                                     conditionalPanel(condition = "input.showLegend == 'TRUE'",
                                                      selectInput("legendPos", "Legend Position", choices = c("topright", "topleft", "bottomleft", "bottomright"))
                                     ),
                                     selectInput("labelClusters", "Label Clusters", choices = c("Yes"="TRUE", "No"="FALSE"), selected = "Yes"),
                                     conditionalPanel(condition = "input.labelClusters == 'TRUE'",
                                                      numericInput("labelClustersSize", "Size of Cluster Labels", value = 2, max = 5, min = 1)
                                                      
                                     ),
                                     textInput("plotTitle", "Enter Plot Title")
                              )
                            ),
                            actionButton("showEmbeddingPlotButton", "Find Embedding"),
                            hr(),
                            conditionalPanel(condition = "output.showEmbeddingButtonClicked",
                                             br(),
                                             fluidRow(
                                               column(8,
                                                      withSpinner(plotOutput("embeddingPlot", height = "50vh"), type = 3, color.background = "#ecf0f5")
                                               ),
                                               column(4,
                                                      h4("Plot Download Options"),
                                                      numericInput("embeddingDownloadHeight", "Plot height (in inches):", value = 15),
                                                      numericInput("embeddingDownloadWidth", "Plot width (in inches):", value = 15),
                                                      radioButtons("embeddingDownloadAs", "Download File Type:", choices = list("PDF" = ".pdf", "SVG" = ".svg", "PNG" = ".png", "JPEG" = ".jpeg"), selected = ".pdf"),
                                                      downloadButton("downloadEmbedding", "Download Plot"),
                                                      br(),
                                                      actionButton("goToVelocity", "Next Step", class = "next", icon = icon("arrow-circle-right"))
                                               )
                                             )
                            )
                    ),
                    tabItem(tabName = "velocityEstimationTab",
                            fluidRow(
                              column(3,
                                     numericInput("ematClusterAvg", "Minimum Average Expresion Magnitude for Exonic Reads", value = 0.2, min = 0.01, max = 15, step = 0.01),
                                     numericInput("nmatClusterAvg", "Minimum Average Expresion Magnitude for Intronic Reads", value = 0.05, min = 0.01, max = 15, step = 0.01)
                              ),
                              column(3,
                                     numericInput("numKnnCells1", "Number of kNN Cells to Use", value = 25, max = 200, min = 5),
                                     numericInput(inputId = "quantile", label = "% Quantile to use for Gamma Fit", value = 0.02, min = 0.01, max = 99.0, step = 0.01),
                                     selectInput("scaleVelocity", "Select Velocity Scale to Use", choices = c("log", "sqrt", "rank", "linear"), selected = "sqrt")
                              ),
                              column(3,
                                     numericInput("neighborhoodVelocity", "Neighborhood Size", value = 100, min = 1, max = 500),
                                     numericInput("pointSize", "Dot Point Size", value = 0.8, min = 0.1, max = 0.9, step = 0.1),
                                     numericInput("arrowScale", "Arrow Scale Multiplier", value = 2, min = 1, max = 10)
                              )
                            ),
                            actionButton("calculateVelocitiesButton", "Calculate Velocities"),
                            hr(),
                            conditionalPanel(condition = "output.calculateVelocitiesButtonClicked",
                                             br(),
                                             fluidRow(
                                               column(8,
                                                      withSpinner(plotOutput("velocityOnEmbeddingPlot", height = "50vh"), type = 3, color.background = "#ecf0f5")
                                               ),
                                               column(4,
                                                      h4("Plot Download Options"),
                                                      numericInput("velocitiesDownloadHeight", "Plot height (in inches):", value = 10),
                                                      numericInput("velocitiesDownloadWidth", "Plot width (in inches):", value = 15),
                                                      radioButtons("velocitiesDownloadAs", "Download File Type:", choices = list("PDF" = ".pdf", "SVG" = ".svg", "PNG" = ".png", "JPEG" = ".jpeg"), selected = ".pdf"),
                                                      downloadButton("downloadVelocities", "Download Plot"),
                                                      br(),
                                                      actionButton("goToGeneVisualization", "Next Step", class = "next", icon = icon("arrow-circle-right"))
                                               )
                                             )
                            )
                    ),
                    tabItem(tabName = "visualizationGene",
                            selectizeInput("gene", label="Select a Gene to See Visualization", choices= NULL, multiple=FALSE),
                            numericInput("numKnnCells2", "Number of kNN Cells to Use", value = 100, max = 200, min = 5),
                            actionButton("visualizeGene", "Visualize Gene"),
                            hr(),
                            conditionalPanel(condition = "output.visualizeGeneButtonClicked",
                                             br(),
                                             fluidRow(
                                               column(8,
                                                      withSpinner(plotOutput("geneVisual", height = "30vh"), type = 3, color.background = "#ecf0f5")
                                               ),
                                               column(4,
                                                      h4("Plot Download Options"),
                                                      numericInput("geneDownloadHeight", "Plot height (in inches):", value = 5),
                                                      numericInput("geneDownloadWidth", "Plot width (in inches):", value = 15),
                                                      radioButtons("geneDownloadAs", "Download File Type:", choices = list("PDF" = ".pdf", "SVG" = ".svg", "PNG" = ".png", "JPEG" = ".jpeg"), selected = ".pdf"),
                                                      downloadButton("downloadGene", "Download Plot"),
                                                      br(),
                                                      actionButton("goToGeneTrajectoryModeling", "Next Step", class = "next", icon = icon("arrow-circle-right"))
                                               )
                                             )
                            )
                    ),
                    tabItem(tabName = "cellTrajectoryTab",
                            p("Note: This step can take several minutes."),
                            fluidRow(
                              column(3,
                                     selectInput("scaleTrajectory", "Select Velocity Scale to Use", choices = c("log", "sqrt"), selected = "log"),
                                     numericInput("neighborhoodTrajectory", "Neighborhood Size", value = 100, min = 1, max = 500),
                                     checkboxInput("PCAReductionBoolean", "Use PCA Dimensional Reduction", value = TRUE),
                                     uiOutput("pcaReduction")
                              ),
                              column(3,
                                     numericInput("sigma", "Sigma to Use in Calculating Transition Probability", value = 2.5, min = 0.1, max = 100),
                                     numericInput("diffusionStep", "Number of Diffusion Steps to Take Forward", value = 400, min = 5, max = 1000),
                                     numericInput("trajectoryClusters", "Number of Trajectory Clusters to Show Median Paths for", value = 15, min = 1, max = 100)
                              )
                            ),
                            actionButton("findTrajectory", "Find Trajectory"),
                            hr(),
                            conditionalPanel(condition = "output.findTrajectoryClicked",
                                             br(),
                                             fluidRow(
                                               column(8,
                                                      withSpinner(plotOutput("trajectory", height = "30vh"), type = 3, color.background = "#ecf0f5")
                                               ),
                                               column(4,
                                                      h4("Plot Download Options"),
                                                      numericInput("trajectoryDownloadHeight", "Plot height (in inches):", value = 5),
                                                      numericInput("trajectoryDownloadWidth", "Plot width (in inches):", value = 15),
                                                      radioButtons("trajectoryDownloadAs", "Download File Type:", choices = list("PDF" = ".pdf", "SVG" = ".svg", "PNG" = ".png", "JPEG" = ".jpeg"), selected = ".pdf"),
                                                      downloadButton("downloadTrajectory", "Download Plot")
                                               )
                                             )
                            )
                    )
                  )
                )
  )
)
