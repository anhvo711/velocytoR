#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


options(shiny.maxRequestSize = 800*1024^2)

GoToTab <- function(name){
    #updateTabItems(session, "tabs", name)
    
    shinyjs::show(selector = paste0("a[data-value=\"",name,"\"]"))
    
    shinyjs::runjs("window.scrollTo(0, 0)")
}

# Define server logic
server <- function(input, output, session){
    values = reactiveValues(
        ldat = NULL,
        emat = NULL,
        nmat = NULL,
        smat = NULL,
        pagodaObject = NULL,
        emb = NULL,
        clusterlabel = NULL,
        celldist = NULL,
        rvel = NULL,
        cellcolors = NULL,
        velocity = NULL
    )
    
    observe ({
        shinyjs::hide(selector = "a[data-value=\"normalizeFilterTab\"]")
        shinyjs::hide(selector = "a[data-value=\"pagodaProcessTab\"]")
        shinyjs::hide(selector = "a[data-value=\"velocityEstimationTab\"]")
        shinyjs::hide(selector = "a[data-value=\"visualizationGene\"]")
        shinyjs::hide(selector = "a[data-value=\"cellTrajectoryTab\"]")
    })
    
    #DATA INPUT TAB ---------------------------------------------------------------------------------------
    
    observeEvent(input$inFile, {
        withProgress(message = "Reading loom file", {
            values$ldat = read.loom.matrices(input$inFile$datapath)
        })
    })
    
    histogramPlot = function() {
        hist(log10(colSums(values$emat)))
    }

    output$colSumsHistPlot <- renderPlot({
        req(input$inFile)
        values$emat <- values$ldat$spliced # exonic read (spliced) expression matrix
        updateSelectInput(session, "filterSpecGenes", choices = rownames(values$emat))
        values$nmat <- values$ldat$unspliced; # intronic read (unspliced) expression matrix
        #values$smat <- values$ldat$spanning; # spanning read (intron+exon) expression matrix
        histogramPlot()
    })
    
    #when file is chosen, show colSumsHistPlot with spinner
    output$fileUploaded <- reactive({
        return(!is.null(values$ldat))
    })
    outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)

    #button clicked to go to normalize & filter tab
    observeEvent(input$goToNormalize, {
        GoToTab("normalizeFilterTab")
        shinyjs::show(selector = "a[data-value=\"normalizeFilterTab\"]")
    })
    
    #plot download
    output$downloadHistogram <- downloadHandler(
        filename <- function() {
            paste("velocytoR_histogram", input$histogramDownloadAs, sep="")
        },
        content <- function(file) {
            if (input$histogramDownloadAs == ".svg") {
                svg(file, width = input$histogramDownloadWidth, height = input$histogramDownloadHeight)
            } else if (input$histogramDownloadAs == ".pdf") {
                pdf(file, width = input$histogramDownloadWidth, height = input$histogramDownloadHeight)
            } else if (input$histogramDownloadAs == ".jpeg") {
                jpeg(file, width = input$histogramDownloadWidth, height = input$histogramDownloadHeight, units = "in", res = 300)
            } else {
                png(file, width = input$histogramDownloadWidth, height = input$histogramDownloadHeight, units = "in", res = 300)
            }
            histogramPlot()
            dev.off()
        },
        contentType = "image"
    )
    
    #NORMALIZE & FILTER TAB --------------------------------------------------------------------------------

    #when filter button is clicked, filter and make new Pagoda2 object
    observeEvent(input$filterButton, {
        values$emat <- values$ldat$spliced
        values$emat <- values$emat[,colSums(values$emat)>=1e3]
        if (!is.null(input$filterSpecGenes)) {
                    values$emat = values$emat[!rownames(values$emat) %in% c(input$filterSpecGenes), ]
                    print("filtered by genes")
        }
        if (input$regexFilter != "") {
                    values$emat = values$emat[!grepl(input$regexFilter, rownames(values$emat)), ]
                    print("filteres by regex")
        }
        values$pagodaObject <- Pagoda2$new(values$emat, trim=10, log.scale=T)
    })
    
    #show adjustedVariancePlot output with spinner when we have a pagoda object
    output$filterButtonClicked <- reactive({
        return(!is.null(values$pagodaObject))
    })
    outputOptions(output, 'filterButtonClicked', suspendWhenHidden=FALSE)
    
    variancePlot <- function() {
        values$pagodaObject$adjustVariance(plot=T, do.par=T, gam.k=10)
    }
    
    output$adjustedVariancePlot <- renderPlot({
        req(values$pagodaObject)
        input$filterButton
        variancePlot()
    })
    
    #plot download
    output$downloadVariance <- downloadHandler(
        filename <- function() {
            paste("velocytoR_adjustedVariance", input$varianceDownloadAs, sep="")
        },
        content <- function(file) {
            if (input$histogramDownloadAs == ".svg") {
                svg(file, width = input$varianceDownloadWidth, height = input$varianceDownloadHeight)
            } else if (input$histogramDownloadAs == ".pdf") {
                pdf(file, width = input$varianceDownloadWidth, height = input$varianceDownloadHeight)
            } else if (input$histogramDownloadAs == ".jpeg") {
                jpeg(file, width = input$varianceDownloadWidth, height = input$varianceDownloadHeight, units = "in", res = 300)
            } else {
                png(file, width = input$varianceDownloadWidth, height = input$varianceDownloadHeight, units = "in", res = 300)
            }
            variancePlot()
            dev.off()
        },
        contentType = "image"
    )
    
    observeEvent(input$goToProcessing, {
        GoToTab("pagodaProcessTab")
        shinyjs::show(selector = "a[data-value=\"pagodaProcessTab\"]")
        
    })
    
    
    #PAGODA2 PROCESSING TAB ------------------------------------------------------------------------------
    
    #button clicked, show embeddingPlot output wit spinner while processing
    observeEvent(input$showEmbeddingPlotButton, {
        withProgress(message = "Processing...", value = 40, {
            values$pagodaObject$calculatePcaReduction(nPcs=input$numPC, n.odgenes=input$numODGenes, maxit=300)
            incProgress(10)
            values$pagodaObject$makeKnnGraph(k=input$numK, type='PCA', center=input$centerBoolean, distance='cosine')
            incProgress(10)
            values$pagodaObject$getKnnClusters(method=multilevel.community, type='PCA', name='multilevel')
            incProgress(10)
            values$pagodaObject$getEmbedding(type='PCA', embeddingType='tSNE' ,perplexity=50, verbose=T)
            incProgress(10)
            par(mfrow=c(1,2))
        })
    })
    
    output$showEmbeddingButtonClicked <- reactive({
        input$showEmbeddingPlotButton
        return(input$showEmbeddingPlotButton > 0)
    })
    outputOptions(output, 'showEmbeddingButtonClicked', suspendWhenHidden=FALSE)
    
    embed <- function() {
        values$pagodaObject$plotEmbedding(type='PCA', embeddingType='tSNE', show.legend=input$showLegend, 
                                          mark.clusters=input$labelClusters, min.group.size=10, shuffle.colors=F, mark.cluster.cex=input$labelClustersSize,
                                          alpha=0.3, main=input$plotTitle, legend.x = input$legendPos)
    }
    
    output$embeddingPlot <- renderPlot({
        input$showEmbeddingPlotButton
        embed()
    })
    
    #plot download
    output$downloadEmbedding <- downloadHandler(
        filename <- function() {
            paste("velocytoR_embeddingPlot", input$embeddingDownloadAs, sep="")
        },
        content <- function(file) {
            if (input$embeddingDownloadAs == ".svg") {
                svg(file, width = input$embeddingDownloadWidth, height = input$embeddingDownloadHeight)
            } else if (input$embeddingDownloadAs == ".pdf") {
                pdf(file, width = input$embeddingDownloadWidth, height = input$embeddingDownloadHeight)
            } else if (input$embeddingDownloadAs == ".jpeg") {
                jpeg(file, width = input$embeddingDownloadWidth, height = input$embeddingDownloadHeight, units = "in", res = 300)
            } else {
                png(file, width = input$embeddingDownloadWidth, height = input$embeddingDownloadHeight, units = "in", res = 300)
            }
            embed()
            dev.off()
        },
        contentType = "image"
    )
    
    observeEvent(input$goToVelocity, {
        GoToTab("velocityEstimationTab")
        shinyjs::show(selector = "a[data-value=\"velocityEstimationTab\"]")
        
    })
    
    #VELOCITIES TAB ----------------------------------------------------------------------------------------
    #when button is clicked, do premliminary processing
    observeEvent(input$calculateVelocitiesButton, {
        if (input$calculateVelocitiesButton == 1) {
            withProgress(message = "Processing...", value = 20, {
                values$emat <- values$emat[,rownames(values$pagodaObject$counts)]
                incProgress(10)
                values$nmat <- values$nmat[,rownames(values$pagodaObject$counts)]; # restrict to cells that passed p2 filter
                incProgress(10)
                values$clusterlabel <- values$pagodaObject$clusters$PCA$multilevel # take the cluster factor that was calculated by p2
                incProgress(10)
                values$cellcolors <- pagoda2:::fac2col(values$clusterlabel)
                incProgress(10)
                values$emb <- values$pagodaObject$embeddings$PCA$tSNE
                incProgress(10)
                values$celldist <- as.dist(1-armaCor(t(values$pagodaObject$reductions$PCA)))
            })
        }
        withProgress(message = "Processing...", value = 20, {
            incProgress(10)
            values$emat <- filter.genes.by.cluster.expression(values$emat, values$clusterlabel, min.max.cluster.average = input$ematClusterAvg)
            incProgress(10)
            values$nmat <- filter.genes.by.cluster.expression(values$nmat, values$clusterlabel, min.max.cluster.average = input$nmatClusterAvg)
            incProgress(10)
            x = input$quantile
            if (!is.null(values$rvel)){
                values$rvel <- gene.relative.velocity.estimates(values$emat, values$nmat, deltaT=1, kCells=input$numKnnCells1, 
                                                                cell.dist=values$celldist, fit.quantile= x, old.fit = values$rvel, n.cores = 16)
            } else {
                values$rvel <- gene.relative.velocity.estimates(values$emat, values$nmat, deltaT=1, kCells=input$numKnnCells1, 
                                                            cell.dist=values$celldist, fit.quantile= x, n.cores = 16)
            }
            
            updateSelectInput(session, "gene", choices = rownames(values$rvel$deltaE))
        })
    })
    
    #allows for velocityOnEmbeddingPlot to be shown with spinner
    output$calculateVelocitiesButtonClicked <- reactive({
        return(!is.null(values$clusterlabel))
    })
    outputOptions(output, 'calculateVelocitiesButtonClicked', suspendWhenHidden=FALSE)
    
    velocitiesPlot <- function() {
       values$velocity = show.velocity.on.embedding.cor(values$emb, values$rvel , n=input$neighborhoodVelocity, scale=input$scaleVelocity, cell.colors=ac(values$cellcolors,alpha=0.5), 
                                       cex=input$pointSize, arrow.scale=input$arrowScale, show.grid.flow=TRUE, min.grid.cell.mass=0.5, grid.n=40, arrow.lwd=1,
                                       do.par=T,cell.border.alpha = 0.1,cc= values$velocity$cc, n.cores = 16)
    }
    
    output$velocityOnEmbeddingPlot <- renderPlot({
        req(values$rvel)
        input$calculateVelocitiesButton #dependency on the button
        
        withProgress(message = "Plotting velocities...", value = 30, {
            isolate(velocitiesPlot())
        })
    })
    
    #plot download
    output$downloadVelocities <- downloadHandler(
        filename <- function() {
            paste("velocytoR_velocitiesPlot", input$velocitiesDownloadAs, sep="")
        },
        content <- function(file) {
            if (input$velocitiesDownloadAs == ".svg") {
                svg(file, width = input$velocitiesDownloadWidth, height = input$velocitiesDownloadHeight)
            } else if (input$velocitiesDownloadAs == ".pdf") {
                pdf(file, width = input$velocitiesDownloadWidth, height = input$velocitiesDownloadHeight)
            } else if (input$velocitiesDownloadAs == ".jpeg") {
                jpeg(file, width = input$velocitiesDownloadWidth, height = input$velocitiesDownloadHeight, units = "in", res = 300)
            } else {
                png(file, width = input$velocitiesDownloadWidth, height = input$velocitiesDownloadHeight, units = "in", res = 300)
            }
            withProgress(message = "Prepping plot for download...", {
                velocitiesPlot()
            })
            dev.off()
        },
        contentType = "image"
    )
    
    observeEvent(input$goToGeneVisualization, {
        GoToTab("visualizationGene")
        shinyjs::show(selector = "a[data-value=\"visualizationGene\"]")
    })
    
    
    #GENE VISUALIZATION -------------------------------------------------------------------------------------------------
    
    #when button is clicked, show plot output space
    output$visualizeGeneButtonClicked <- reactive({
        input$visualizeGene
        return(input$visualizeGene > 0)
    })
    outputOptions(output, 'visualizeGeneButtonClicked', suspendWhenHidden=FALSE)
    
    geneVis <- function() {
        gene.relative.velocity.estimates(values$emat, values$nmat, deltaT=1, kCells=input$numKnnCells2, kGenes=1, cell.emb=values$emb,
                                         cell.dist=values$celldist, fit.quantile= 0.02, show.gene=input$gene, old.fit=values$rvel, do.par=T, n.cores = 16)    
    }
    
    output$geneVisual <- renderPlot({
        input$visualizeGene
        isolate(geneVis())
    })
    
    #plot download
    output$downloadGene <- downloadHandler(
        filename <- function() {
            paste("velocytoR_geneVisualization", input$geneDownloadAs, sep="")
        },
        content <- function(file) {
            if (input$geneDownloadAs == ".svg") {
                svg(file, width = input$geneDownloadWidth, height = input$geneDownloadHeight)
            } else if (input$geneDownloadAs == ".pdf") {
                pdf(file, width = input$geneDownloadWidth, height = input$geneDownloadHeight)
            } else if (input$geneDownloadAs == ".jpeg") {
                jpeg(file, width = input$geneDownloadWidth, height = input$geneDownloadHeight, units = "in", res = 300)
            } else {
                png(file, width = input$geneDownloadWidth, height = input$geneDownloadHeight, units = "in", res = 300)
            }
            withProgress(message = "Prepping plot for download...", {
                geneVis()
            })
            dev.off()
        },
        contentType = "image"
    )
    
    observeEvent(input$goToGeneTrajectoryModeling, {
        GoToTab("cellTrajectoryTab")
        shinyjs::show(selector = "a[data-value=\"cellTrajectoryTab\"]")
    })
    
    #TRAJECTORY MODELING -------------------------------------------------------------------------------------------------
    #when button is clicked, show plot output space
    calcTrajectory <- function() {
        if (input$PCAReductionBoolean == TRUE){
            npc = input$numPC
        } else {
            npc = "NA"
        }
            show.velocity.on.embedding.eu(values$emb, values$rvel, n=input$neighborhoodTrajectory, scale=input$scaleTrajectory, cell.colors=ac(values$cellcolors,alpha=0.5), 
                                          cex=0.8, nPcs=npc, sigma=input$sigma, show.trajectories=TRUE, diffusion.steps=input$diffusionStep, 
                                          n.trajectory.clusters=input$trajectoryClusters, ntop.trajectories=1, embedding.knn=T, control.for.neighborhood.density=TRUE,n.cores = 16)
    }
    
    output$trajectory <- renderPlot({
        req(input$findTrajectory > 0)
        input$findTrajectory
        withProgress(message = "Finding trajectory...", {
            isolate(calcTrajectory())
        })
    })
    
    #creates output variable that can be seen by conditional panel to show plot and download section
    output$findTrajectoryClicked <- reactive({
        input$findTrajectory
        return(input$findTrajectory > 0)
    })
    outputOptions(output, "findTrajectoryClicked", suspendWhenHidden=FALSE)
    
    #shows numPC option is PCA reduction checkbox is checked
    output$pcaReduction <- renderUI({
        if (input$PCAReductionBoolean) {
            numericInput("numPC", "Number of PCs to Project Cells onto", value = 30, min = 1, max = 100)
        }
    })
    
    #plot download
    output$downloadTrajectory <- downloadHandler(
        filename <- function() {
            paste("velocytoR_cellTrajectoryModeling", input$trajectoryDownloadAs, sep="")
        },
        content <- function(file) {
            if (input$trajectoryDownloadAs == ".svg") {
                svg(file, width = input$trajectoryDownloadWidth, height = input$trajectoryDownloadHeight)
            } else if (input$trajectoryDownloadAs == ".pdf") {
                pdf(file, width = input$trajectoryDownloadWidth, height = input$trajectoryDownloadHeight)
            } else if (input$trajectoryDownloadAs == ".jpeg") {
                jpeg(file, width = input$trajectoryDownloadWidth, height = input$trajectoryDownloadHeight, units = "in", res = 300)
            } else {
                png(file, width = input$trajectoryDownloadWidth, height = input$trajectoryDownloadHeight, units = "in", res = 300)
            }
            withProgress(message = "Prepping plot for download...", {
                calcTrajectory()
            })
            dev.off()
        },
        contentType = "image"
    )
}


