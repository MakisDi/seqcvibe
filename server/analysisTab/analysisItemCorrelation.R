correlationTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    currentHeatmap <- allReactiveVars$currentHeatmap
        
    performRnaCorrelation <- eventReactive(input$performRnaCorrelation,{
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        meta <- currentMetadata$final
        samples <- as.character(meta$sample_id)
        D <- loadedData[[s]][[d]]
        genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
        
        # Determine which genes to use
        switch(input$rnaCorrelationGeneList,
            select = {
                g <- input$selectCorrelationGeneName
            },
            custom = {
                g <- input$rnaCorrelationCustomList
                g <- strsplit(g,split="\n")[[1]]
                genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
                m <- match(g,names(genes))
                na <- which(is.na(m))
                if (length(na)>0)
                    m <- m[-na]
                g <- genes[m]
            },
            all = {
                bad <- apply(D$norm,1,function(x) { return(all(x==0)) })
                g <- genes[-which(bad)]
            }
        )
        
        # Determine the measurements table
        switch(input$rnaCorrelationMeasureRadio,
            counts = {
                tab <- D$norm[g,samples,drop=FALSE]
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    tab <- rbind(tab,A)
                }
            },
            rpkm = {
                tab <- round(edgeR::rpkm(
                    D$counts[g,samples,drop=FALSE],
                    gene.length=D$length[g],
                    lib.size=unlist(D$libsize[samples])
                    ),digits=6)
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    A <- round(edgeR::rpkm(A,
                        gene.length=currentCustomRnaTables$lengths,
                        lib.size=unlist(D$libsize[samples])
                    ),digits=6)
                    tab <- rbind(tab,A)
                }
            },
            rpgm = {
                tab <- round(D$norm[g,
                    samples,drop=FALSE]/D$length[g],
                        digits=6)
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    A <- round(A/currentCustomRnaTables$lengths,
                        digits=6)
                    tab <- rbind(tab,A)
                }
            }
        )
        
        # We now have the tab, transformations based on other selections
        switch(input$rnaCorrelationScaleRadio,
            natural = {
                tab <- tab
            },
            log2 = {
                tab <- round(log2(tab+1),digits=6)
            }
        )
        
    })
    
    return(list(
        performRnaCorrelation=performRnaCorrelation
    ))
}

correlationTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}

correlationTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentCorrelation <- allReactiveVars$currentCorrelation
    
    output$correlationOutput <- renderUI({
        if (is.null(currentCorrelation$data)) {
            output$correlation <- renderPlot({
                currentCorrelation$entry
            })
            plotOutput("correlation",height="640px")
        }
        else {
        }
    })
    
    output$exportCorrelationPDF <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("rna_correlation_",tt,".pdf", sep='')
        },
        content=function(con) {
        }
    )
    
    output$exportCorrelationPNG <- downloadHandler(
        filename=function() {
            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
            paste("rna_correlation_",tt,".png", sep='')
        },
        content=function(con) {
        }
    )
}

correlationTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentCorrelation <- allReactiveVars$currentHeatmap
    
    correlationTabPanelReactiveEvents <- 
        correlationTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    performRnaClustering <- 
        correlationTabPanelReactiveEvents$performRnaClustering
    
    correlationTabPanelReactiveExprs <- 
        correlationTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    correlationTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
        
    observe({
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"selectCorrelationGeneName",
                choices=NULL,
                server=TRUE
            )
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            g <- isolate({input$selectClusteringGeneName})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"selectCorrelationGeneName",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
    
    observe({
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"rnaCorrelationRefGene",
                choices=NULL,
                server=TRUE
            )
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            g <- isolate({input$selectClusteringGeneName})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"rnaCorrelationRefGene",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
}

