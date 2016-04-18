correlationTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    #currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    currentCorrelation <- allReactiveVars$currentCorrelation
        
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
        if (!is.null(currentMetadata$final$alt_id))
			colnames(tab) <- as.character(currentMetadata$final$alt_id)
        
        switch(input$rnaCorrelateWhat,
			samples = {
				currentCorrelation$corMatrix <- cor(tab,
					method=input$rnaCorrelationMethod)
			},
			allgenes = {
				currentCorrelation$corMatrix <- cor(t(tab),
					method=input$rnaCorrelationMethod)
			},
			refgene = {
			}
		)
        
        currentCorrelation$opts$method <- input$rnaCorrelationMethod
        currentCorrelation$opts$colors <- c(
			input$corrColourHigh,
			input$corrColourNo,
			input$corrColourLow
        )
        
        
        
        #n <- dim(cor.mat)[1]
        #labs <- matrix(NA,n,n)
        #for (i in 1:n)
        #    for (j in 1:n)
        #        labs[i,j] <- sprintf("%.2f",cor.mat[i,j])
        #if (n <= 5)
        #    notecex <- 1.2
        #else if (n > 5 & n < 10)
        #    notecex <- 0.9
        #else
        #    notecex <- 0.7
        #heatmap.2(
		#	cor.mat,
		#	col=colorRampPalette(c("yellow","grey","blue")),
        #    revC=TRUE,
        #    trace="none",
        #    symm=TRUE,
        #    Colv=TRUE,
        #    cellnote=labs,
        #    keysize=1,
        #    density.info="density",
        #    notecex=notecex,
        #    cexCol=0.9,
        #    cexRow=0.9,
        #    font.lab=2
        #)
        
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
        if (is.null(currentCorrelation$corMatrix)) {
            output$correlation <- renderPlot({
                currentCorrelation$entry
            })
            plotOutput("correlation",height="640px")
        }
        else {
			n <- dim(currentCorrelation$corMatrix)[1]
			labs <- matrix(NA,n,n)
			for (i in 1:n)
				for (j in 1:n)
					labs[i,j] <- sprintf("%.2f",
						currentCorrelation$corMatrix[i,j])
			output$heatmap <- renderD3heatmap({
				d3heatmap(
					currentCorrelation$corMatrix,
					dendrogram="both",
					Rowv=TRUE,
					Colv=TRUE,
					colors=colorRampPalette(currentCorrelation$opts$colors),
					revC=TRUE,
					cellnote=labs,
					symm=TRUE
				)
			})
			div(
				class="heatmap-container",
				d3heatmapOutput("heatmap",height="800px")
			)
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
    
    performRnaCorrelation <- 
        correlationTabPanelReactiveEvents$performRnaCorrelation
    
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
            g <- isolate({input$selectCorrelationGeneName})
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
            g <- isolate({input$selectCorrelationGeneName})
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
    
    observe({
        tryCatch({
            shinyjs::disable("performRnaCorrelation")
            performRnaCorrelation()
        },error=function(e) {
            print(e)
        },
        finally={
			shinyjs::enable("performRnaCorrelation")
		})
    })
}

