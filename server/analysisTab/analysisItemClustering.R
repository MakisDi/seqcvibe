clusteringTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
	currentMetadata <- allReactiveVars$currentMetadata
	currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
	currentRnaDeTable <- allReactiveVars$ddcurrentRnaDeTable
	currentHeatmap <- allReactiveVars$currentHeatmap
		
	performRnaClustering <- eventReactive(input$performRnaClustering,{
		s <- currentMetadata$source
        d <- currentMetadata$dataset
        meta <- currentMetadata$final
        samples <- as.character(meta$sample_id)
		D <- loadedData[[s]][[d]]
		genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
		
		# Determine which genes to use
		switch(input$rnaClusterGeneList,
			select = {
				g <- input$selectClusteringGeneName
			},
			custom = {
				g <- input$rnaClusteringCustomList
				g <- strsplit(g,split="\n")[[1]]
				genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
				m <- match(g,names(genes))
				na <- which(is.na(m))
				if (length(na)>0)
					m <- m[-na]
				g <- genes[m]
			},
			degenes = {
				deTable <- currentRnaDeTable$totalTable
				if (is.null(deTable))
					return()
				g <- rownames(deTable)	
			}
		)
		
		# Determine the measurements table
		switch(input$rnaClusteringMeasureRadio,
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
		switch(input$rnaClusteringScaleRadio,
			natural = {
				tab <- tab
			},
			log2 = {
				tab <- round(log2(tab+1),digits=6)
			}
		)
		
		switch(input$rnaClusteringVariablesRadio,
			replicates = {
				X <- tab
				if (!is.null(meta$alt_id))
					colnames(X) <- as.character(meta$alt_id)
			},
			averages = {
				cc <- unique(as.character(meta$class))
				classList <- vector("list",length(cc))
				names(classList) <- cc
				for (cl in cc)
					classList[[cl]] <- 
						as.character(meta$sample_id[which(meta$class==cl)])
				X <- do.call("cbind",lapply(classList,function(x,tab,s,v) {
                    makeStat(x,tab,s,v)
                },tab,"mean","rpkm"))
                colnames(X) <- names(classList)
			},
			fc = {
				cc <- unique(as.character(meta$class))
				classList <- vector("list",length(cc))
				names(classList) <- cc
				for (cl in cc)
					classList[[cl]] <- 
						as.character(meta$sample_id[which(meta$class==cl)])
				control <- input$rnaClusteringControl
				treatments <- setdiff(cc,control)
				contrast <- c(treatments,control)
				contrast <- paste(contrast,collapse="_vs_")
				X <- round(makeFoldChange(contrast,classList,tab,
                    input$rnaClusteringScaleRadio),6)
			}
		)
		
		rownames(X) <- names(genes[match(rownames(X),genes)])
		currentHeatmap$colors <- c(
			input$heatmapColourExtremeDown,
			input$heatmapColourMildDown,
			input$heatmapColourMiddle,
			input$heatmapColourMildUp,
			input$heatmapColourExtremeUp
		)
		currentHeatmap$data <- X
	})
	
	return(list(
		performRnaClustering=performRnaClustering
	))
}

clusteringTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}

clusteringTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
	currentMetadata <- allReactiveVars$currentMetadata
	currentHeatmap <- allReactiveVars$currentHeatmap
	currentRnaDeTable <- allReactiveVars$currentRnaDeTable
	
	output$rnaClusteringFcControl <- renderUI({
        if (is.null(currentMetadata$final))
            list(
                div(style="font-weight:bold","Select control condition"),
                helpText("Please create a dataset first from the 'Data ",
                    "selector' menu on the top")
            )
        else {
            cls <- unique(as.character(currentMetadata$final$class))
            if (length(cls==2))
				list(
					div(style="font-weight:bold","Select control condition"),
					div(class="small",
					helpText("Fold changes cannot be used for clustering with ",
						"only two conditions in the dataset as this would ",
						"produce a 1-column matrix which cannot be clustered.")
					)
				)
			else
				selectInput(
					inputId="rnaClusteringControl",
					label="Select control condition",
					choices=cls
				)
        }
    })
	
	output$heatmapOutput <- renderUI({
		if (is.null(currentHeatmap$data)) {
			output$heatmap <- renderPlot({
				currentHeatmap$heatmap
			})
			plotOutput("heatmap",height="800px")
		}
		else {
			output$heatmap <- renderD3heatmap({
				d3heatmap(
					currentHeatmap$data,
					dendrogram=input$rnaClusterWhat,
					Rowv=TRUE,
					Colv=TRUE,
					distfun=dist,
					hclustfun=hclust,
					k_row=as.numeric(isolate(input$kRowInput)),
					k_col=as.numeric(isolate(input$kColInput)),
					colors=currentHeatmap$colors
				)
				
			})
			d3heatmapOutput("heatmap",height="800px")
		}
	})
}

clusteringTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    
    clusteringTabPanelReactiveEvents <- 
        clusteringTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    performRnaClustering <- 
		clusteringTabPanelReactiveEvents$performRnaClustering
    
    clusteringTabPanelReactiveExprs <- 
        clusteringTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    clusteringTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
        
    observe({
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"selectClusteringGeneName",
                choices=NULL,
                server=TRUE
            )
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            g <- isolate({input$selectClusteringGeneName})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"selectClusteringGeneName",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
    
    observe({
        if (isEmpty(currentMetadata$final))
            shinyjs::disable("performRnaClustering")
        else {
			if (isEmpty(input$selectClusteringGeneName) 
				&& isEmpty(input$rnaClusteringCustomList))
				shinyjs::disable("performRnaClustering")
			else if (input$rnaClusterGeneList=="degenes" 
				&& is.null(currentRnaDeTable$totalTable))
				shinyjs::disable("performRnaClustering")
			else if (input$rnaClusteringVariablesRadio=="fc" 
				&& length(unique(currentMetadata$final$class))<=2)
				shinyjs::disable("performRnaClustering")
			else {
				kr <- as.numeric(input$kRowInput)
				kc <- as.numeric(input$kRowInput)
				if (is.na(kr) || kr<=0 || is.na(kc) || kc<=0) {
					output$rnaClusteringSettingsError <- renderUI({
						div(class="error-message",paste("The number of ",
							"row and column clusters must be a positive ",
							"integer!",sep=""))
					})
					shinyjs::disable("performRnaClustering")
				}
				else {
					output$rnaClusteringSettingsError <- renderUI({div()})
					shinyjs::enable("performRnaClustering")
				}
			}
		}
    })
    
	observe({
        tryCatch({
            shinyjs::disable("performRnaClustering")
            performRnaClustering()
        },error=function(e) {
            print(e)
        },
        finally={
            if (isEmpty(currentMetadata$final))
				shinyjs::disable("performRnaClustering")
			else {
				if (isEmpty(input$selectClusteringGeneName) 
					&& isEmpty(input$rnaClusteringCustomList))
					shinyjs::disable("performRnaClustering")
				else if (input$rnaClusterGeneList=="degenes" 
					&& is.null(currentRnaDeTable$totalTable))
					shinyjs::disable("performRnaClustering")
				else if (input$rnaClusteringVariablesRadio=="fc" 
					&& length(unique(currentMetadata$final$class))<=2)
					shinyjs::disable("performRnaClustering")
				else
					shinyjs::enable("performRnaClustering")
			}
        })
    })
}

