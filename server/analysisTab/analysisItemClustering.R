clusteringTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}

clusteringTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {

}

clusteringTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
	currentHeatmap <- allReactiveVars$currentHeatmap
	currentRnaDeTable <- allReactiveVars$currentRnaDeTable
		
	output$heatmapOutput <- renderUI({
		if (is.null(currentRnaDeTable$totalTable))
			output$heatmap <- renderPlot({
				currentHeatmap$heatmap
			})
		else
			output$heatmap <- renderD3heatmap({
				currentHeatmap$heatmap
			})
	})
}

clusteringTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    
    clusteringTabPanelReactiveEvents <- 
        clusteringTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
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
}

