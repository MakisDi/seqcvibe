diffExprTabPanelEventReactive <- function(input,output,session,
	allReactiveVars,allReactiveMsgs) {
	currentMetadata <- allReactiveVars$currentMetadata
	
	return(NULL)
}

diffExprTabPanelReactive <- function(input,output,session,
	allReactiveVars,allReactiveMsgs) {
	currentMetadata <- allReactiveVars$currentMetadata
	
	return(NULL)
}

diffExprTabPanelRenderUI <- function(output,session,allReactiveVars,
	allReactiveMsgs) {
	currentMetadata <- allReactiveVars$currentMetadata
	maPlots <- allReactiveVars$maPlots
	
	output$checkboxBiotypeListRna <- renderUI({
		bts <- getBiotypes(currentMetadata$genome)
		lapply(bts,function(b) {
			checkboxInput(
				inputId=b,
				label=b,
				value=FALSE
			)
		})
	})
	
	output$setDeChrs <- renderUI({
		selectizeInput(
			inputId="customDeChr",
			label="Filter by chromosome", 
			choices=c("Show all",getValidChromosomes("hg19"))
		)
	})
	
	output$checkboxBiotypeListAnalyzedRna <- renderUI({
		bts <- getBiotypes(currentMetadata$genome)
		lapply(bts,function(b) {
			checkboxInput(
				inputId=b,
				label=b,
				value=FALSE
			)
		})
	})
	
	output$rnaDeMAPlot <- renderPlot({
		maPlots$maPlot
	})
}

diffExprTabPanelObserve <- function(input,output,session,
	allReactiveVars,allReactiveMsgs) {
	currentMetadata <- allReactiveVars$currentMetadata
	
	diffExprTabPanelReactiveEvents <- 
		dataSelectorTabPanelEventReactive(input,output,session,
			allReactiveVars,allReactiveMsgs)
	
	diffExprTabPanelReactiveExprs <- 
		dataSelectorTabPanelReactive(input,output,session,
			allReactiveVars,allReactiveMsgs)
	
	diffExprTabPanelRenderUI(output,session,allReactiveVars,
		allReactiveMsgs)
	
	observe({
	})
}

################################################################################

