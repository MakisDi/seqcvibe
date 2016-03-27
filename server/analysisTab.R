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

