clusteringTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
}

clusteringTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {

}

clusteringTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
}

clusteringTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    clusteringTabPanelReactiveEvents <- 
        clusteringTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    clusteringTabPanelReactiveExprs <- 
        clusteringTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    clusteringTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
}

