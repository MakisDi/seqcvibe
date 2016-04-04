# server.R

# Load init script
source("config/init.R")

shinyServer(
    function(input,output,session) {
        source("server/reactiveVars.R",local=TRUE)
        source("server/dataSelectorTab.R",local=TRUE)
        source("server/signalViewerTab.R",local=TRUE)
        source("server/expressionViewerTab.R",local=TRUE)
        source("server/analysisTab.R",local=TRUE)
        source("server/genomeBrowserTab.R",local=TRUE)
        
        # Make %#^%$^%$@( globals visible AND changeable
        makeReactiveBinding("loadedGenomes")
        makeReactiveBinding("loadedData")
        
        # Initialize all the reactive variables used...
        allReactiveVars <- initReactiveVars()
        # ...and reactive messages
        allReactiveMsgs <- initReactiveMsgs()
        
        # Data selector
        dataSelectorTabPanelObserve(input,output,session,allReactiveVars,
            allReactiveMsgs)
        
        # Signal viewer - Gene signal
        geneSignalTabPanelObserve(input,output,session,allReactiveVars,
            allReactiveMsgs)
       
        # Signal viewer - Area signal
        areaSignalTabPanelObserve(input,output,session,allReactiveVars,
            allReactiveMsgs)
        
        # Expresion viewer - Known genes
        expressionExplorerTabPanelObserve(input,output,session,allReactiveVars,
            allReactiveMsgs)
        
        # Expresion viewer - Calculator
        expressionCalculatorTabPanelObserve(input,output,session,
            allReactiveVars,allReactiveMsgs)        
        
        # Analysis - Differential expression
        diffExprTabPanelObserve(input,output,session,allReactiveVars,
            allReactiveMsgs)
        
        # Genome browser
        genomeBrowserTabPanelObserve(input,output,session,allReactiveVars,
            allReactiveMsgs)
        
    }
)
