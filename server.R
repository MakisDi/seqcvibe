# server.R

# Load init script
source("config/init.R")
#source("server/reactiveVars.R")
#source("server/dataSelectorTab.R")

shinyServer(
    function(input,output,session) {
        source("server/reactiveVars.R",local=TRUE)
		source("server/dataSelectorTab.R",local=TRUE)
		source("server/signalViewerTab.R",local=TRUE)
		source("server/expressionViewerTab.R",local=TRUE)
		source("server/analysisTab.R",local=TRUE)
        
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
        output$genomeBrowser <- renderUI({
            tags$iframe(
                src=loadJBrowse(
                    source=as.character(currentMetadata$source),
                    dataset=as.character(currentMetadata$dataset),
                    config=currentMetadata$final,
                    org=currentMetadata$genome
                ),
                name="JBrowse",seamless=NA,
                height="800px",width="100%"
            )
        })
    }
)
