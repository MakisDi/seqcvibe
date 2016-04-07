initReactiveVars <- function() {
    currentMetadata <- reactiveValues(
        source=sources[1],
        dataset=datasets[1],
        class=classes,
        genome=genomes[1],
        metadata=NULL,
        final=NULL
    )
    
    currentGenes <- reactiveValues(
        genes=NULL,
        coords=list(
            chromosome=NULL,
            start=NULL,
            end=NULL,
            strand=NULL,
            name=NULL
        )
    )
    
    customRegions <- reactiveValues(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    currentOpts <- reactiveValues(
        flank=c(2000,2000),
        sumStat="mean",
        trim=0.1,
        colours=NULL
    )
    
    genePlots <- reactiveValues(
        geneProfile=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label="Gene profiles will\nbe displayed here"),
                aes(x=x,y=y,label=label),size=10) +
            theme(
                axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank()
            ),
        rendered=TRUE
    )
    
    customArea <- reactiveValues(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    customRegionsInArea <- reactiveValues(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    areaPlots <- reactiveValues(
        areaProfile=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label="Area profiles will\nbe displayed here"),
                aes(x=x,y=y,label=label),size=10) +
            theme(
                axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank()
            ),
        rendered=FALSE
    )
    
    currentTables <- reactiveValues()
    
    customRnaRegions <- reactiveValues(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    currentCustomRnaTables <- reactiveValues(
        tables=list(),
        lengths=NULL
    )
    
    maPlots <- reactiveValues(
        maPlot=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label="Resulting MA plots will\nbe displayed here"),
                aes(x=x,y=y,label=label),size=10) +
            theme(
                axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank()
            ),
        rendered=TRUE
    )
    
    currentPipelineOutput <- reactiveValues(
        annotation=NULL,
        counts=NULL,
        flags=NULL,
        classList=NULL,
        contrastList=NULL,
        pValue=NULL,
        fdr=NULL
    )
    
    currentRnaDeTable <- reactiveValues(
        totalTable=NULL,
        tableFilters=list(
            p=0.05,
            fdr=0.05,
            scale="natural",
            fc=c(0.5,2),
            bt=NULL,
            genes=NULL,
            chr=NULL
        )
    )
    
    currentMdsTables <- reactiveValues()
    
    return(list(
        currentMetadata=currentMetadata,
        currentGenes=currentGenes,
        customRegions=customRegions,
        currentOpts=currentOpts,
        genePlots=genePlots,
        customArea=customArea,
        customRegionsInArea=customRegionsInArea,
        areaPlots=areaPlots,
        currentTables=currentTables,
        customRnaRegions=customRnaRegions,
        currentCustomRnaTables=currentCustomRnaTables,
        maPlots=maPlots,
        currentPipelineOutput=currentPipelineOutput,
        currentRnaDeTable=currentRnaDeTable,
        currentMdsTables=currentMdsTables
    ))
}

initReactiveMsgs <- function() {
    dataSelectorMessages <- reactiveValues(
        messages=list(
            list(
                type="INFO",
                msg=paste(getTime("INFO"),"Welcome to the data selector of ",
                    "BigSeqCVis! This is an info message. Make your ",
                    "selections on the left.")
            )
        )
    )
    
    geneExplorerMessages <- reactiveValues(
        messages=list(
            list(
                type="INFO",
                msg=paste(getTime("INFO"),"Welcome to the gene selector of ",
                    "BigSeqCVis! This is an info message. Make your ",
                    "selections on the left.")
            )
        )
    )
    
    areaExplorerMessages <- reactiveValues(
        messages=list(
            list(
                type="INFO",
                msg=paste(getTime("INFO"),"Welcome to the area explorer ",
                    "of BigSeqCVis! This is an info message. Make your ",
                    "selections on the left.")
            )
        )
    )
    
    return(list(
        dataSelectorMessages=dataSelectorMessages,
        geneExplorerMessages=geneExplorerMessages,
        areaExplorerMessages=areaExplorerMessages
    ))
}

## Messages boilerplate
#dataSelectorMessages <- reactiveValues(
#   messages=list(
#       list(
#           type="INFO",
#           msg=paste(getTime("INFO"),"Example info message.")
#       ),
#       list(
#           type="SUCCESS",
#           msg=paste(getTime("SUCCESS"),"Example success message.")
#       ),
#       list(
#           type="WARNING",
#           msg=paste(getTime("WARNING"),"Example warning message.")
#       ),
#       list(
#           type="ERROR",
#           msg=paste(getTime("ERROR"),"Example error message.")
#       )
#   )
#)
