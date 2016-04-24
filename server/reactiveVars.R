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
        displayTables=list(),
        lengths=NULL
    )
    
    maPlots <- reactiveValues(
        #maPlot=data.frame(A=1,M=1,Status=1,
        #    Gene="Resulting MA plots will be displayed here"),
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
        maData=NULL,
        maColours=list(
            Up="#B40000",
            Down="#00B400",
            Neutral="#6B6B6B"
        ),
        maZoom=list(
            x=NULL,
            y=NULL
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
    
    currentHeatmap <- reactiveValues(
        entry=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label="Resulting heatmap will\nbe displayed here"),
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
        timeout=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label=paste("Clustering operation took too long\nto complete",
                    " andwas aborted.\nConsider lowering the number of genes",
                    "\nand/or conditions to prevent this.",sep="")),
                aes(x=x,y=y,label=label),size=9) +
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
        error=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label=paste("Clustering operation resulted in an error\nmost ",
                    "probably because of memory reasons.\nConsider lowering ",
                    "the number of genes\nand/or conditions to prevent this.",
                    sep="")),
                aes(x=x,y=y,label=label),size=9,color="red2") +
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
        data=NULL,
        opts=list(
            dendrogram="both",
            distfun="dist",
            hclustfun="hclust",
            Rowv=TRUE,
            Colv=TRUE,
            k_row=1,
            k_col=1,
            colors="RdYlBu"
        )
    )
    
    currentCorrelation <- reactiveValues(
        entryCor=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label="Resulting figures will\nbe displayed here"),
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
        errorCor=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label=paste("Correlation analysis has failed.\nThe most ",
                    "probable reason is that the\nselected gene set does not ",
                    "show\nenough variability to produce a\ncorrelation ",
                    "matrix. Select another gene\nset and try again.",sep="")),
                aes(x=x,y=y,label=label),size=9,color="red2") +
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
        entryMds=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label="MDS sample will\nbe displayed here"),
                aes(x=x,y=y,label=label),size=5) +
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
        warnMds=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label="MDS failed! Try\nchanging some parameters."),
                aes(x=x,y=y,label=label),size=5,color="orange") +
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
        errorMds=ggplot(data=data.frame(x=1:100,y=1:100)) + 
            geom_text(data=data.frame(x=50,y=50,
                label="Error creating MDS!\nSee error on the right"),
                aes(x=x,y=y,label=label),size=5,color="red2") +
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
        corMatrix=NULL,
        datMatrix=NULL,
        mdsRsq=NULL,
        refGene=NULL,
        what="samples",
        opts=list(
            method="pearson",
            symm=FALSE,
            colors=c("#FFFF00","#BEBEBE","#0000FF")
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
        currentMdsTables=currentMdsTables,
        currentHeatmap=currentHeatmap,
        currentCorrelation=currentCorrelation
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

clearReactiveVars <- function(allReactiveVars) {
    allReactiveVars$currentGenes$genes <- NULL
    allReactiveVars$currentGenes$coords <- list(
        chromosome=NULL,
        start=NULL,
        end=NULL,
        strand=NULL,
        name=NULL
    )
    
    allReactiveVars$customRegions$chromosome <- NULL
    allReactiveVars$customRegions$start <- NULL
    allReactiveVars$customRegions$end <- NULL
    allReactiveVars$customRegions$strand <- NULL
    allReactiveVars$customRegions$name <- NULL
    
    allReactiveVars$genePlots$geneProfile <- 
        ggplot(data=data.frame(x=1:100,y=1:100)) + 
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
            )
    allReactiveVars$genePlots$rendered <- TRUE
    
    allReactiveVars$customArea$chromosome <- NULL
    allReactiveVars$customArea$start <- NULL
    allReactiveVars$customArea$end <- NULL
    allReactiveVars$customArea$strand <- NULL
    allReactiveVars$customArea$name <- NULL
    
    allReactiveVars$customArea$customRegionsInArea <- NULL
    allReactiveVars$customArea$chromosome <- NULL
    allReactiveVars$customArea$start <- NULL
    allReactiveVars$customArea$end <- NULL
    allReactiveVars$customArea$strand <- NULL
    allReactiveVars$customArea$name <- NULL
    
    allReactiveVars$areaPlots$areaProfile <-
        ggplot(data=data.frame(x=1:100,y=1:100)) + 
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
            )
    allReactiveVars$areaPlots$rendered=FALSE
    
    allReactiveVars$currentTables <- reactiveValues()
    
    allReactiveVars$customRnaRegions$chromosome <- NULL
    allReactiveVars$customRnaRegions$start <- NULL
    allReactiveVars$customRnaRegions$end <- NULL
    allReactiveVars$customRnaRegions$strand <- NULL
    allReactiveVars$customRnaRegions$name <- NULL
    
    allReactiveVars$currentCustomRnaTables$tables <- list()
    allReactiveVars$currentCustomRnaTables$lengths <- NULL
    
    allReactiveVars$maPlots$maPlot <- 
        ggplot(data=data.frame(x=1:100,y=1:100)) + 
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
            )
    allReactiveVars$maPlots$maData <- NULL
    allReactiveVars$maPlots$maColours <- list(
        Up <- "#B40000",
        Down <- "#00B400",
        Neutral <- "#6B6B6B"
    )
    allReactiveVars$maPlots$maZoom <- list(
        x <- NULL,
        y <- NULL
    )
    allReactiveVars$maPlots$rendered <- TRUE
    
    allReactiveVars$currentPipelineOutput$annotation <- NULL
    allReactiveVars$currentPipelineOutput$counts <- NULL
    allReactiveVars$currentPipelineOutput$flags <- NULL
    allReactiveVars$currentPipelineOutput$classList <- NULL
    allReactiveVars$currentPipelineOutput$contrastList <- NULL
    allReactiveVars$currentPipelineOutput$pValue <- NULL
    allReactiveVars$currentPipelineOutput$fdr <- NULL
    
    allReactiveVars$currentRnaDeTable$totalTable <- NULL
    allReactiveVars$currentRnaDeTable$tableFilters <- list(
        p=0.05,
        fdr=0.05,
        scale="natural",
        fc=c(0.5,2),
        bt=NULL,
        genes=NULL,
        chr=NULL
    )
    
    allReactiveVars$currentHeatmap$data <- NULL
    allReactiveVars$currentHeatmap$opts <- list(
        dendrogram="both",
        distfun="dist",
        hclustfun="hclust",
        Rowv=TRUE,
        Colv=TRUE,
        k_row=1,
        k_col=1,
        colors="RdYlBu"
    )
        
    allReactiveVars$currentCorrelation$corMatrix <- NULL
    allReactiveVars$currentCorrelation$datMatrix <- NULL
    allReactiveVars$currentCorrelation$what <- "samples"
    allReactiveVars$currentCorrelation$opts <- list(
        method="pearson",
        symm=FALSE,
        colors=c("#FFFF00","#BEBEBE","#0000FF")
    )
        
    allReactiveVars$currentMdsTables <- reactiveValues()
    
    return(allReactiveVars)
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
