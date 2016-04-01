diffExprTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentPipelineOutput <- allReactiveVars$currentPipelineOutput
    
    runPipeline <- eventReactive(input$performDeAnalysis,{
        require(metaseqR)
            
        # Set up count data table with embedded annotation
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        meta <- currentMetadata$final
        samples <- as.character(meta$sample_id)
        dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
        ann <- data.frame(
            chromosome=as.character(seqnames(dbGene)),
            start=start(dbGene),
            end=end(dbGene),
            gene_id=as.character(dbGene$gene_id),
            gc_content=rep(0,length(dbGene)),
            strand=as.character(strand(dbGene)),
            gene_name=as.character(dbGene$gene_name),
            biotype=as.character(dbGene$biotype)
        )
        rownames(ann) <- ann$gene_id
        M <- loadedData[[s]][[d]]$counts[rownames(ann),samples]
        D <- cbind(ann,M)
        
        # Set up sample and contrast list for metaseqr
        cc <- unique(as.character(M$class))
        classList <- vector("list",length(cc))
        names(classList) <- cc
        for (cl in cc)
            classList[[cl]] <- 
                as.character(meta$sample_id[which(meta$class==cl)])
        control <- input$rnaDePipelineControl
        treatments <- setdiff(cc,control)
        contrast <- c(control,treatments)
        contrast <- paste(contrast,collapse="_vs_")
        
        # Set up normalization
        if (input$rnaDePipeline=="deseq")
            norm <- "deseq"
        else
            norm <- "edger"
        
        # Set up filters
        exprFilt <- list(median=TRUE,mean=FALSE,quantile=NA,known=NA)
        switch(rnaDeGeneFilter,
            median = {
                exprFilt <- list(
                    median=TRUE,
                    mean=FALSE,
                    quantile=NA,
                    known=NA
                )
            },
            mean = {
                exprFilt <- list(
                    median=FALSE,
                    mean=TRUE,
                    quantile=NA,
                    known=NA
                )
            },
            quantile = {
                qq <- as.numeric(input$rnaDeQuantileFilter)
                if (!is.na(qq))
                    exprFilt <- list(
                        median=FALSE,
                        mean=FALSE,
                        quantile=qq,
                        known=NA
                    )
            },
            known = {
                if (!isEmpty(input$rnaDeKnownFilter))
                    exprFilt <- list(
                        median=FALSE,
                        mean=FALSE,
                        quantile=NA,
                        known=input$rnaDeKnownFilter
                    )
            }
        )
        btFilt <- NULL
        if (input$rnaDeBiotypeFilter) {
            bts <- getBiotypes(currentMetadata$genome)
            btStatus <- lapply(bts,function(b) {
                return(output[[b]]$value)
            })
            names(btStatus) <- bts
        }
        geneFilters=list(
            length=list(
                length=as.numeric(input$rnaDeGeneLengthFilter)
            ),
            expression=exprFilt, 
            biotype=btFilt
        )
        
        pipOutput <- tryCatch(
            metaseqr(
                counts="loaded_matrix_with_annotation",
                contrast=contrast,
                annotation="embedded",
                id.col=4,
                gc.col=5,
                name.col=7,
                bt.col=8,
                org=currentMetadata$genome,
                count.type="gene",
                normalization=norm,
                statistics=input$rnaDePipeline,
                adjust.method=input$rnaDeMTC,
                fig.format="png",
                export.where=file.path(tempdir(),paste("analysis_",
                        format(Sys.time(),format="%Y%m%d%H%M%S"),sep="")),
                qc.plots="mds",
                exon.filters=NULL,
                gene.filters=geneFilters,
                when.apply.filter=input$rnaDeNormalizeWhen,
                export.what=c("annotation","p.value","adj.p.value","counts",
                    "flags"),
                export.scale="natural",
                export.values="normalized",
                export.stats="mean",
                save.gene.model=FALSE,
                report=FALSE,
                out.list=TRUE
            ),
        error=function(e) {
            print(e)
        },
        finally="")
        
        currentPipelineOutput$annotation <- pipOutput$gene.data
        currentPipelineOutput$counts <- pipOutput$norm.counts
        currentPipelineOutput$flags <- pipOutput$flags
        currentPipelineOutput$classList <- pipOutput$sample.list
        currentPipelineOutput$contastList <- pipOutput$contrast
        currentPipelineOutput$pValue <- pipOutput$p.value
        currentPipelineOutput$fdr <- pipOutput$fdr
    })
    
    
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
        if (!is.null(currentMetadata$final)) {
            shinyjs::enable("rnaDePipelineControl")
            cls <- unique(as.character(currentMetadata$final$class))
            updateSelectInput(session,"rnaDePipelineControl",choices=cls)
        }
    })
    
    observe({
        geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
        g <- isolate({input$rnaDeKnownFilter})
        i <- grep(paste0("^",g),geneNames,perl=TRUE)
        if (length(i)>0) {
            updateSelectizeInput(session,"rnaDeKnownFilter",
                choices=geneNames[i],
                selected=g,
                server=TRUE
            )
        }
    })
    
    observe({
    })
}

################################################################################

