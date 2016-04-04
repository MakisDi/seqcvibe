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
        cc <- unique(as.character(meta$class))
        classList <- vector("list",length(cc))
        names(classList) <- cc
        for (cl in cc)
            classList[[cl]] <- 
                as.character(meta$sample_id[which(meta$class==cl)])
        control <- input$rnaDePipelineControl
        treatments <- setdiff(cc,control)
        contrast <- c(treatments,control)
        contrast <- paste(contrast,collapse="_vs_")
        
        # Set up normalization
        if (input$rnaDePipeline=="deseq")
            norm <- "deseq"
        else
            norm <- "edger"
        
        # Set up filters
        exprFilt <- list(median=TRUE,mean=FALSE,quantile=NA,known=NA)
        switch(input$rnaDeGeneFilter,
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
            btFilt <- lapply(bts,function(b) {
                return(input[[b]])
            })
            names(btFilt) <- bts
        }
        geneFilters=list(
            length=list(
                length=as.numeric(input$rnaDeGeneLengthFilterValue)
            ),
            avg.reads=list(
                average.per.bp=100,
                quantile=0.25
            ),
            expression=exprFilt, 
            biotype=btFilt
        )
        
        pipOutput <- tryCatch(
            metaseqr(
                counts=D,
                sample.list=classList,
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
        
        currentPipelineOutput$annotation <- pipOutput$complete$gene.data
        currentPipelineOutput$counts <- pipOutput$complete$norm.counts
        currentPipelineOutput$flags <- pipOutput$complete$flags
        currentPipelineOutput$classList <- pipOutput$complete$sample.list
        currentPipelineOutput$contrastList <- pipOutput$complete$contrast
        currentPipelineOutput$pValue <- pipOutput$complete$p.value
        currentPipelineOutput$fdr <- pipOutput$complete$fdr
    })
    
    return(list(
        runPipeline=runPipeline
    ))
}

diffExprTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentPipelineOutput <- allReactiveVars$currentPipelineOutput
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    
    rnaDeTotalTable <- reactive({
        s <- unique(as.character(currentMetadata$source))
        d <- unique(as.character(currentMetadata$dataset))
        
        ann <- currentPipelineOutput$annotation
        counts <- currentPipelineOutput$counts
        flags <- currentPipelineOutput$flags
        classList <- currentPipelineOutput$classList
        contrastList <- currentPipelineOutput$contrastList
        pValue <- currentPipelineOutput$pValue
        fdr <- currentPipelineOutput$fdr
        
        if (!is.null(counts)) {
            p <- pValue[[contrastList[1]]][,1,drop=FALSE]
            fdr <- fdr[[contrastList[1]]][,1,drop=FALSE]
            
            filters <- currentRnaDeTable$tableFilters
            filterInds <- list(
                stat=rownames(counts),
                fold=rownames(counts),
                chr=rownames(counts),
                bt=rownames(counts)
            )
            
            # Gradual application of filters for faster rendering
            # 1. Chromosomes
            tmp <- rownames(counts)
            if (!is.null(filters$chr))
                tmp <- 
                    names(which(as.character(ann$chromosome) %in% filters$chr))
            if (length(tmp)>0)
                filterInds$chr <- tmp
            # 2. Biotypes
            if (!is.null(filters$bt))
                tmp <- names(which(as.character(ann$biotype) %in% filters$bt))
            if (length(tmp)>0)
                filterInds$bt <- tmp
            
            # 3. Statistical score
            if (input$statThresholdType=="pvalue")
                tmp <- names(which(p[,1]<filters$p))
            else if (input$statThresholdType=="fdr")
                tmp <- names(which(fdr[,1]<filters$fdr))
            if (length(tmp)>0)
                filterInds$stat <- tmp
                
            # Gather so far
            sofar <- Reduce("intersect",filterInds[c("stat","chr","bt")])
            sumTable <- data.frame(
                pvalue=p[sofar,,drop=FALSE],
                fdr=fdr[sofar,,drop=FALSE]
            )
            names(sumTable) <- c("pvalue","fdr")
            ann <- ann[sofar,]
            flags <- flags[sofar,]
            
            # Proceed with count table procesing
            tab <- counts[sofar,,drop=FALSE]
            switch(input$rnaDeValueCompRadio,
                counts = {
                    tab <- tab
                },
                rpkm = {
                    len <- loadedData[[s]][[d]]$length[sofar]
                    libsize=unlist(loadedData[[s]][[d]]$libsize[colnames(tab)])
                    tab <- round(edgeR::rpkm(tab,gene.length=len,
                        lib.size=libsize),digits=6)
                },
                rpgm = {
                    len <- loadedData[[s]][[d]]$length[sofar]
                    for (j in 1:ncol(tab))
                        tab[,j] <- round(tab[,j]/len,digits=6)
                }
            )
            switch(input$rnaDeValueScaleRadio,
                natural = {
                    tab <- tab
                },
                log2 = {
                    tab <- round(log2(tab+1),digits=6)
                }
            )
            
            fcMat <- round(makeFoldChange(contrastList[1],classList,tab,
                input$rnaDeValueScaleRadio),6)
            tmp <- names(which(apply(fcMat,1,function(x,f) {
                return(any(x<=f[1] | x>=f[2]))
            },filters$fc)))
            if (length(tmp)>0)
                filterInds$fold <- tmp
            else
                filterInds$fold <- sofar
                
            tab <- tab[filterInds$fold,,drop=FALSE]
            fcMat <- fcMat[filterInds$fold,,drop=FALSE]
            
            avgMatrix <- do.call("cbind",lapply(classList,function(x,tab,s,v) {
                makeStat(x,tab,s,v)
            },tab,input$rnaDeValueAverageRadio,input$rnaDeValueCompRadio))
            colnames(avgMatrix) <- paste(names(classList),
                input$rnaDeValueAverageRadio)
            
            stdMatrix <- do.call("cbind",lapply(classList,function(x,tab,s,v) {
                makeStat(x,tab,s,v)
            },tab,input$rnaDeValueDeviationRadio,input$rnaDeValueCompRadio))
            colnames(stdMatrix) <- paste(names(classList),
                input$rnaDeValueDeviationRadio,input$rnaDeValueCompRadio)
            
            totalTable <- cbind(
                ann[filterInds$fold,],
                sumTable[filterInds$fold,],
                fcMat,
                as.data.frame(avgMatrix),
                as.data.frame(stdMatrix),
                as.data.frame(flags[filterInds$fold,])
            )
            
            currentRnaDeTable$totalTable <- totalTable
        }
    })
    
    handleRnaDeAnalysisSummarySelection <- reactive({
        observeEvent(input$clearRnaDeSummarySelection,{
            proxy <- dataTableProxy("rnaDeAnalysisSummaryTable")
            selectRows(proxy,NULL)
        })
        
        observeEvent(input$invertRnaDeSummarySelection,{
            N <- input$rnaDeAnalysisSummaryTable_rows_all
            sel <- input$rnaDeAnalysisSummaryTable_rows_selected
            if (length(sel)>0) {
                N <- N[-sel]
                proxy <- dataTableProxy("rnaDeAnalysisSummaryTable")
                selectRows(proxy,N)
            }
        })
    })
    
    handleRnaDeAnalysisAnnotationSelection <- reactive({
        observeEvent(input$clearRnaDeAnnotationSelection,{
            proxy <- dataTableProxy("rnaDeAnalysisAnnotationTable")
            selectRows(proxy,NULL)
        })
        
        observeEvent(input$invertRnaDeAnnotationSelection,{
            N <- input$rnaDeAnalysisAnnotationTable_rows_all
            sel <- input$rnaDeAnalysisAnnotationTable_rows_selected
            if (length(sel)>0) {
                N <- N[-sel]
                proxy <- dataTableProxy("rnaDeAnalysisAnnotationTable")
                selectRows(proxy,N)
            }
        })
    })
    
    handleRnaDeAnalysisFlagsSelection <- reactive({
        observeEvent(input$clearRnaDeFlagsSelection,{
            proxy <- dataTableProxy("rnaDeAnalysisFlagsTable")
            selectRows(proxy,NULL)
        })
        
        observeEvent(input$invertRnaDeFlagsSelection,{
            N <- input$rnaDeAnalysisFlagsTable_rows_all
            sel <- input$rnaDeAnalysisFlagsTable_rows_selected
            if (length(sel)>0) {
                N <- N[-sel]
                proxy <- dataTableProxy("rnaDeAnalysisFlagsTable")
                selectRows(proxy,N)
            }
        })
    })
    
    handleRnaDeAnalysisAllSelection <- reactive({
        observeEvent(input$clearRnaDeAllSelection,{
            proxy <- dataTableProxy("rnaDeAnalysisAllTable")
            selectRows(proxy,NULL)
        })
        
        observeEvent(input$invertRnaDeAllSelection,{
            N <- input$rnaDeAnalysisAllTable_rows_all
            sel <- input$rnaDeAnalysisAllTable_rows_selected
            if (length(sel)>0) {
                N <- N[-sel]
                proxy <- dataTableProxy("rnaDeAnalysisAllTable")
                selectRows(proxy,N)
            }
        })
    })
    
    handleRnaDeAnalysisSummaryDownload <- reactive({
        output$exportRnaDeSummarySelection <- 
            downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste("diffexpr_summary_",tt,".txt", sep='')
                },
                content=function(con) {
                    sel <- input$rnaDeAnalysisSummaryTable_rows_selected
                    if (length(sel)>0) {
                        totalTable <- currentRnaDeTable$totalTable
                        res <- totalTable[,c("gene_name","pvalue","fdr")]
                        fcInd <- grep("_vs_",names(totalTable))
                        avgInd <- grep("mean|median",colnames(totalTable),
                            perl=TRUE)
                        devInd <- grep("sd|mad|IQR",colnames(totalTable),
                            perl=TRUE)
                        res <- cbind(res,totalTable[,c(fcInd,avgInd,devInd)]) 
                        write.table(res[sel,],file=con,sep="\t",quote=FALSE,
                            row.names=FALSE)
                    }
                }
            )
            
        output$exportRnaDeSummaryAll <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("diffexpr_summary_",tt,".txt", sep='')
            },
            content=function(con) {
                totalTable <- currentRnaDeTable$totalTable
                res <- totalTable[,c("gene_name","pvalue","fdr")]
                fcInd <- grep("_vs_",names(totalTable))
                avgInd <- grep("mean|median",colnames(totalTable),
                    perl=TRUE)
                devInd <- grep("sd|mad|IQR",colnames(totalTable),
                    perl=TRUE)
                res <- cbind(res,totalTable[,c(fcInd,avgInd,devInd)])
                write.table(res,file=con,sep="\t",quote=FALSE,row.names=FALSE)
            }
        )
    })
    
    handleRnaDeAnalysisAnnotationDownload <- reactive({
        output$exportRnaDeAnnotationSelection <- 
            downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste("diffexpr_annotation_",tt,".txt", sep='')
                },
                content=function(con) {
                    sel <- input$rnaDeAnalysisAnnotationTable_rows_selected
                    if (length(sel)>0) {
                        totalTable <- currentRnaDeTable$totalTable
                        res <- totalTable[,c("chromosome","start","end",
                            "gene_id","strand","gene_name","biotype")]
                        write.table(res[sel,],file=con,sep="\t",quote=FALSE,
                            row.names=FALSE)
                    }
                }
            )
            
        output$exportRnaDeAnnotationAll <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("diffexpr_annotation_",tt,".txt", sep='')
            },
            content=function(con) {
                totalTable <- currentRnaDeTable$totalTable
                res <- totalTable[,c("chromosome","start","end","gene_id",
                    "strand","gene_name","biotype")]
                write.table(res,file=con,sep="\t",quote=FALSE,row.names=FALSE)
            }
        )
    })
    
    handleRnaDeAnalysisFlagsDownload <- reactive({
        output$exportRnaDeFlagsSelection <- 
            downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste("diffexpr_flags_",tt,".txt", sep='')
                },
                content=function(con) {
                    sel <- input$rnaDeAnalysisFlagsTable_rows_selected
                    if (length(sel)>0) {
                        totalTable <- currentRnaDeTable$totalTable
                        res <- totalTable[,c("gene_name","LN","MD","MN","QN",
                            "KN","BT")] 
                        write.table(res[sel,],file=con,sep="\t",quote=FALSE,
                            row.names=FALSE)
                    }
                }
            )
            
        output$exportRnaDeFlagsAll <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("diffexpr_flags_",tt,".txt", sep='')
            },
            content=function(con) {
                totalTable <- currentRnaDeTable$totalTable
                res <- totalTable[,c("gene_name","LN","MD","MN","QN","KN","BT")]
                write.table(res,file=con,sep="\t",quote=FALSE,row.names=FALSE)
            }
        )
    })
    
    handleRnaDeAnalysisAllDownload <- reactive({
        output$exportRnaDeAllSelection <- 
            downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste("diffexpr_summary_",tt,".txt", sep='')
                },
                content=function(con) {
                    sel <- input$rnaDeAnalysisAllTable_rows_selected
                    if (length(sel)>0) {
                        totalTable <- currentRnaDeTable$totalTable
                        write.table(totalTable[sel,],file=con,sep="\t",
                            quote=FALSE,row.names=FALSE)
                    }
                }
            )
            
        output$exportRnaDeAllAll <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("diffexpr_summary_",tt,".txt", sep='')
            },
            content=function(con) {
                totalTable <- currentRnaDeTable$totalTable
                write.table(totalTable,file=con,sep="\t",quote=FALSE,
                    row.names=FALSE)
            }
        )
    })
    
    # Now the sliders
    statSliderUpdate <- reactive({
        if (input$statThresholdType=="pvalue") {
            index <- as.integer(input$pvalue)
            # Dirty hack for the 1st loading...
            if (index==0) index <- 11
            currentRnaDeTable$tableFilters$p <- statScoreValues[index]
        }
        else if (input$statThresholdType=="fdr") {
            index <- as.integer(input$fdr)
            currentRnaDeTable$tableFilters$fdr <- statScoreValues[index]
        }
    })
    
    foldChangeSliderUpdate <- reactive({
        if (input$foldThresholdType=="natural") {
            currentRnaDeTable$tableFilters$scale <- "natural"
            currentRnaDeTable$tableFilters$fc <- input$fcNatural
        }
        else if (input$foldThresholdType=="log2") {
            currentRnaDeTable$tableFilters$scale <- "log2"
            currentRnaDeTable$tableFilters$fc <- input$fcLog
        }
    })
    
    
    return(list(
        rnaDeTotalTable=rnaDeTotalTable,
        handleRnaDeAnalysisSummarySelection=handleRnaDeAnalysisSummarySelection,
        handleRnaDeAnalysisAnnotationSelection=
            handleRnaDeAnalysisAnnotationSelection,
        handleRnaDeAnalysisFlagsSelection=handleRnaDeAnalysisFlagsSelection,
        handleRnaDeAnalysisAllSelection=handleRnaDeAnalysisAllSelection,
        handleRnaDeAnalysisSummaryDownload=handleRnaDeAnalysisSummaryDownload,
        handleRnaDeAnalysisAnnotationDownload=
            handleRnaDeAnalysisAnnotationDownload,
        handleRnaDeAnalysisFlagsDownload=handleRnaDeAnalysisFlagsDownload,
        handleRnaDeAnalysisAllDownload=handleRnaDeAnalysisAllDownload,
        statSliderUpdate=statSliderUpdate,
        foldChangeSliderUpdate=foldChangeSliderUpdate
    ))
}

diffExprTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentRnaDeTable <- allReactiveVars$currentRnaDeTable
    maPlots <- allReactiveVars$maPlots
     
    output$rnaDeAnalysisSummary <- renderUI({
        if (is.null(currentRnaDeTable$totalTable))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please run a differential expression analysis ",
                    "first."))
            )
        else {
            totalTable <- currentRnaDeTable$totalTable
            res <- totalTable[,c("gene_name","pvalue","fdr")]
            fcInd <- grep("_vs_",names(totalTable))
            avgInd <- grep("mean|median",colnames(totalTable),perl=TRUE)
            devInd <- grep("sd|mad|IQR",colnames(totalTable),perl=TRUE)
            res <- cbind(res,totalTable[,c(fcInd,avgInd,devInd)])
            output$rnaDeAnalysisSummaryTable <- 
                DT::renderDataTable(
                    res,
                    class="display compact",
                    rownames=FALSE,
                    options=list(
                        searchHighlight=TRUE,
                        pageLength=10,
                        lengthMenu=c(10,20,50,100)
                    )
                )
            list(
                div(
                    class="small table-container",
                    DT::dataTableOutput("rnaDeAnalysisSummaryTable"),
                    br(),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="clearRnaDeSummarySelection",
                            label="Clear selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="invertRnaDeSummarySelection",
                            label="Invert selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeSummarySelection",
                            label="Export selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeSummaryAll",
                            label="Export all",
                            class="btn-xs"
                        )
                    )
                )
            )
        }
    })
    
    output$rnaDeAnalysisAnnotation <- renderUI({
        if (is.null(currentRnaDeTable$totalTable))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please run a differential expression analysis ",
                    "first."))
            )
        else {
            totalTable <- currentRnaDeTable$totalTable
            res <- totalTable[,c("chromosome","start","end","gene_id",
                "strand","gene_name","biotype")]
            output$rnaDeAnalysisAnnotationTable <- 
                DT::renderDataTable(
                    res,
                    class="display compact",
                    rownames=FALSE,
                    options=list(
                        searchHighlight=TRUE,
                        pageLength=10,
                        lengthMenu=c(10,20,50,100)
                    )
                )
            list(
                div(
                    class="small table-container",
                    DT::dataTableOutput("rnaDeAnalysisAnnotationTable"),
                    br(),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="clearRnaDeAnnotationSelection",
                            label="Clear selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="invertRnaDeAnnotationSelection",
                            label="Invert selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeAnnotationSelection",
                            label="Export selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeAnnotationAll",
                            label="Export all",
                            class="btn-xs"
                        )
                    )
                )
            )
        }
    })
    
    output$rnaDeAnalysisFlags <- renderUI({
        if (is.null(currentRnaDeTable$totalTable))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please run a differential expression analysis ",
                    "first."))
            )
        else {
            totalTable <- currentRnaDeTable$totalTable
            res <- totalTable[,c("gene_name","LN","MD","MN","QN","KN","BT")]
            output$rnaDeAnalysisFlagsTable <- 
                DT::renderDataTable(
                    res,
                    class="display compact",
                    rownames=FALSE,
                    options=list(
                        searchHighlight=TRUE,
                        pageLength=10,
                        lengthMenu=c(10,20,50,100)
                    )
                )
            list(
                div(
                    class="small table-container",
                    DT::dataTableOutput("rnaDeAnalysisFlagsTable"),
                    br(),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="clearRnaDeFlagsSelection",
                            label="Clear selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="invertRnaDeFlagsSelection",
                            label="Invert selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeFlagsSelection",
                            label="Export selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeFlagsAll",
                            label="Export all",
                            class="btn-xs"
                        )
                    )
                )
            )
        }
    })
    
    output$rnaDeAnalysisAll <- renderUI({
        if (is.null(currentRnaDeTable$totalTable))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please run a differential expression analysis ",
                    "first."))
            )
        else {
            totalTable <- currentRnaDeTable$totalTable
            res <- totalTable[,c("gene_name","pvalue","fdr")]
            fcInd <- grep("_vs_",names(totalTable))
            avgInd <- grep("mean|median",colnames(totalTable),perl=TRUE)
            devInd <- grep("sd|mad|IQR",colnames(totalTable),perl=TRUE)
            res <- cbind(res,totalTable[,c(fcInd,avgInd,devInd)])
            output$rnaDeAnalysisAllTable <-
                DT::renderDataTable(
                    res,
                    class="display compact",
                    rownames=FALSE,
                    options=list(
                        searchHighlight=TRUE,
                        pageLength=10,
                        lengthMenu=c(10,20,50,100)
                    )
                )
            list(
                div(
                    class="small table-container",
                    DT::dataTableOutput("rnaDeAnalysisAllTable"),
                    br(),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="clearRnaDeAllSelection",
                            label="Clear selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="invertRnaDeAllSelection",
                            label="Invert selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeAllSelection",
                            label="Export selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaDeAllAll",
                            label="Export all",
                            class="btn-xs"
                        )
                    )
                )
            )
        }
    })
    
    output$rnaDePipelineControl <- renderUI({
        if (is.null(currentMetadata$final))
            list(
                div(style="font-weight:bold","Select control condition"),
                helpText("Please create a dataset first from the 'Data ",
                    "selector' menu on the top")
            )
        else {
            cls <- unique(as.character(currentMetadata$final$class))
            selectInput(
                inputId="rnaDePipelineControl",
                label="Select control condition",
                choices=cls
            )
        }
    })
    
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
                inputId=paste(b,"asFilter",sep="_"),
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
        diffExprTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    runPipeline <- diffExprTabPanelReactiveEvents$runPipeline
    
    diffExprTabPanelReactiveExprs <- 
        diffExprTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    rnaDeTotalTable <- diffExprTabPanelReactiveExprs$rnaDeTotalTable
    handleRnaDeAnalysisSummarySelection <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisSummarySelection
    handleRnaDeAnalysisAnnotationSelection <- 
            diffExprTabPanelReactiveExprs$handleRnaDeAnalysisAnnotationSelection
    handleRnaDeAnalysisFlagsSelection <-
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisFlagsSelection
    handleRnaDeAnalysisAllSelection <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisAllSelection
    handleRnaDeAnalysisSummaryDownload <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisSummaryDownload
    handleRnaDeAnalysisAnnotationDownload <-
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisAnnotationDownload
    handleRnaDeAnalysisFlagsDownload <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisFlagsDownload
    handleRnaDeAnalysisAllDownload <- 
        diffExprTabPanelReactiveExprs$handleRnaDeAnalysisAllDownload
    statSliderUpdate <- diffExprTabPanelReactiveExprs$statSliderUpdate
    foldChangeSliderUpdate <- 
        diffExprTabPanelReactiveExprs$foldChangeSliderUpdate
    
    diffExprTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
    
    observe({
        rnaDeTotalTable()
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
        geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
        g <- isolate({input$rnaDeKnownFilter})
        i <- grep(paste0("^",g),geneNames,perl=TRUE)
        if (length(i)>0) {
            updateSelectizeInput(session,"rnaDeShowSpecificGenes",
                choices=geneNames[i],
                selected=g,
                server=TRUE
            )
        }
    })
    
    observe({
        if (isEmpty(currentMetadata$final))
            shinyjs::disable("performDeAnalysis")
        else
            shinyjs::enable("performDeAnalysis")
    })
    
    observe({
        handleRnaDeAnalysisSummarySelection()
        handleRnaDeAnalysisAnnotationSelection()
        handleRnaDeAnalysisFlagsSelection()
        handleRnaDeAnalysisAllSelection()
        handleRnaDeAnalysisSummaryDownload()
        handleRnaDeAnalysisAnnotationDownload()
        handleRnaDeAnalysisFlagsDownload()
        handleRnaDeAnalysisAllDownload()
        statSliderUpdate()
        foldChangeSliderUpdate()
    })
    
    observe({
        tryCatch({
            shinyjs::disable("performDeAnalysis")
            runPipeline()
        },error=function(e) {
            shinyjs::enable("performDeAnalysis")
            print(e)
        },
        finally={
            shinyjs::enable("performDeAnalysis")
        })
    })
}

################################################################################

