mdsPcaTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentPipelineOutput <- allReactiveVars$currentPipelineOutput
    currentCustomRnaTables <- allReactiveVars$currentCustomRnaTables
    currentDimRed <- allReactiveVars$currentDimRed
        
    performRnaMdsPca <- eventReactive(input$performRnaMdsPca,{
        if (is.null(currentMetadata$final)) {
            output$rnaMdsPcaSettingsError <- renderUI({
                div(class="error-message",paste("You must create a ",
                    "dataset first!",sep=""))
            })
            return()
        }
        output$rnaMdsPcaSettingsError <- renderUI({div()})
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        meta <- currentMetadata$final
        samples <- as.character(meta$sample_id)
        D <- loadedData[[s]][[d]]
        genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
        
        # Determine which genes to use
        switch(input$rnaMdsPcaGeneList,
            select = {
                g <- input$selectMdsPcaGeneName
                if (input$rnaCorrelateWhat=="refgene" 
                    && !isEmpty(input$rnaMdsPcaRefGene))
                    g <- unique(c(input$rnaMdsPcaRefGene,g))
            },
            custom = {
                g <- input$rnaMdsPcaCustomList
                g <- strsplit(g,split="\n")[[1]]
                genes <- loadedGenomes[[currentMetadata$genome]]$geneNames
                m <- match(g,names(genes))
                na <- which(is.na(m))
                if (length(na)>0)
                    m <- m[-na]
                g <- genes[m]
                if (input$rnaCorrelateWhat=="refgene" 
                    && !isEmpty(input$rnaMdsPcaRefGene))
                    g <- unique(c(input$rnaMdsPcaRefGene,g))
            },
            expr = {
                if (is.null(currentPipelineOutput$counts)) {
                    output$rnaMdsPcaSettingsError <- renderUI({
                        div(class="error-message",paste("You must create a ",
                            "dataset first!",sep=""))
                    })
                    return()
                }
                else {
                    output$rnaMdsPcaSettingsError <- renderUI({div()})
                    # Code to select the genes
                }
            },
            all = {
                bad <- apply(D$norm,1,function(x) { return(all(x==0)) })
                g <- genes[-which(bad)]
                if (input$rnaCorrelateWhat=="refgene" 
                    && !isEmpty(input$rnaMdsPcaRefGene)) {
                        if (!(input$rnaMdsPcaRefGene %in% g))
                            g <- unique(c(input$rnaMdsPcaRefGene,g))
                    }
            }
        )
        
        # Determine the measurements table
        switch(input$rnaMdsPcaMeasureRadio,
            counts = {
                tab <- D$norm[g,samples,drop=FALSE]
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    tab <- rbind(tab,A)
                }
            },
            rpkm = {
                tab <- round(edgeR::rpkm(
                    D$counts[g,samples,drop=FALSE],
                    gene.length=D$length[g],
                    lib.size=unlist(D$libsize[samples])
                    ),digits=6)
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    A <- round(edgeR::rpkm(A,
                        gene.length=currentCustomRnaTables$lengths,
                        lib.size=unlist(D$libsize[samples])
                    ),digits=6)
                    tab <- rbind(tab,A)
                }
            },
            rpgm = {
                tab <- round(D$norm[g,
                    samples,drop=FALSE]/D$length[g],
                        digits=6)
                if (!is.null(currentCustomRnaTables$lengths)) {
                    A <- do.call("cbind",currentCustomRnaTables$tables)
                    A <- A[,colnames(tab),drop=FALSE]
                    A <- round(A/currentCustomRnaTables$lengths,
                        digits=6)
                    tab <- rbind(tab,A)
                }
            }
        )
        # We now have the tab, transformations based on other selections
        switch(input$rnaMdsPcaScaleRadio,
            natural = {
                tab <- tab
            },
            log2 = {
                tab <- round(log2(tab+1),digits=6)
            }
        )
        if (!is.null(currentMetadata$final$alt_id))
            colnames(tab) <- as.character(currentMetadata$final$alt_id)
        
        mds.obj <- pca.obj <- d <- NULL
        switch(input$rnaDimRedMethod,
            mds = {
                distFun <- distFuns()
                dfun=distFun[[input$rnaMdsDistMethod]]
                tryCatch({
					d <- dfun(t(tab))
					mds.obj <- cmdscale(d,eig=TRUE,
						k=as.numeric(input$rnaMdsKDim))
				},
				warning=function(w) {
					currentDimRed$mdsPlot <- 
						ggmessage(paste("An unexpected warning occured:\n",
							"Try changing your settings."),type="warning")
				},
				error=function(e) {
					currentDimRed$mdsPlot <- 
						ggmessage(paste("An unexpected error occured:\n",
							"Try changing your settings."),type="error")
				},
				finally="")
            },
            pca = {
            }
        )
        
        currentDimRed$datMatrix <- tab
        currentDimRed$mdsObj <- mds.obj
        currentDimRed$pcaObj <- pca.obj
        currentDimRed$mdsGof$dist <- d
        currentDimRed$opts$method <- input$rnaMdsPcaMethod
    })
    
    return(list(
        performRnaMdsPca=performRnaMdsPca
    ))
}

mdsPcaTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentDimRed <- allReactiveVars$currentDimRed
    
    updateMdsPcaPlotColours <- reactive({
        c <- currentMetadata$class
        lapply(c,function(x) {
            observeEvent(input[[paste("rnaMdsPcaColour_",x,sep="")]],{
                newc <- input[[paste("rnaMdsPcaColour_",x,sep="")]]
                if (is.null(currentDimRed$opts$colors)) {
                    currentDimRed$opts$colors <- character(length(c))
                    names(currentDimRed$opts$colors) <- c
                    currentDimRed$opts$colors[x] <- newc
                }
                else {
                    if (newc!=currentDimRed$opts$colors[x])
                        currentDimRed$opts$colors[x] <- newc
                }
            })
        })
    })
    
    updateMdsPlot <- reactive({
        if (is.null(currentDimRed$mdsObj))
            currentDimRed$mdsPlot <-
                ggmessage("Resulting MDS plots will\nbe displayed here")
        else {
            mds.obj <- currentDimRed$mdsObj
            cc <- unique(as.character(currentMetadata$final$class))
            classList <- vector("list",length(cc))
            names(classList) <- cc
            for (cl in cc)
                classList[[cl]] <- 
                    as.character(
                        currentMetadata$final$sample_id[which(
                            currentMetadata$final$class==cl)])
            classes <- as.factor(rep(names(classList),lengths(classList)))
            i <- as.numeric(input$rnaDimRedXAxis)
            j <- as.numeric(input$rnaDimRedYAxis)
            if (!is.na(i) && !is.na(j)) {
                nd <- as.dist(0.5*(1-cor(t(mds.obj$points))))
                currentDimRed$mdsGof$rsq <- 
                    cor(c(currentDimRed$mdsGof$dist),c(nd))^2
                currentDimRed$mdsGof$gof <- mds.obj$GOF
                for.ggplot <- data.frame(
                    x=mds.obj$points[,i],
                    y=mds.obj$points[,j],
                    Condition=classes
                )
                rownames(for.ggplot) <- rownames(mds.obj$points)
                classColors <- currentDimRed$opts$colors
                names(classColors) <- cc
                if (nrow(for.ggplot)<=100) {
                    psize <- 3
                    tsize <- 4
                }
                else {
                    psize <- 2
                    tsize <- 3
                }
                mds <- ggplot() +
                    geom_point(data=for.ggplot,mapping=aes(x=x,y=y,
                        colour=Condition),size=psize) +
                    xlab(paste("\nPrincipal Coordinate",i)) +
                    ylab(paste("Principal Coordinate",j)) +
                    theme_bw() +
                    theme(
                        axis.title.x=element_text(size=14),
                        axis.title.y=element_text(size=14),
                        axis.text.x=element_text(size=12,face="bold"),
                        axis.text.y=element_text(size=12,face="bold"),
                        legend.position="bottom",
                        legend.text=element_text(size=14),
                        legend.key=element_blank()
                    ) +
                    scale_color_manual(values=classColors) +
                    scale_fill_manual(values=classColors)
                if (input$rnaMdsPcaTogglePointNames) {
                    require(ggrepel)
                    mds <- mds +
                        geom_text_repel(data=for.ggplot,mapping=aes(x=x,y=y,
                            label=rownames(for.ggplot)),size=tsize)
                }
                currentDimRed$mdsPlot <- mds
            }
        }           
    })
    
    handleRnaMdsPcaSelection <- reactive({
        observeEvent(input$clearRnaMdsPcaSelection,{
            proxy <- dataTableProxy("rnaMdsPcaGeneSampleTable")
            selectRows(proxy,NULL)
        })
        
        observeEvent(input$invertRnaMdsPcaSelection,{
            N <- input$rnaMdsPcaGeneSampleTable_rows_all
            sel <- input$rnaMdsPcaGeneSampleTable_rows_selected
            if (length(sel)>0) {
                N <- N[-sel]
                proxy <- dataTableProxy("rnaMdsPcaGeneSampleTable")
                selectRows(proxy,N)
            }
        })
    })
    
    handleRnaMdsPcaDownload <- reactive({
        output$exportRnaMdsPcaSelection <- 
            downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste("mdspca_",tt,".txt", sep='')
                },
                content=function(con) {
                    sel <- input$rnaMdsPcaGeneSampleTable_rows_selected
                    if (length(sel)>0) {
                        gNames <- as.character(loadedGenomes[[
                            currentMetadata$genome]]$dbGene[rownames(
                            currentDimRed$datMatrix)]$gene_name)
                        res <- data.frame(
                            gene_name=gNames,
                            currentDimRed$datMatrix
                        )
                        write.table(res[sel,],file=con,sep="\t",quote=FALSE,
                            row.names=FALSE)
                    }
                }
            )
            
        output$exportRnaMdsPcaAll <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("mdspca_",tt,".txt", sep='')
            },
            content=function(con) {
                gNames <- as.character(loadedGenomes[[
                    currentMetadata$genome]]$dbGene[rownames(
                    currentDimRed$datMatrix)]$gene_name)
                res <- data.frame(
                    gene_name=gNames,
                    currentDimRed$datMatrix
                )
                write.table(res,file=con,sep="\t",quote=FALSE,row.names=FALSE)
            }
        )
    })
    
    return(list(
        updateMdsPcaPlotColours=updateMdsPcaPlotColours,
        updateMdsPlot=updateMdsPlot,
        handleRnaMdsPcaSelection=handleRnaMdsPcaSelection,
        handleRnaMdsPcaDownload=handleRnaMdsPcaDownload
    ))
}

mdsPcaTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentDimRed <- allReactiveVars$currentDimRed
    
    output$rnaMdsPcaGenesSamples <- renderUI({
        if (is.null(currentMetadata$final))
            div(
                style="display:inline-block; margin:5px;",
                h4("Please create a dataset first.")
            )
        else if (is.null(currentDimRed$datMatrix)) {
            div(
                style="display:inline-block; margin:5px;",
                h4("Please run a variance projection analysis first.")
            )
        }
        else {
            output$rnaMdsPcaGeneSampleTable <- DT::renderDataTable(
                if (is.null(currentDimRed$datMatrix))
                    data.frame(
                        name=character(0),
                        value=numeric(0)
                    )
                else {
                    s <- currentMetadata$source
                    d <- currentMetadata$dataset
                    gNames <- as.character(loadedGenomes[[
                        currentMetadata$genome]]$dbGene[rownames(
                        currentDimRed$datMatrix)]$gene_name)
                    data.frame(
                        gene_name=gNames,
                        currentDimRed$datMatrix
                    )
                },
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
                    DT::dataTableOutput("rnaMdsPcaGeneSampleTable"),
                    br(),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="clearRnaMdsPcaSelection",
                            label="Clear selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        actionButton(
                            inputId="invertRnaMdsPcaSelection",
                            label="Invert selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaMdsPcaSelection",
                            label="Export selection",
                            class="btn-xs"
                        )
                    ),
                    div(
                        style="display:inline-block; margin:5px;",
                        downloadButton(
                            outputId="exportRnaMdsPcaAll",
                            label="Export all",
                            class="btn-xs"
                        )
                    )
                )
            )
        }
    })
    
    output$rnaMdsPcaPlotColours <- renderUI({
        if (!is.null(currentMetadata$final)) {
            c <- unique(as.character(currentMetadata$final$class))
            c(
                list(h4("Class colours")),
                lapply(1:length(c),function(i,c) {
                    colourInput(
                        inputId=paste("rnaMdsPcaColour_",c[i],sep=""),
                        label=paste("Select colour for",c[i]),
                        value=baseColours[i]
                    )
                },c)
            )
        }
    })
    
    output$rnaMdsPcaDisplay <- renderText({
        if (is.null(currentDimRed$mdsObj))
            paste("R-square goodness of fit: ","MDS goodness of fit: ",
                sep="   ")
        else
            paste(
                "R-square goodness of fit: ",
                paste(round(100*currentDimRed$mdsGof$rsq,3),"%",sep=""),
                "   ",
                "MDS goodness of fit: ",
                round(currentDimRed$mdsObj$GOF[1],5),
                round(currentDimRed$mdsObj$GOF[2],5)
            )
    })
    
    output$rnaMdsPlot <- renderPlot({
        currentDimRed$mdsPlot
    })
    
    output$rnaPcaScreePlot <- renderPlot({
        currentDimRed$pcaScreePlot
    })
    
    output$rnaPcaScoresPlot <- renderPlot({
        currentDimRed$pcaScoresPlot
    })
    
    output$rnaPcaLoadingsPlot <- renderPlot({
        currentDimRed$pcaLoadingsPlot
    })
    
    output$rnaPcaBiplotPlot <- renderPlot({
        currentDimRed$pcaBiplotPlot
    })
}

mdsPcaTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentDimRed <- allReactiveVars$currentDimRed
    
    mdsPcaTabPanelReactiveEvents <- 
        mdsPcaTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    performRnaMdsPca <- mdsPcaTabPanelReactiveEvents$performRnaMdsPca
    
    mdsPcaTabPanelReactiveExprs <- 
        mdsPcaTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    updateMdsPcaPlotColours <- 
        mdsPcaTabPanelReactiveExprs$updateMdsPcaPlotColours
    updateMdsPlot <- mdsPcaTabPanelReactiveExprs$updateMdsPlot
    handleRnaMdsPcaSelection <- 
        mdsPcaTabPanelReactiveExprs$handleRnaMdsPcaSelection
    handleRnaMdsPcaDownload <- 
        mdsPcaTabPanelReactiveExprs$handleRnaMdsPcaSelection
            
    mdsPcaTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
        
    observe({
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        if (is.null(currentMetadata$final))
            updateSelectizeInput(session,"selectMdsPcaGeneName",
                choices=NULL,
                server=TRUE
            )
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            g <- isolate({input$selectCorrelationGeneName})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"selectMdsPcaGeneName",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
    
    observe({
        if (isEmpty(currentMetadata$final)) {
            shinyjs::disable("performRnaMdsPca")
            updateSelectizeInput(session,"rnaMdsKDim",
                choices=NULL)
        }
        else {
            shinyjs::enable("performRnaMdsPca")
            updateSelectizeInput(session,"rnaMdsKDim",
                choices=2:(nrow(currentMetadata$final)-1))
        }
    })
    
    observe({
        if (isEmpty(currentDimRed$mdsObj)) {
            updateSelectizeInput(session,"rnaDimRedXAxis",
                choices=NULL)
            updateSelectizeInput(session,"rnaDimRedYAxis",
                choices=NULL)
        }
        else {
            choices <- as.character(1:as.numeric(input$rnaMdsKDim))
            names(choices) <- paste("PC",1:as.numeric(input$rnaMdsKDim))
            updateSelectizeInput(session,"rnaDimRedXAxis",
                choices=choices,selected="1")
            updateSelectizeInput(session,"rnaDimRedYAxis",
                choices=choices,selected="2")
        }
    })
    
    observe({
        updateMdsPcaPlotColours()
        updateMdsPlot()
        handleRnaMdsPcaSelection()
        handleRnaMdsPcaDownload()
    })
    
    observe({   
    })
    
    observe({
        tryCatch({
            shinyjs::disable("performRnaMdsPca")
            performRnaMdsPca()
        },error=function(e) {
            #print(e)
        },
        finally={
            shinyjs::enable("performRnaMdsPca")
        })
    })
}

