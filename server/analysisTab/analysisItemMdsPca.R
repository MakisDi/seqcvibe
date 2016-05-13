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
        #if (!is.null(currentMetadata$final$alt_id))
        #    colnames(tab) <- as.character(currentMetadata$final$alt_id)
        
        mds.obj <- pca.obj <- d <- NULL
        switch(input$rnaDimRedMethod,
            mds = {
                if (input$rnaMdsDistMethod %in% c("pearson","spearman")
                    && (nrow(tab)<2 || ncol(tab)<2)) {
                    output$rnaMdsPcaSettingsError <- renderUI({
                        div(class="error-message",paste("At least two genes ",
                            "and two samples are required for correlation ",
                            "similarity metrics!",sep=""))
                    })
                    currentDimRed$mdsPlot <- 
                        ggmessage(paste("An unexpected error occured.\n",
                            "Try changing your settings."),type="error")
                    return()
                }
                else {
                    output$rnaMdsPcaSettingsError <- renderUI({div()})
                    distFun <- distFuns()
                    dfun=distFun[[input$rnaMdsDistMethod]]
                    tryCatch({
                        d <- dfun(t(tab))
                        mds.obj <- cmdscale(d,eig=TRUE,
                            k=as.numeric(input$rnaMdsKDim))
                    },
                    warning=function(w) {
                        currentDimRed$mdsPlot <- 
                            ggmessage(paste("An unexpected warning occured.\n",
                                "Try changing your settings."),type="warning")
                    },
                    error=function(e) {
                        currentDimRed$mdsPlot <- 
                            ggmessage(paste("An unexpected error occured.\n",
                                "Try changing your settings."),type="error")
                    },
                    finally="")
                }
            },
            pca = {
                output$rnaMdsPcaSettingsError <- renderUI({div()})
                tryCatch({
                    pca.obj <- prcomp(t(tab),center=input$rnaPcaDoCenter,
                        scale.=input$rnaPcaDoScale,tol=1e-9)
                },
                warning=function(w) {
                    currentDimRed$pcaScreePlot <- currentDimRed$pcaScoresPlot <-
                        currentDimRed$pcaLoadingsPlot <- 
                        currentDimRed$pcaRankedLoadingsPlot <- pcaBiplotPlot <-
                            ggmessage(paste("An unexpected warning occured.\n",
                                "Try changing your settings."),type="warning")
                },
                error=function(e) {
                    currentDimRed$pcaScreePlot <- currentDimRed$pcaScoresPlot <-
                        currentDimRed$pcaLoadingsPlot <- 
                        currentDimRed$pcaRankedLoadingsPlot <- pcaBiplotPlot <-
                            ggmessage(paste("An unexpected error occured.\n",
                                "Try changing your settings."),type="error")
                },
                finally="")
            }
        )
        
        currentDimRed$datMatrix <- currentDimRed$selMatrix <- tab
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
        mds.obj <- currentDimRed$mdsObj
        if (is.null(mds.obj))
            return()
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
            currentDimRed$mdsData <- for.ggplot
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
                labs <- rownames(for.ggplot)
                if (!is.null(currentMetadata$final$alt_id)) {
                    alt_id <- as.character(currentMetadata$final$alt_id)
                    names(alt_id) <- 
                        as.character(currentMetadata$final$sample_id)
                    # TCGA names hack
                    labs <- gsub("\\.","-",alt_id[rownames(for.ggplot)])
                }
                mds <- mds +
                    geom_text_repel(data=for.ggplot,mapping=aes(x=x,y=y,
                        label=labs),size=tsize)
            }
            currentDimRed$mdsPlot <- mds
        }
    })
    
    updateMdsTable <- reactive({
        if (is.null(currentDimRed$mdsData))
            return()
        co <- NULL
        # If from click
        tabClick <- 
            nearPoints(currentDimRed$mdsData,input$rnaMdsPlotClick,
                xvar="x",yvar="y",allRows=TRUE)
        selcl <- which(tabClick$selected_)
        if (length(selcl)>0)
            co <- rownames(tabClick[selcl,])
        
        # If from brush
        tabBrush <- 
            brushedPoints(currentDimRed$mdsData,input$rnaMdsPlotBrush,
                xvar="x",yvar="y",allRows=TRUE)
        selbr <- which(tabBrush$selected_)
        if (length(selbr)>0)
            co <- rownames(tabBrush[selbr,])
        if (isEmpty(co))
            currentDimRed$selMatrix <- currentDimRed$datMatrix
        else
            currentDimRed$selMatrix <- currentDimRed$datMatrix[,co,drop=FALSE]
    })
    
    updatePcaScreePlot <- reactive({
        pca.obj <- currentDimRed$pcaObj
        if (is.null(pca.obj))
            return()
        for.ggplot <- data.frame(
            x=1:length(pca.obj$sdev),
            y=pca.obj$sdev^2
        )
        brs <- 1:nrow(for.ggplot)
        scree <- ggplot() +
            geom_bar(data=for.ggplot,mapping=aes(x=x,y=y),stat="identity",
                fill="deepskyblue",color="deepskyblue3",width=0.5) +
            geom_point(data=for.ggplot,mapping=aes(x=x,y=y),colour="darkblue",
                size=5) +
            geom_line(data=for.ggplot,mapping=aes(x=x,y=y),colour="red2",size=1) +
            xlab(paste("\nPC number")) +
            ylab(paste("Eigenvalue")) +
            theme_bw() +
            theme(
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14),
                axis.text.x=element_text(size=12,face="bold"),
                axis.text.y=element_text(size=12,face="bold")
            ) + 
            scale_x_continuous(breaks=brs)
        currentDimRed$pcaScreePlot <- scree
    })
    
    updatePcaScoresPlot <- reactive({
        pca.obj <- currentDimRed$pcaObj
        if (is.null(pca.obj))
            return()
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
            pvarx <- round(100*pca.obj$sdev[i]^2/sum(pca.obj$sdev^2),1)
            pvary <- round(100*pca.obj$sdev[j]^2/sum(pca.obj$sdev^2),1)
            for.ggplot <- data.frame(
                x=pca.obj$x[,i],
                y=pca.obj$x[,j],
                Condition=classes
            )
            rownames(for.ggplot) <- rownames(pca.obj$x)
            currentDimRed$pcaScoresData <- for.ggplot
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
            pcaScores <- ggplot() +
                geom_point(data=for.ggplot,mapping=aes(x=x,y=y,
                    colour=Condition),size=psize) +
                xlab(paste("\nPrincipal Component ",i," (",pvarx,
                    "% of variance explained)",sep="")) +
                ylab(paste("Principal Component ",j," (",pvary,
                    "% of variance explained)\n",sep="")) +
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
                labs <- rownames(for.ggplot)
                if (!is.null(currentMetadata$final$alt_id)) {
                    alt_id <- as.character(currentMetadata$final$alt_id)
                    names(alt_id) <- 
                        as.character(currentMetadata$final$sample_id)
                    # TCGA names hack
                    labs <- gsub("\\.","-",alt_id[rownames(for.ggplot)])
                }
                pcaScores <- pcaScores +
                    geom_text_repel(data=for.ggplot,mapping=aes(x=x,y=y,
                        label=labs),size=tsize)
            }
            currentDimRed$pcaScoresPlot <- pcaScores
        }
    })
    
    updatePcaLoadingsPlot <- reactive({
        pca.obj <- currentDimRed$pcaObj
        if (is.null(pca.obj))
            return()
        i <- as.numeric(input$rnaDimRedXAxis)
        j <- as.numeric(input$rnaDimRedYAxis)
        if (!is.na(i) && !is.na(j)) {
            labs <- loadedGenomes[[currentMetadata$genome]]$dbGene[
                rownames(pca.obj$rotation)]$gene_name
            pvarx <- round(100*pca.obj$sdev[i]^2/sum(pca.obj$sdev^2),1)
            pvary <- round(100*pca.obj$sdev[j]^2/sum(pca.obj$sdev^2),1)
            for.ggplot <- data.frame(
                x=pca.obj$rotation[,i],
                y=pca.obj$rotation[,j],
                labs=labs
            )
            rownames(for.ggplot) <- rownames(pca.obj$rotation)
            currentDimRed$pcaLoadingsData <- for.ggplot
            if (nrow(for.ggplot)<=100) {
                psize <- 2.5
                tsize <- 3.5
            }
            else {
                psize <- 1.5
                tsize <- 2.5
            }
            pcaLoadings <- ggplot() +
                geom_point(data=for.ggplot,mapping=aes(x=x,y=y),
                    colour="deepskyblue",size=psize) +
                xlab(paste("\nPrincipal Component ",i," (",pvarx,
                    "% of variance explained)",sep="")) +
                ylab(paste("Principal Component ",j," (",pvary,
                    "% of variance explained)\n",sep="")) +
                theme_bw() +
                theme(
                    axis.title.x=element_text(size=14),
                    axis.title.y=element_text(size=14),
                    axis.text.x=element_text(size=12,face="bold"),
                    axis.text.y=element_text(size=12,face="bold"),
                    legend.position="bottom",
                    legend.text=element_text(size=14),
                    legend.key=element_blank()
                )
            if (input$rnaMdsPcaTogglePointNames) {
                if (nrow(for.ggplot)>1000)
                    currentDimRed$tooManyGenes <- TRUE
                else {
                    currentDimRed$tooManyGenes <- FALSE 
                    require(ggrepel)
                    pcaLoadings <- pcaLoadings +
                        geom_text_repel(data=for.ggplot,mapping=aes(x=x,y=y,
                            label=labs),size=tsize)
                }
            }
            else
                currentDimRed$tooManyGenes <- FALSE
            if (input$rnaMdsPcaToggleGeneSearch 
                && !is.null(input$rnaMdsPcaSearchGeneNames)) {
                pcaLoadings <- pcaLoadings +
                    geom_point(data=for.ggplot[input$rnaMdsPcaSearchGeneNames,],
                        mapping=aes(x=x,y=y),colour="red2",size=4.5)
            }
            currentDimRed$pcaLoadingsPlot <- pcaLoadings
        }
    })
    
    updatePcaRankedLoadingsPlot <- reactive({
        pca.obj <- currentDimRed$pcaObj
        if (is.null(pca.obj))
            return()
        i <- as.numeric(input$rnaDimRedXAxis)
        if (!is.na(i)) {
            labs <- loadedGenomes[[currentMetadata$genome]]$dbGene[
                rownames(pca.obj$rotation)]$gene_name
            pvarx <- round(100*pca.obj$sdev[i]^2/sum(pca.obj$sdev^2),1)
            y <- pca.obj$rotation[,i]
            y <- sort(y,index.return=TRUE)
            for.ggplot <- data.frame(
                x=1:length(y$x),
                y=y$x,
                labs=labs[y$ix]
            )
            rownames(for.ggplot) <- rownames(pca.obj$rotation)[y$ix]
            currentDimRed$pcaRankedLoadingsData <- for.ggplot
            if (nrow(for.ggplot)<=100) {
                psize <- 2.5
                tsize <- 3.5
            }
            else {
                psize <- 1.5
                tsize <- 2.5
            }
            pcaRankedLoadings <- ggplot() +
                geom_point(data=for.ggplot,mapping=aes(x=x,y=y),
                    colour="deepskyblue",size=psize) +
                xlab(paste("\nPrincipal Component ",i," (",pvarx,
                    "% of variance explained)",sep="")) +
                ylab("Ranked loadings") +
                theme_bw() +
                theme(
                    axis.title.x=element_text(size=14),
                    axis.title.y=element_text(size=14),
                    axis.text.x=element_text(size=12,face="bold"),
                    axis.text.y=element_text(size=12,face="bold"),
                    legend.position="bottom",
                    legend.text=element_text(size=14),
                    legend.key=element_blank()
                )
            if (input$rnaMdsPcaTogglePointNames) {
                if (nrow(for.ggplot)>1000)
                    currentDimRed$tooManyGenes <- TRUE
                else {
                    currentDimRed$tooManyGenes <- FALSE
                    require(ggrepel)
                    pcaRankedLoadings <- pcaRankedLoadings +
                        geom_text_repel(data=for.ggplot,mapping=aes(x=x,y=y,
                            label=labs),size=tsize)
                }
            }
            else
                currentDimRed$tooManyGenes <- FALSE
            if (input$rnaMdsPcaToggleGeneSearch 
                && !is.null(input$rnaMdsPcaSearchGeneNames)) {
                pcaRankedLoadings <- pcaRankedLoadings +
                    geom_point(data=for.ggplot[input$rnaMdsPcaSearchGeneNames,],
                        mapping=aes(x=x,y=y),colour="red2",size=4.5)
            }
            currentDimRed$pcaRankedLoadingsPlot <- pcaRankedLoadings
        }
    })
    
    updatePcaBiplotPlot <- reactive({
        pca.obj <- currentDimRed$pcaObj
        if (is.null(pca.obj))
            return()
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
            pvarx <- round(100*pca.obj$sdev[i]^2/sum(pca.obj$sdev^2),1)
            pvary <- round(100*pca.obj$sdev[j]^2/sum(pca.obj$sdev^2),1)
            #for.ggplot <- data.frame(
            #    x=pca.obj$x[,i],
            #    y=pca.obj$x[,j],
            #    Condition=classes
            #)
            #rownames(for.ggplot) <- rownames(pca.obj$x)
            #currentDimRed$pcaScoresData <- for.ggplot
            classColors <- currentDimRed$opts$colors
            names(classColors) <- cc
            #if (nrow(for.ggplot)<=100) {
            #    psize <- 3
            #    tsize <- 4
            #}
            #else {
            #    psize <- 2
            #    tsize <- 3
            #}
            
            pcaBiplot <- ggbiplot(
                pcobj=pca.obj,
                choices=c(i,j),
                #groups=classes,
                ellipse=TRUE
            )
                
            pcaBiplot <- pcaBiplot + 
                xlab(paste("\nPrincipal Component ",i," (",pvarx,
                    "% of variance explained)",sep="")) +
                ylab(paste("Principal Component ",j," (",pvary,
                    "% of variance explained)\n",sep="")) +
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
            #if (input$rnaMdsPcaTogglePointNames) {
            #    require(ggrepel)
            #    labs <- rownames(for.ggplot)
            #    if (!is.null(currentMetadata$final$alt_id)) {
            #        alt_id <- as.character(currentMetadata$final$alt_id)
            #        names(alt_id) <- 
            #            as.character(currentMetadata$final$sample_id)
            #        # TCGA names hack
            #        labs <- gsub("\\.","-",alt_id[rownames(for.ggplot)])
            #    }
            #    pcaBiplot <- pcaBiplot +
            #        geom_text_repel(data=for.ggplot,mapping=aes(x=x,y=y,
            #            label=labs),size=tsize)
            #}
            currentDimRed$pcaBiplotPlot <- pcaBiplot
        }
    })
    
    updatePcaScoresTable <- reactive({
        if (is.null(currentDimRed$pcaScoresData))
            return()
        co <- NULL
        # If from click
        tabClick <- 
            nearPoints(currentDimRed$pcaScoresData,input$rnaPcaScoresPlotClick,
                xvar="x",yvar="y",allRows=TRUE)
        selcl <- which(tabClick$selected_)
        if (length(selcl)>0)
            co <- rownames(tabClick[selcl,])
        
        # If from brush
        tabBrush <- 
            brushedPoints(currentDimRed$pcaScoresData,
                input$rnaPcaScoresPlotBrush,xvar="x",yvar="y",allRows=TRUE)
        selbr <- which(tabBrush$selected_)
        if (length(selbr)>0)
            co <- rownames(tabBrush[selbr,])
        if (isEmpty(co))
            currentDimRed$selMatrix <- currentDimRed$datMatrix
        else
            currentDimRed$selMatrix <- currentDimRed$datMatrix[,co,drop=FALSE]
    })
    
    updatePcaLoadingsTable <- reactive({
        if (is.null(currentDimRed$pcaLoadingsData))
            return()
        ro <- NULL
        if (input$rnaMdsPcaToggleGeneSearch 
                && !is.null(input$rnaMdsPcaSearchGeneNames)) {
            ro <- input$rnaMdsPcaSearchGeneNames
        }
        else {
            # If from click
            tabClick <- 
                nearPoints(currentDimRed$pcaLoadingsData,
                    input$rnaPcaLoadingsPlotClick,xvar="x",yvar="y",
                        allRows=TRUE)
            selcl <- which(tabClick$selected_)
            if (length(selcl)>0)
                ro <- rownames(tabClick[selcl,])
            
            # If from brush
            tabBrush <- 
                brushedPoints(currentDimRed$pcaLoadingsData,
                    input$rnaPcaLoadingsPlotBrush,xvar="x",yvar="y",
                        allRows=TRUE)
            selbr <- which(tabBrush$selected_)
            if (length(selbr)>0)
                ro <- rownames(tabBrush[selbr,])
        }
        if (isEmpty(ro))
            currentDimRed$selMatrix <- currentDimRed$datMatrix
        else
            currentDimRed$selMatrix <- 
                currentDimRed$datMatrix[ro,,drop=FALSE]
    })
    
    updatePcaRankedLoadingsTable <- reactive({
        if (is.null(currentDimRed$pcaRankedLoadingsData))
            return()
        ro <- NULL
        if (input$rnaMdsPcaToggleGeneSearch 
                && !is.null(input$rnaMdsPcaSearchGeneNames)) {
            ro <- input$rnaMdsPcaSearchGeneNames
        }
        else {
            # If from click
            tabClick <- 
                nearPoints(currentDimRed$pcaRankedLoadingsData,
                    input$rnaPcaRankedLoadingsPlotClick,xvar="x",yvar="y",
                        allRows=TRUE)
            selcl <- which(tabClick$selected_)
            if (length(selcl)>0)
                ro <- rownames(tabClick[selcl,])
            
            # If from brush
            tabBrush <- 
                brushedPoints(currentDimRed$pcaRankedLoadingsData,
                    input$rnaPcaRankedLoadingsPlotBrush,xvar="x",yvar="y",
                        allRows=TRUE)
            selbr <- which(tabBrush$selected_)
            if (length(selbr)>0)
                ro <- rownames(tabBrush[selbr,])
        }
        if (isEmpty(ro))
            currentDimRed$selMatrix <- currentDimRed$datMatrix
        else
            currentDimRed$selMatrix <- 
                currentDimRed$datMatrix[ro,,drop=FALSE]
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
                            currentDimRed$selMatrix)]$gene_name)
                        res <- data.frame(
                            gene_name=gNames,
                            currentDimRed$selMatrix
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
                    currentDimRed$selMatrix)]$gene_name)
                res <- data.frame(
                    gene_name=gNames,
                    currentDimRed$selMatrix
                )
                write.table(res,file=con,sep="\t",quote=FALSE,row.names=FALSE)
            }
        )
    })
    
    return(list(
        updateMdsPcaPlotColours=updateMdsPcaPlotColours,
        updateMdsPlot=updateMdsPlot,
        updateMdsTable=updateMdsTable,
        updatePcaScreePlot=updatePcaScreePlot,
        updatePcaScoresPlot=updatePcaScoresPlot,
        updatePcaLoadingsPlot=updatePcaLoadingsPlot,
        updatePcaRankedLoadingsPlot=updatePcaRankedLoadingsPlot,
        updatePcaBiplotPlot=updatePcaBiplotPlot,
        updatePcaScoresTable=updatePcaScoresTable,
        updatePcaLoadingsTable=updatePcaLoadingsTable,
        updatePcaRankedLoadingsTable=updatePcaRankedLoadingsTable,
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
        else if (is.null(currentDimRed$selMatrix)) {
            div(
                style="display:inline-block; margin:5px;",
                h4("Please run a variance projection analysis first.")
            )
        }
        else {
            output$rnaMdsPcaGeneSampleTable <- DT::renderDataTable(
                if (is.null(currentDimRed$selMatrix))
                    data.frame(
                        name=character(0),
                        value=numeric(0)
                    )
                else {
                    s <- currentMetadata$source
                    d <- currentMetadata$dataset
                    gNames <- as.character(loadedGenomes[[
                        currentMetadata$genome]]$dbGene[rownames(
                        currentDimRed$selMatrix)]$gene_name)
                    dispMatrix <- currentDimRed$selMatrix
                    if (!is.null(currentMetadata$final$alt_id)) {
                        alt_id <- as.character(currentMetadata$final$alt_id)
                        names(alt_id) <- colnames(currentDimRed$selMatrix)
                        colnames(dispMatrix) <- 
                            gsub("\\.","-",alt_id[colnames(dispMatrix)])
                    }
                    data.frame(
                        gene_name=gNames,
                        dispMatrix
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
    
    output$rnaPcaRankedLoadingsPlot <- renderPlot({
        currentDimRed$pcaRankedLoadingsPlot
    })
    
    output$rnaPcaBiplotPlot <- renderPlot({
        currentDimRed$pcaBiplotPlot
    })
    
    output$helpTooManyRnaMdsPcaGeneNames <- renderUI({
        if (currentDimRed$tooManyGenes)
            list(
                helpText(paste("Point name plotting suppressed due to too ",
                    "many points present on the plot. Please use other ",
                    "filters to check for specific genes"))
            )
        else 
            list(div())
    })
}

mdsPcaTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentDimRed <- allReactiveVars$currentDimRed
    customRnaRegions <- allReactiveVars$customRnaRegions
    
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
    updateMdsTable <- mdsPcaTabPanelReactiveExprs$updateMdsTable
    updatePcaScreePlot <- mdsPcaTabPanelReactiveExprs$updatePcaScreePlot
    updatePcaScoresPlot <- mdsPcaTabPanelReactiveExprs$updatePcaScoresPlot
    updatePcaLoadingsPlot <- mdsPcaTabPanelReactiveExprs$updatePcaLoadingsPlot
    updatePcaRankedLoadingsPlot <- 
        mdsPcaTabPanelReactiveExprs$updatePcaRankedLoadingsPlot
    updatePcaBiplotPlot <- mdsPcaTabPanelReactiveExprs$updatePcaBiplotPlot
    updatePcaScoresTable <- mdsPcaTabPanelReactiveExprs$updatePcaScoresTable
    updatePcaLoadingsTable <- mdsPcaTabPanelReactiveExprs$updatePcaLoadingsTable
    updatePcaRankedLoadingsTable <- 
        mdsPcaTabPanelReactiveExprs$updatePcaRankedLoadingsTable
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
            g <- isolate({input$selectMdsPcaGeneName})
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
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        if (is.null(currentMetadata$final) || is.null(currentDimRed$datMatrix))
            updateSelectizeInput(session,"rnaMdsPcaSearchGeneNames",
                choices=NULL,
                server=TRUE
            )
        else {
            if (!is.null(currentDimRed$datMatrix)) {
                G <- rownames(currentDimRed$datMatrix)
                nG <- loadedGenomes[[currentMetadata$genome]]$dbGene[
                    G]$gene_name
                nG <- as.character(nG)
                names(G) <- nG
                g <- isolate({input$rnaMdsPcaSearchGeneNames})
                i <- grep(paste0("^",g),nG,perl=TRUE)
                if (length(i)>0) {
                    updateSelectizeInput(session,"rnaMdsPcaSearchGeneNames",
                        choices=G[i],
                        selected=g,
                        server=TRUE
                    )
                }
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
        if (isEmpty(currentDimRed$mdsObj) && isEmpty(currentDimRed$pcaObj)) {
            updateSelectizeInput(session,"rnaDimRedXAxis",
                choices=NULL)
            updateSelectizeInput(session,"rnaDimRedYAxis",
                choices=NULL)
        }
        else if (!isEmpty(currentDimRed$mdsObj) 
            && isEmpty(currentDimRed$pcaObj)) {
            choices <- as.character(1:as.numeric(input$rnaMdsKDim))
            names(choices) <- paste("PC",1:as.numeric(input$rnaMdsKDim))
            updateSelectizeInput(session,"rnaDimRedXAxis",
                choices=choices,selected="1")
            updateSelectizeInput(session,"rnaDimRedYAxis",
                choices=choices,selected="2")
        }
        else if (isEmpty(currentDimRed$mdsObj) 
            && !isEmpty(currentDimRed$pcaObj)) {
            npca <- length(currentDimRed$pcaObj$sdev)
            choices <- as.character(1:npca)
            names(choices) <- paste("PC",1:npca)
            updateSelectizeInput(session,"rnaDimRedXAxis",
                choices=choices,selected="1")
            updateSelectizeInput(session,"rnaDimRedYAxis",
                choices=choices,selected="2")
        }
    })
    
    observe({
        updateMdsPcaPlotColours()
        updateMdsPlot()
        updateMdsTable()
        updatePcaScreePlot()
        updatePcaScoresPlot()
        updatePcaLoadingsPlot()
        updatePcaRankedLoadingsPlot()
        updatePcaBiplotPlot()
        updatePcaScoresTable()
        updatePcaLoadingsTable()
        updatePcaRankedLoadingsTable()
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

