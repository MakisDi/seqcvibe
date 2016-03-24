# server.R

# Load init script
source("config/init.R")

shinyServer(
    function(input,output,session) {
        # Make %#^%$^%$@( globals visible AND changeable
        makeReactiveBinding("loadedGenomes")
        makeReactiveBinding("loadedData")
        
        ########################################################################

        # First main panel - Data selector
        
        currentMetadata <- reactiveValues(
            source=sources[1],
            dataset=datasets[1],
            class=classes,
            genome=genomes[1],
            metadata=NULL,
            final=NULL
        )
        
        dataSelectorMessages <- reactiveValues(
            messages=list(
                list(
                    type="INFO",
                    msg=paste(getTime("INFO"),"Welcome to the data selector ",
                        "of BigSeqCVis! This is an info message. Make your ",
                        "selections on the left.")
                )
            )
        )
        
        upateCurrentSource <- eventReactive(input$dataSource,{
            currentMetadata$source <- input$dataSource
        })
        
        updateCurrentMetadata <- eventReactive(input$dataDataset,{
            s <- currentMetadata$source
            d <- input$dataDataset
            currentMetadata$dataset <- d
            currInd <- which(as.character(metadata$source)==s 
                & as.character(metadata$dataset)==d)
            currentMetadata$class <- 
                unique(as.character(metadata$class[currInd]))
            currentMetadata$metadata <- metadata[currInd,]
            currentMetadata$genome <- as.character(metadata$genome[currInd])[1]
            if (!is.na(currentMetadata$genome)
                && is.null(loadedGenomes[[currentMetadata$genome]]$dbGene)) {
                load(file.path("genome",currentMetadata$genome,"gene.rda"))
                loadedGenomes[[currentMetadata$genome]]$dbGene <<- gene
                #gn <- names(loadedGenomes[[currentMetadata$genome]]$dbGene)
                gn <- names(gene)
                #names(gn) <- as.character(elementMetadata(
                #    loadedGenomes[[currentMetadata$genome]]$dbGene)$gene_name)
                names(gn) <- as.character(elementMetadata(gene)$gene_name)
                loadedGenomes[[currentMetadata$genome]]$geneNames <<- gn
                dataSelectorMessages <- updateMessages(
                    dataSelectorMessages,
                    type="SUCCESS",
                    msg=paste(getTime("SUCCESS"),"Genome ",
                        currentMetadata$genome," genes loaded!", sep="")
                )
            }
            if (!isEmpty(s) && !isEmpty(d))
                dataSelectorMessages <- updateMessages(
                    dataSelectorMessages,
                    type="INFO",
                    msg=paste(getTime("INFO"),"Current source: ",s,
                        ", Current dataset: ",d,", Classes: ", 
                        paste(currentMetadata$class,collapse=", "),
                        ", Total samples: ",nrow(currentMetadata$metadata),
                        sep="")
                )
        })
        
        createDataset <- eventReactive(input$createDataset,{
            s <- currentMetadata$source
            d <- currentMetadata$dataset
            c <- currentMetadata$class
            m <- currentMetadata$metadata
            currentMetadata$final <- do.call("rbind",lapply(c,function(x) {
                if (input[[paste("sampleSelectType_",x,sep="")]]=="all")
                    return(m[m$class==x,])
                else if (input[[paste("sampleSelectType_",x,sep="")]]=="none") {
                    #m <- match(x,currentMetadata$class)
                    #currentMetadata$class <- currentMetadata$class[-m]
                    return(NULL)
                }
                else {
                    sel <- 
                        input[[paste("classTable_",x,"_rows_selected",sep="")]]
                    if (length(sel)==0)
                        return(m[m$class==x,])
                    else {
                        mm <- m[m$class==x,]
                        return(mm[sel,])
                    }
                }
            }))
            
            msgpart <- paste(lapply(c,function(x,m) {
                paste("Number of samples in ",x,": ",nrow(m[m$class==x,]))
            },currentMetadata$final),collapse=", ")
            dataSelectorMessages <- updateMessages(
                dataSelectorMessages,
                type="SUCCESS",
                msg=paste(getTime("SUCCESS"),"Dataset created! Source: ",
                    currentMetadata$source,", Dataset: ",
                    currentMetadata$dataset,", Classes: ",
                    paste(currentMetadata$class,collapse=", "),", ",msgpart,
                    ", Total samples: ",nrow(currentMetadata$final),sep="")
            )
        })
        
        clearDataset <- eventReactive(input$clearDataset,{
            currentMetadata$final <- NULL
            dataSelectorMessages <- updateMessages(
                dataSelectorMessages,
                type="WARNING",
                msg=paste(getTime("WARNING"),"Dataset cleared!")
            )
        })
        
        loadDataset <- eventReactive(input$loadDataset,{
            s <- input$dataSource
            d <- input$dataDataset
            c <- currentMetadata$class
            if (!isEmpty(s) && !isEmpty(d)) {
                if (is.null(loadedData[[s]][[d]])) {
                    dataSelectorMessages <- updateMessages(
                        dataSelectorMessages,
                        type="INFO",
                        msg=paste(getTime("INFO"),
                            "Loading dataset ",d," from ",s,sep="")
                    )                    
                    load(dataFiles[[s]][[d]])
                    # ^$#$#@@#$^%$$!!!!
                    loadedData[[s]][[d]] <<- b2c.out
                    dataSelectorMessages <- updateMessages(
                        dataSelectorMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),
                            "Loaded ",d," from ",s,"!",sep="")
                    )
                }
                else
                    dataSelectorMessages <- updateMessages(
                        dataSelectorMessages,
                        type="INFO",
                        msg=paste(getTime("INFO"),
                            "Dataset ",d," from ",s," already loaded.",sep="")
                    )
                }
        })
        
        classDataTables <- reactive({
            s <- input$dataSource
            d <- input$dataDataset
            c <- currentMetadata$class
            lapply(c,function(x,s,d) {
                tab <- metadata[which(as.character(metadata$source)==s 
                    & as.character(metadata$dataset)==d
                    & as.character(metadata$class)==x),
                    c("sample_id","alt_id","norm_factor","library_strategy")]
                output[[paste("classTable",x,sep="_")]] <- 
                    DT::renderDataTable(
                        tab,
                        class="display compact",
                        rownames=FALSE,
                        options=list(
                            searchHighlight=TRUE,
                            pageLength=10,
                            lengthMenu=c(10,20,50,100)
                        )
                    )
            },s,d)
        })
        
        handleSampleSelection <- reactive({
            s <- currentMetadata$source
            d <- currentMetadata$dataset
            c <- currentMetadata$class
            lapply(c,function(x) {
                observeEvent(input[[paste("clearSelection_",x,sep="")]],{
                    proxy <- dataTableProxy(paste("classTable_",x,sep=""))
                    selectRows(proxy,NULL)
                })
            })
            lapply(c,function(x) {
                N <- 1:nrow(metadata[which(as.character(metadata$source)==s 
                    & as.character(metadata$dataset)==d
                    & as.character(metadata$class)==x),])
                observeEvent(input[[paste("invertSelection_",x,sep="")]],{
                    sel <- input[[paste("classTable_",x,
                        "_rows_selected",sep="")]]
                    if (length(sel)>0) {
                        N <- N[-sel]
                        proxy <- dataTableProxy(paste("classTable_",x,sep=""))
                        selectRows(proxy,N)
                    }
                })
            })
        })
        
        updateDataSelectorMessages <- reactive({
            output$dataSelectorMessages <- renderUI({
                lapply(dataSelectorMessages$messages,function(x) {
                    switch(x$type,
                        INFO = {
                            div(class="info-box",x$msg)
                        },
                        SUCCESS = {
                            div(class="success-box",x$msg)
                        },
                        WARNING = {
                            div(class="warn-box",x$msg)
                        },
                        ERROR = {
                            div(class="error-box",x$msg)
                        }
                    )
                })
            }) 
        })
        
        output$currentDatasetTable <- DT::renderDataTable(
            if (is.null(currentMetadata$final))
                data.frame(
                    sample_id=character(0),
                    alt_id=character(0),
                    class=character(0),
                    norm_factor=numeric(0),
                    quality=integer(0)
                )
            else    
                as.data.frame(currentMetadata$final[,
                    c("sample_id","alt_id","class","norm_factor")]),
            class="display compact",
            rownames=FALSE,
            options=list(
                searchHighlight=TRUE,
                pageLength=10,
                lengthMenu=c(10,20,50,100)
            )
        )
        
        # Init the data source
        output$dataSource <- renderUI({
            selectInput(
                inputId="dataSource",
                label="Select data source",
                choices=unique(as.character(metadata$source))
            )
        })
        
        # Init the dataset
        output$dataDataset <- renderUI({
            selectInput(
                inputId="dataDataset",
                label="Select dataset",
                choices=unique(as.character(metadata$dataset[
                    which(as.character(metadata$source)==input$dataSource)]))
            )
        })
        
        # Fill genome
        output$dataGenome <- renderUI({
            disabled(textInput(
                inputId="dataGenome",
                label="Genome",
                value=currentMetadata$genome
            ))
        })
        
        output$dataSelectHint <- renderUI({
            s <- input$dataSource
            d <- input$dataDataset
            if (!isEmpty(s) && !isEmpty(d)) {
                if (!is.null(loadedData[[s]][[d]]))
                    list(
                        h4("Select samples for selected dataset"),
                        helpText(paste("Select 'Custom samples' for ",
                            "customized sample selections or 'All samples' to ",
                            "select all samples from each class."))
                    )
            }
        }) 
        
        output$dataCustomSamples <- renderUI({
            s <- input$dataSource
            d <- input$dataDataset
            c <- unique(as.character(metadata$class[
                    which(as.character(metadata$source)==s 
                        & as.character(metadata$dataset)==d)]))
            if (!isEmpty(s) && !isEmpty(d)) {
            
            if (!is.null(loadedData[[s]][[d]])) {
            classDataTables()
            lapply(c,function(x,s,d) {
                list(
                    radioButtons(
                        inputId=paste("sampleSelectType",x,sep="_"),
                        label=paste("Select",x,"samples"),
                        inline=TRUE,
                        choices=list(
                            "All samples"="all",
                            "Custom samples"="custom",
                            "No samples"="none"
                        )
                    ),
                    conditionalPanel(
                        condition=paste("input.sampleSelectType_",x,
                            "=='custom'",sep=""),
                        div(
                            class="small table-container",
                            DT::dataTableOutput(paste("classTable",x,sep="_")),
                            div(
                                class="pull-left",
                                style="display:block; margin:5px;",
                                actionButton(
                                    paste("clearSelection_",x,sep=""),
                                    paste("Clear selection"),
                                    class="btn-xs"
                                )
                            ),
                            div(
                                class="pull-left",
                                style="display:block; margin:5px;",
                                actionButton(
                                    paste("invertSelection_",x,sep=""),
                                    paste("Invert selection"),
                                    class="btn-xs"
                                )
                            )
                        )
                    ),
                    hr()
                )
            },s,d) }}
        })
        
        observe({
            loadDataset()
        })
        
        observe({
            updateDataSelectorMessages()
            upateCurrentSource()
            updateCurrentMetadata()
            handleSampleSelection()
            createDataset()
            clearDataset()
        })
        
        observe({
            s <- currentMetadata$source
            d <- currentMetadata$dataset
            if (is.null(loadedData[[s]][[d]]))
                shinyjs::disable("createDataset")
            else
                shinyjs::enable("createDataset")
        })
        
        ########################################################################
        
        # Second main panel - Gene explorer
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
        
        geneExplorerMessages <- reactiveValues(
            messages=list(
                list(
                    type="INFO",
                    msg=paste(getTime("INFO"),"Welcome to the gene selector ",
                        "of BigSeqCVis! This is an info message. Make your ",
                        "selections on the left.")
                )
            )
        )
        
        updateCurrentGenes <- eventReactive(input$geneGeneName,{
            g <- input$geneGeneName
            dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
            if (length(g)>=length(currentGenes$genes)) {
                currentGenes$genes <- g
                if (!isEmpty(g)) {
                    mm <- g[length(g)]
                    currentGenes$coords$chromosome <- 
                        c(currentGenes$coords$chromosome,
                            as.character(seqnames(dbGene[mm]))[1])
                    currentGenes$coords$start <- 
                        c(currentGenes$coords$start,start(dbGene[mm]))
                    currentGenes$coords$end <- 
                        c(currentGenes$coords$end,end(dbGene[mm]))
                    currentGenes$coords$strand <- 
                        c(currentGenes$coords$strand,
                            as.character(strand(dbGene[mm]))[1])
                    currentGenes$coords$name <- 
                        c(currentGenes$coords$name,
                            as.character(dbGene[mm]$gene_name))
                    geneExplorerMessages <- updateMessages(
                        geneExplorerMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),"Gene ",g[length(g)],
                            " (",as.character(dbGene[mm]$gene_name),") ",
                            "added to the known gene list! Number of genes ",
                            "to plot is now ",length(currentGenes$genes),sep="")
                    )
                }
            }
            else if (length(g)<length(currentGenes$genes) || isEmpty(g)) {
                if (isEmpty(g)) {
                    currentGenes$genes <- NULL
                    currentGenes$coords$chromosome <- NULL
                    currentGenes$coords$start <- NULL
                    currentGenes$coords$end <- NULL
                    currentGenes$coords$strand <- NULL
                    currentGenes$coords$name <- NULL
                    geneExplorerMessages <- updateMessages(
                        geneExplorerMessages,
                        type="WARNING",
                        msg=paste(getTime("WARNING"),"All genes removed from ",
                            "the known gene list!","sep=")
                    )
                }
                else {
                    prev <- currentGenes$genes
                    removed <- setdiff(prev,g)
                    ii <- match(removed,prev)
                    gii <- currentGenes$coords$name[ii]
                    currentGenes$genes <- g
                    currentGenes$coords$chromosome <- 
                            currentGenes$coords$chromosome[-ii]
                    currentGenes$coords$start <- currentGenes$coords$start[-ii]
                    currentGenes$coords$end <- currentGenes$coords$end[-ii]
                    currentGenes$coords$strand <- 
                        currentGenes$coords$strand[-ii]
                    currentGenes$coords$name <- currentGenes$coords$name[-ii]
                    geneExplorerMessages <- updateMessages(
                        geneExplorerMessages,
                        type="WARNING",
                        msg=paste(getTime("WARNING"),"Gene ",removed," (",gii,
                            ") removed from the known gene list! Number of ",
                            "gene to plot is now ",length(currentGenes$genes),
                            sep="")
                    )
                }
            }
        })
        
        updateFlanks <- eventReactive({
            input$upstreamFlank
            input$downstreamFlank
        },{
            flankValidate <- c(
                upstream=FALSE,
                downstream=FALSE
            )
            errMsg <- c(
                upstream="",
                downstream=""
            )
            u <- input$upstreamFlank
            d <- input$downstreamFlank
            if (u!=currentOpts$flank[1]) {
                u <- as.integer(u)
                if (is.na(suppressWarnings(u)) || isEmpty(u) || u<=0) {
                    flankValidate["upstream"] <- TRUE
                    errMsg["upstream"] <- paste("Upstream flanking region",
                        "must be a positive integer!")
                }
            }
            if (d!=currentOpts$flank[2]) {
                d <- as.integer(d)
                if (is.na(suppressWarnings(d)) || isEmpty(d) || d<=0) {
                    flankValidate["downstream"] <- TRUE
                    errMsg["downstream"] <- paste("Downstream flanking region",
                        "must be a positive integer!")
                }
            }
            if (any(flankValidate)) {
                output$regionFlankError <- renderUI({
                    div(class="error-message",paste(errMsg,collapse="\n"))
                })
            }
            else {
                output$regionFlankError <- renderUI({div()})
                if (u!=currentOpts$flank[1]) {
                    currentOpts$flank[1] <- u
                    geneExplorerMessages <- updateMessages(
                        geneExplorerMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),"Upstream flanking region",
                            " changed! New upstream flanking: ",u,sep="")
                    )
                }
                if (d!=currentOpts$flank[2]) {
                    currentOpts$flank[2] <- d
                    geneExplorerMessages <- updateMessages(
                        geneExplorerMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),"Downstream flanking ",
                            "region changed! New downstream flanking: ",d,
                            sep="")
                    )
                }
            }
        })
        
        addTheCustomRegion <- eventReactive(input$addCustomRegion,{
            regionValidate <- c(
                name=FALSE,
                strand=FALSE,
                chrom=FALSE,
                start=FALSE,
                end=FALSE,
                startBeforeEnd=FALSE,
                regionExistsName=FALSE,
                regionExistsLocus=FALSE
            )
            errMsg <- c(
                name="",
                strand="",
                chrom="",
                start="",
                end="",
                startBeforeEnd="",
                regionExistsName="",
                regionExistsLocus=""
            )
            if (input$customName=="") {
                regionValidate["name"] <- TRUE
                errMsg["name"] <- "A region name is required!"
            
            }
            if (input$customStrand=="") {
                regionValidate["strand"] <- TRUE
                errMsg["strand"] <- "A region strand is required!"
            }
            if (input$customChr=="") {
                regionValidate["chrom"] <- TRUE
                errMsg["chrom"] <- "A region chromosome is required!"
            
            }
            if (input$customStart=="" || as.integer(input$customStart)<=0
                || is.na(suppressWarnings(as.numeric(input$customStart)))) {
                regionValidate["start"] <- TRUE
                errMsg["start"] <- 
                    "Region start base must be a valid positive integer!"
            }
            if (input$customEnd=="" || as.numeric(input$customEnd)<=0
                || is.na(suppressWarnings(as.integer(input$customEnd)))) {
                regionValidate["end"] <- TRUE
                errMsg["end"] <- 
                    "Region end base must be a valid positive integer!"
            }
            if (!(regionValidate["start"] || regionValidate["end"])) {
                # Check start before end
                s <- as.integer(input$customStart)
                e <- as.integer(input$customEnd)
                if (e-s<=0) {
                    regionValidate["startBeforeEnd"] <- TRUE
                    errMsg["startBeforeEnd"] <- paste("Region start must be ",
                        "smaller than region end and the region length must ",
                        "be greater than 1!",sep="")
                }
                
                # Check the new region does not already exists
                cc <- customRegions$chromosome
                ss <- customRegions$start
                ee <- customRegions$end
                tt <- customRegions$strand
                if (input$customName %in% customRegions$name) {
                    regionValidate["regionExistsName"] <- TRUE
                    errMsg["regionExistsName"] <- paste("Custom regions must ",
                        "have each a unique name! ",input$customName," ",
                        "already exists.",sep="")
                }
                if (length(ss)>0) {
                    for (i in 1:length(ss)) {
                        if (input$customChr==cc[i] && input$customStart==ss[i] 
                            && input$customEnd==ee[i] 
                            && input$customStrand==tt[i])
                        regionValidate["regionExistsLocus"] <- TRUE
                        errMsg["regionExistsLocus"] <- 
                            "Custom regions must have unique coordinates!"
                    }
                }
            }
            
            if (any(regionValidate)) {
                output$customRegionError <- renderUI({
                    div(class="error-message",paste(errMsg,collapse="\n"))
                })
                return("")
            }
            else {
                output$customRegionError <- renderUI({div()})
                customRegions$chromosome <- c(customRegions$chromosome,
                    input$customChr)
                customRegions$start <- c(customRegions$start,input$customStart)
                customRegions$end <- c(customRegions$end,input$customEnd)
                customRegions$strand <- c(customRegions$strand,
                    input$customStrand)
                customRegions$name <- c(customRegions$name,input$customName)
                
                geneExplorerMessages <- updateMessages(
                    geneExplorerMessages,
                    type="SUCCESS",
                    msg=paste(getTime("SUCCESS"),"Custom region ",
                        input$customName,
                        " added to the custom regions to plot! Number of ",
                        "custom regions to plot is now ",
                        length(customRegions$name),sep="")
                )
            }
        })
        
        removeTheCustomRegion <- eventReactive(input$removeCustomRegion,{
            output$customRegionError <- renderUI({div()})
            j <- as.integer(input$customRegionList_rows_selected)
            if (length(j)>0) {
                ii <- customRegions$name[j]
                customRegions$chromosome <- customRegions$chromosome[-j]
                customRegions$start <- customRegions$start[-j]
                customRegions$end <- customRegions$end[-j]
                customRegions$strand <- customRegions$strand[-j]
                customRegions$name <- customRegions$name[-j]
                geneExplorerMessages <- updateMessages(
                    geneExplorerMessages,
                    type="WARNING",
                    msg=paste(getTime("WARNING"),"Custom region ",ii," ",
                        " removed from the custom regions to plot! Number of ",
                        "custom regions to plot is now ",
                        length(customRegions$name),sep="")
                )
            }
        })
        
        createGeneProfile <- eventReactive(input$createGeneProfile,{
            if (is.null(currentMetadata$final)) {
                output$geneExplorerError <- renderUI({
                    div(class="error-message",paste("You must create a ",
                        "dataset first!",sep=""))
                })
            }
            else {
                output$geneExplorerError <- renderUI({div()})
                dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
                isolate(input$geneType)
                knownGenes <- as.list(currentGenes$genes)
                names(knownGenes) <- currentGenes$genes
                customGenes <- as.list(makeGRangesFromDataFrame(
                    df=data.frame(
                        chromosome=as.character(customRegions$chromosome),
                        start=as.numeric(customRegions$start),
                        end=as.numeric(customRegions$end),
                        strand=as.character(customRegions$strand),
                        gene_id=as.character(customRegions$name)
                    ),
                    keep.extra.columns=TRUE
                ))
                names(customGenes) <- as.character(customRegions$name)
                theGenes <- c(knownGenes,customGenes)
                progress <- shiny::Progress$new()
                progress$initialize(
                    session,
                    min=0,
                    max=length(unique(as.character(
                        currentMetadata$finalclass)))*length(theGenes)*2+2
                )
                progress$set(message="",value=0)
                on.exit(progress$close())
                
                updateProgress <- function(value=NULL,detail=NULL) {
                    if (is.null(value)) {
                        value <- progress$getValue()
                        value <- value + 1
                    }
                    progress$set(value=value,detail=detail)
                }
                
                ggProf <- tryCatch(
                    getProfile(
                        gene=theGenes,
                        flank=as.integer(currentOpts$flank),
                        source=as.character(currentMetadata$source),
                        dataset=as.character(currentMetadata$dataset),
                        #class=as.character(currentMetadata$class),
                        class=unique(as.character(currentMetadata$final$class)),
                        sumStat=as.character(currentOpts$sumStat),
                        config=currentMetadata$final,
                        dbGene=dbGene,
                        trim=as.numeric(currentOpts$trim),
                        messageContainer=geneExplorerMessages,
                        progressFun=updateProgress,
                        rc=RC
                    ),
                    #warning=function(w) {
                    #    geneExplorerMessages <- updateMessages(
                    #        geneExplorerMessages,
                    #        type="WARNING",
                    #        msg=paste(getTime("WARNING"),"A warning occured ",
                    #            "while calculating profiles from bigWig ",
                    #            "files. Please contact the administrator ",
                    #            "stating the following warning: ",w,sep="")
                    #    )
                    #},
                    error=function(e) {
                        geneExplorerMessages <- updateMessages(
                            geneExplorerMessages,
                            type="WARNING",
                            msg=paste(getTime("WARNING"),"An error occured ",
                                "while calculating profiles from bigWig ",
                                "files. Will try with BAM if present. Please ",
                                "contact the administrator stating the ",
                                "following error: ",e,sep="")
                        )
                        tryCatch(
                            getProfile(
                                gene=theGenes,
                                flank=as.integer(currentOpts$flank),
                                source=as.character(currentMetadata$source),
                                dataset=as.character(currentMetadata$dataset),
                                #class=as.character(currentMetadata$class),
                                class=unique(as.character(
                                    currentMetadata$final$class)),
                                sumStat=as.character(currentOpts$sumStat),
                                config=currentMetadata$final,
                                dbGene=dbGene,
                                trim=as.numeric(currentOpts$trim),
                                fromBam=TRUE,
                                messageContainer=geneExplorerMessages,
                                progressFun=updateProgress,
                                rc=RC
                            ),
                            #warning=function(w) {
                            #    geneExplorerMessages <- updateMessages(
                            #        geneExplorerMessages,
                            #        type="WARNING",
                            #        msg=paste(getTime("WARNING"),"A warning ",
                            #            "occured while calculating profiles ",
                            #            "from BAM files. Please contact the ",
                            #            "administrator stating the following ",
                            #            "warning: ",w,sep="")
                            #    )
                            #},
                            error=function(e) {
                                geneExplorerMessages <- updateMessages(
                                    geneExplorerMessages,
                                    type="ERROR",
                                    msg=paste(getTime("ERROR"),"An error ",
                                        "occured while calculating profiles ",
                                        "from BAM files. Please contact the ",
                                        "administrator stating the following ",
                                        "error: ",e,sep="")
                                )
                            },finally=""
                        )
                    },finally=""
                )
                
                ggProf <- ggProf + 
                    facet_wrap(~ Locus,scales="free") +
                    scale_color_manual(values=currentOpts$colours) + 
                    scale_fill_manual(values=currentOpts$colours)
                genePlots$geneProfile <- ggProf
                genePlots$rendered <- TRUE
            }
        })
        
        updateClassColours <- reactive({
            c <- currentMetadata$class
            currentOpts$colours <- baseColours[c]
        })
        
        updateGeneColours <- reactive({
            c <- currentMetadata$class
            lapply(c,function(x) {
                observeEvent(input[[paste("geneColour_",x,sep="")]],{
                    newc <- input[[paste("geneColour_",x,sep="")]]
                    if (!is.null(currentOpts$colours[x]) 
                        && newc!=currentOpts$colours[x]) {
                        currentOpts$colours[x] <- newc
                        geneExplorerMessages <- updateMessages(
                            geneExplorerMessages,
                            type="SUCCESS",
                            msg=paste(getTime("SUCCESS"),"Class ",x," colour ",
                            "changed! New colour: ",newc,sep="")
                        )
                    }
                })
            })
        })
        
        updateGeneExplorerMessages <- reactive({
            output$geneExplorerMessages <- renderUI({
                lapply(geneExplorerMessages$messages,function(x) {
                    switch(x$type,
                        INFO = {
                            div(class="info-box",x$msg)
                        },
                        SUCCESS = {
                            div(class="success-box",x$msg)
                        },
                        WARNING = {
                            div(class="warn-box",x$msg)
                        },
                        ERROR = {
                            div(class="error-box",x$msg)
                        }
                    )
                })
            }) 
        })
        
        updateGeneNames <- reactive({
            g <- isolate({input$geneGeneName})
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"geneGeneName",
                    choices=geneNames[i],
                    selected=g,
                    server=TRUE)
            }
        })
        
        # Very dirty hack to update reactive values holding gene coordinates as
        # there is a problem when selectize input get empty, see also this:
        # http://stackoverflow.com/questions/26803536/shiny-how-to-update-a-reactivevalues-object
        dirtyHack <- reactive({
            currentGenes$genes <- NULL
            currentGenes$coords$chromosome <- NULL
            currentGenes$coords$start <- NULL
            currentGenes$coords$end <- NULL
            currentGenes$coords$strand <- NULL
            currentGenes$coords$name <- NULL
            geneExplorerMessages <- updateMessages(
                geneExplorerMessages,
                type="WARNING",
                msg=paste(getTime("WARNING"),"All genes removed from ",
                    "the known gene list!",sep="")
            )
            return(data.frame(
                chromosome=character(),
                start=integer(),
                end=integer(),
                strand=character(),
                name=character()
            ))
        })
         
        output$knownGeneList <- DT::renderDataTable(
            if (isEmpty(input$geneGeneName))
                isolate(dirtyHack())
            else
                data.frame(
                    chromosome=currentGenes$coords$chromosome,
                    start=currentGenes$coords$start,
                    end=currentGenes$coords$end,
                    strand=currentGenes$coords$strand,
                    name=currentGenes$coords$name
                ),
            class="display compact",
            rownames=FALSE,
            options=list(
                dom='t',
                paging=FALSE
            )
        )
        
        output$customRegionList <- DT::renderDataTable(
            data.frame(
                chromosome=customRegions$chromosome,
                start=customRegions$start,
                end=customRegions$end,
                strand=customRegions$strand,
                name=customRegions$name
            ),
            class="display compact",
            rownames=FALSE,
            options=list(
                dom='t',
                paging=FALSE#,
                #columnDefs=list(list(
                #    targets=0,
                #    visible=FALSE,
                #    searchable=FALSE
                #))
            )
        )
            
        output$setChrs <- renderUI({
            selectizeInput(
                inputId="customChr",
                label="", 
                choices=c("",
                    getValidChromosomes(currentMetadata$genome)
                ),
                options=list(
                    placeholder="Chrom"
                )
            )
        })
        
        output$geneExplorerColours <- renderUI({
            if (!is.null(currentMetadata$final)) {
                c <- unique(as.character(currentMetadata$final$class))
                lapply(1:length(c),function(i,c) {
                    colourInput(
                        inputId=paste("geneColour_",c[i],sep=""),
                        label=paste("Select colour for",c[i]),
                        value=baseColours[i]
                    )
                },c)
            }
        })
        
        output$geneExplorerMessages <- renderUI({
            lapply(geneExplorerMessages$messages,function(x) {
                switch(x$type,
                    INFO = {
                        div(class="info-box",x$msg)
                    },
                    SUCCESS = {
                        div(class="success-box",x$msg)
                    },
                    WARNING = {
                        div(class="warn-box",x$msg)
                    },
                    ERROR = {
                        div(class="error-box",x$msg)
                    }
                )
            })
        })
        
        output$geneProfile <- renderPlot({
            genePlots$geneProfile
        })
        
        output$exportGenePDF <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("gene_plot_",tt,".pdf", sep='')
            },
            content=function(con) {
                ggsave(filename=con,plot=genePlots$geneProfile,
                    width=10,height=7)
            }
        )
        
        output$exportGenePNG <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("gene_plot_",tt,".png", sep='')
            },
            content=function(con) {
                ggsave(filename=con,plot=genePlots$geneProfile,
                    width=10,height=7)
            }
        )
        
        observe({
            if (isEmpty(input$geneGeneName) && length(customRegions$name)==0)
                shinyjs::disable("createGeneProfile")
            else
                shinyjs::enable("createGeneProfile")
        })
        
        observe({
            if (input$geneSumStatType!="trimmed") {
                shinyjs::disable("geneTrimPct")
                if (isEmpty(input$geneGeneName) 
                    && length(customRegions$name)==0)
                    shinyjs::disable("createGeneProfile")
                else
                    shinyjs::enable("createGeneProfile")
            }
            else {
                shinyjs::enable("geneTrimPct")
                trimp <- as.numeric(input$geneTrimPct)
                if (trimp<0 || trimp>0.5) {
                    output$geneExplorerError <- renderUI({
                        div(class="error-message",paste("The trimming ",
                            "must be between 0 and 0.5!",sep=""))
                    })
                    shinyjs::disable("createGeneProfile")
                }
                else {
                    output$geneExplorerError <- renderUI({div()})
                    if (isEmpty(input$geneGeneName) 
                        && length(customRegions$name)==0)
                        shinyjs::disable("createGeneProfile")
                    else
                        shinyjs::enable("createGeneProfile")
                    currentOpts$trim <- trimp
                }
            }
        })
        
        observe({
            if (length(customRegions$name)==0)
                shinyjs::disable("removeCustomRegion")
            else
                shinyjs::enable("removeCustomRegion")
        })
        
        observe({
            updateGeneNames()
            updateCurrentGenes()
        })
        
        observe({
            customRegions <- addTheCustomRegion()
            customRegions <- removeTheCustomRegion()
        })
        
        observe({
            updateFlanks()
            updateClassColours()
            updateGeneColours()
            updateGeneExplorerMessages()
        })
        
        observe({
            tryCatch({
                shinyjs::disable("createGeneProfile")
                createGeneProfile()
            },error=function(e) {
                genePlots$rendered <- FALSE
                return()
            },
            finally={
                if (isEmpty(input$geneGeneName) 
                    && length(customRegions$name)==0)
                    shinyjs::disable("createGeneProfile")
                else
                    shinyjs::enable("createGeneProfile")
            })
        })
       
        ########################################################################
        
        # Fourth main panel - Area explorer
        
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
        
        updateAreaExplorerMessages <- reactive({
            output$areaExplorerMessages <- renderUI({
                lapply(areaExplorerMessages$messages,function(x) {
                    switch(x$type,
                        INFO = {
                            div(class="info-box",x$msg)
                        },
                        SUCCESS = {
                            div(class="success-box",x$msg)
                        },
                        WARNING = {
                            div(class="warn-box",x$msg)
                        },
                        ERROR = {
                            div(class="error-box",x$msg)
                        }
                    )
                })
            }) 
        })
        
        updateAreaColours <- reactive({
            c <- currentMetadata$class
            lapply(c,function(x) {
                observeEvent(input[[paste("areaColour_",x,sep="")]],{
                    newc <- input[[paste("areaColour_",x,sep="")]]
                    if (newc!=currentOpts$colours[x]) {
                        currentOpts$colours[x] <- newc
                        areaExplorerMessages <- updateMessages(
                            areaExplorerMessages,
                            type="SUCCESS",
                            msg=paste(getTime("SUCCESS"),"Class ",x," colour ",
                            "changed! New colour: ",newc,sep="")
                        )
                    }
                })
            })
        })
        
        updateCurrentArea <- eventReactive({
            input$customChrA
            input$customStartA
            input$customEndA
        },{
            customArea$chromosome <- input$customChrA
            customArea$start <- input$customStartA
            customArea$end <- input$customEndA
            customArea$name <- paste(customArea$chromosome,":",
                customArea$start,"-",customArea$end,sep="")
            customArea$strand <- "*"
            if (!any(
                isEmpty(input$customChrA),
                isEmpty(input$customStartA),
                isEmpty(input$customEndA),
                as.integer(input$customStartA)>=as.integer(input$customEndA)
            ))
                areaExplorerMessages <- updateMessages(
                    areaExplorerMessages,
                    type="INFO",
                    msg=paste(getTime("INFO"),"Area to plot: ",customArea$name,
                        sep="")
                )
        })
        
        updateFlanksA <- eventReactive({
            input$upstreamFlankA
            input$downstreamFlankA
        },{
            flankValidate <- c(
                upstream=FALSE,
                downstream=FALSE
            )
            errMsg <- c(
                upstream="",
                downstream=""
            )
            u <- input$upstreamFlankA
            d <- input$downstreamFlankA
            if (u!=currentOpts$flank[1]) {
                u <- as.integer(u)
                if (is.na(suppressWarnings(u)) || isEmpty(u) || u<=0) {
                    flankValidate["upstream"] <- TRUE
                    errMsg["upstream"] <- paste("Upstream flanking region",
                        "must be a positive integer!")
                }
            }
            if (d!=currentOpts$flank[2]) {
                d <- as.integer(d)
                if (is.na(suppressWarnings(d)) || isEmpty(d) || d<=0) {
                    flankValidate["downstream"] <- TRUE
                    errMsg["downstream"] <- paste("Downstream flanking region",
                        "must be a positive integer!")
                }
            }
            if (any(flankValidate)) {
                output$areaFlankError <- renderUI({
                    div(class="error-message",paste(errMsg,collapse="\n"))
                })
            }
            else {
                output$areaFlankError <- renderUI({div()})
                if (u!=currentOpts$flank[1]) {
                    currentOpts$flank[1] <- u
                    areaExplorerMessages <- updateMessages(
                        areaExplorerMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),"Upstream flanking region",
                            " changed! New upstream flanking: ",u,sep="")
                    )
                }
                if (d!=currentOpts$flank[2]) {
                    currentOpts$flank[2] <- d
                    areaExplorerMessages <- updateMessages(
                        areaExplorerMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),"Downstream flanking ",
                            "region changed! New downstream flanking: ",d,
                            sep="")
                    )
                }
            }
        })
        
        addTheCustomRegionInArea <- eventReactive(input$addCustomRegionInArea,{
            regionValidate <- c(
                name=FALSE,
                strand=FALSE,
                chrom=FALSE,
                start=FALSE,
                end=FALSE,
                startBeforeEnd=FALSE,
                regionExistsName=FALSE,
                regionExistsLocus=FALSE
            )
            errMsg <- c(
                name="",
                strand="",
                chrom="",
                start="",
                end="",
                startBeforeEnd="",
                regionExistsName="",
                regionExistsLocus=""
            )
            if (input$customNameInArea=="") {
                regionValidate["name"] <- TRUE
                errMsg["name"] <- "A region name is required!"
            
            }
            if (input$customStrandInArea=="") {
                regionValidate["strand"] <- TRUE
                errMsg["strand"] <- "A region strand is required!"
            }
            if (input$customChrInArea=="") {
                regionValidate["chrom"] <- TRUE
                errMsg["chrom"] <- "A region chromosome is required!"
            
            }
            if (input$customStartInArea=="" 
                || as.integer(input$customStartInArea)<=0
                || is.na(suppressWarnings(
                    as.numeric(input$customStartInArea)))) {
                regionValidate["start"] <- TRUE
                errMsg["start"] <- 
                    "Region start base must be a valid positive integer!"
            }
            if (input$customEndInArea=="" 
                || as.numeric(input$customEndInArea)<=0
                || is.na(suppressWarnings(
                    as.integer(input$customEndInArea)))) {
                regionValidate["end"] <- TRUE
                errMsg["end"] <- 
                    "Region end base must be a valid positive integer!"
            }
            if (!(regionValidate["start"] || regionValidate["end"])) {
                # Check start before end
                s <- as.integer(input$customStartInArea)
                e <- as.integer(input$customEndInArea)
                if (e-s<=0) {
                    regionValidate["startBeforeEnd"] <- TRUE
                    errMsg["startBeforeEnd"] <- paste("Region start must be ",
                        "smaller than region end and the region length must ",
                        "be greater than 1!",sep="")
                }
                
                # Check the new region does not already exists
                cc <- customRegionsInArea$chromosome
                ss <- customRegionsInArea$start
                ee <- customRegionsInArea$end
                tt <- customRegionsInArea$strand
                if (input$customNameInArea %in% customRegionsInArea$name) {
                    regionValidate["regionExistsName"] <- TRUE
                    errMsg["regionExistsName"] <- paste("Custom regions must ",
                        "have each a unique name! ",input$customNameInArea," ",
                        "already exists.",sep="")
                }
                if (length(ss)>0) {
                    for (i in 1:length(ss)) {
                        if (input$customChrInArea==cc[i] 
                            && input$customStartInArea==ss[i] 
                            && input$customEndInArea==ee[i] 
                            && input$customStrandInArea==tt[i])
                        regionValidate["regionExistsLocus"] <- TRUE
                        errMsg["regionExistsLocus"] <- 
                            "Custom regions must have unique coordinates!"
                    }
                }
            }
            
            if (any(regionValidate)) {
                output$customRegionErrorInArea <- renderUI({
                    div(class="error-message",paste(errMsg,collapse="\n"))
                })
                return("")
            }
            else {
                output$customRegionErrorInArea <- renderUI({div()})
                customRegionsInArea$chromosome <- c(customRegionsInArea$chr,
                    input$customChrInArea)
                customRegionsInArea$start <- c(customRegionsInArea$start,
                    input$customStartInArea)
                customRegionsInArea$end <- c(customRegionsInArea$end,
                    input$customEndInArea)
                customRegionsInArea$strand <- c(customRegionsInArea$strand,
                    input$customStrandInArea)
                customRegionsInArea$name <- c(customRegionsInArea$name,
                    input$customNameInArea)
                
                areaExplorerMessages <- updateMessages(
                    areaExplorerMessages,
                    type="SUCCESS",
                    msg=paste(getTime("SUCCESS"),"Custom region ",
                        input$customNameInArea,
                        " added to the search of custom area to plot! Number",
                        " of custom regions to plot is now ",
                        length(customRegionsInArea$name),sep="")
                )
            }
        })
        
        removeTheCustomRegionInArea <- eventReactive(
            input$removeCustomRegionInArea,{
            output$customRegionErrorInArea <- renderUI({div()})
            j <- as.integer(input$customRegionListInArea_rows_selected)
            if (length(j)>0) {
                ii <- customRegionsInArea$name[j]
                customRegionsInArea$chromosome <- 
                    customRegionsInArea$chromosome[-j]
                customRegionsInArea$start <- customRegionsInArea$start[-j]
                customRegionsInArea$end <- customRegionsInArea$end[-j]
                customRegionsInArea$strand <- customRegionsInArea$strand[-j]
                customRegionsInArea$name <- customRegionsInArea$name[-j]
                areaExplorerMessages <- updateMessages(
                    areaExplorerMessages,
                    type="WARNING",
                    msg=paste(getTime("WARNING"),"Custom region ",ii," ",
                        " removed from the search of custom area to plot! ",
                        "Number of custom regions to plot is now ",
                        length(customRegionsInArea$name),sep="")
                )
            }
        })
        
        createAreaProfile <- eventReactive(input$createAreaProfile,{
            if (is.null(currentMetadata$final)) {
                output$areaExplorerError <- renderUI({
                    div(class="error-message",paste("You must create a ",
                        "dataset first!",sep=""))
                })
            }
            else {
                if (!is.na(currentMetadata$genome) &&
                    is.null(loadedGenomes[[currentMetadata$genome]]$dbExon)) {
                    load(file.path("genome",currentMetadata$genome,
                        "summarized_exon.rda"))
                    loadedGenomes[[currentMetadata$genome]]$dbExon <<- sexon
                    dataSelectorMessages <- updateMessages(
                        dataSelectorMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),"Genome ",
                            currentMetadata$genome," exons loaded!", sep="")
                    )
                }
                isolate(input$areaTypeRadio)
                output$areaExplorerError <- renderUI({div()})
                dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
                dbExon <- loadedGenomes[[currentMetadata$genome]]$dbExon
                if (input$areaTypeRadio=="gene") {
                    gn <- as.character(input$areaGeneName)
                    fs <- as.integer(input$upstreamFlankA)
                    fd <- as.integer(input$downstreamFlankA)
                    dbgi <- dbGene[gn]
                    refArea <- list(
                        chr=as.character(seqnames(dbgi))[1],
                        start=start(dbgi)-fs,
                        end=end(dbgi)+fd
                    )
                    customTrans <- NULL
                }
                else {
                    refArea <- list(
                        chr=as.character(customArea$chromosome),
                        start=as.integer(customArea$start),
                        end=as.integer(customArea$end)
                    )
                    if (input$customRegionsFromGeneExplorer=="fromge") {
                        customRegionsInArea <- customRegions
                        ii <- which(customRegions$chromosome==refArea$chr)
                        if (length(ii)>0) {
                            customRegionsInArea$chromosome <- 
                                customRegionsInArea$chromosome[ii]
                            customRegionsInArea$start <- 
                                customRegionsInArea$start[ii]
                            customRegionsInArea$end <- 
                                customRegionsInArea$end[ii]
                            customRegionsInArea$strand <- 
                                customRegionsInArea$strand[ii]
                            customRegionsInArea$name <- 
                                customRegionsInArea$name[ii]
                        }
                        else {
                            customRegionsInArea$chromosome <- character(0)
                            customRegionsInArea$start <- integer(0)
                            customRegionsInArea$end <- integer(0)
                            customRegionsInArea$strand <- character(0)
                            customRegionsInArea$name <- character(0)
                        }
                    }
                    if (length(customRegionsInArea$chromosome)>0) {
                        customTrans <- as.list(makeGRangesFromDataFrame(
                            df=data.frame(
                                chromosome=
                                    as.character(
                                        customRegionsInArea$chromosome),
                                start=as.integer(customRegionsInArea$start),
                                end=as.integer(customRegionsInArea$end),
                                strand=as.character(customRegionsInArea$strand),
                                gene_id=as.character(customRegionsInArea$name)
                            ),
                            keep.extra.columns=TRUE
                        ))
                        names(customTrans) <- as.character(customRegions$name)
                    }
                    else
                        customTrans <- NULL
                }
                
                progress <- shiny::Progress$new()
                progress$initialize(
                    session,
                    min=0,
                    max=length(unique(as.character(
                        currentMetadata$finalclass)))*2+1
                )
                progress$set(message="",value=0)
                on.exit(progress$close())
                
                updateProgress <- function(value=NULL,detail=NULL) {
                    if (is.null(value)) {
                        value <- progress$getValue()
                        value <- value + 1
                    }
                    progress$set(value=value,detail=detail)
                }
                
                ggTrack <- tryCatch(
                    getTrack(
                        refArea=refArea,
                        customGene=customTrans,
                        source=as.character(currentMetadata$source),
                        dataset=as.character(currentMetadata$dataset),
                        #class=as.character(currentMetadata$class),
                        class=unique(as.character(currentMetadata$final$class)),
                        sumStat=as.character(currentOpts$sumStat),
                        config=currentMetadata$final,
                        dbGene=dbGene,
                        dbExon=dbExon,
                        trim=as.numeric(currentOpts$trim),
                        classColours=currentOpts$colours,
                        messageContainer=areaExplorerMessages,
                        progressFun=updateProgress,
                        rc=RC
                    ),
                    error=function(e) {
                        areaExplorerMessages <- updateMessages(
                            areaExplorerMessages,
                            type="WARNING",
                            msg=paste(getTime("WARNING"),"An error occured ",
                                "while calculating profiles from bigWig ",
                                "files. Will try with BAM if present. Please ",
                                "contact the administrator stating the ",
                                "following error: ",e,sep="")
                        )
                        tryCatch(
                            getTrack(
                                refArea=refArea,
                                customGene=customTrans,
                                source=as.character(currentMetadata$source),
                                dataset=as.character(currentMetadata$dataset),
                                #class=as.character(currentMetadata$class),
                                class=unique(as.character(
                                    currentMetadata$final$class)),
                                sumStat=as.character(currentOpts$sumStat),
                                config=currentMetadata$final,
                                dbGene=dbGene,
                                dbExon=dbExon,
                                trim=as.numeric(currentOpts$trim),
                                fromBam=TRUE,
                                classColours=currentOpts$colours,
                                messageContainer=areaExplorerMessages,
                                progressFun=updateProgress,
                                rc=RC
                            ),
                            error=function(e) {
                                areaExplorerMessages <- updateMessages(
                                    areaExplorerMessages,
                                    type="ERROR",
                                    msg=paste(getTime("ERROR"),"An error ",
                                        "occured while calculating profiles ",
                                        "from BAM files. Please contact the ",
                                        "administrator stating the following ",
                                        "error: ",e,sep="")
                                )
                            },
                            finally=""
                        )
                    },finally=""
                )
                
                ggTrack <- ggTrack + 
                    scale_color_manual(values=currentOpts$colours) + 
                    scale_fill_manual(values=currentOpts$colours)
                
                areaPlots$areaProfile <- ggTrack
                areaPlots$rendered <- TRUE
            }
        })
        
        output$customRegionListInArea <- DT::renderDataTable(
            data.frame(
                chromosome=customRegionsInArea$chromosome,
                start=customRegionsInArea$start,
                end=customRegionsInArea$end,
                strand=customRegionsInArea$strand,
                name=customRegionsInArea$name
            ),
            class="display compact",
            rownames=FALSE,
            options=list(
                dom='t',
                paging=FALSE
            )
        )
        
        output$setChrsA <- renderUI({
            selectizeInput(
                inputId="customChrA",
                label="", 
                choices=c("",
                    getValidChromosomes("hg19")
                ),
                options=list(
                    placeholder="Chrom"
                )
            )
        })
        
        output$areaExplorerColours <- renderUI({
            if (!is.null(currentMetadata$final)) {
                #c <- currentMetadata$class
                c <- unique(as.character(currentMetadata$final$class))
                lapply(1:length(c),function(i,c) {
                    colourInput(
                        inputId=paste("areaColour_",c[i],sep=""),
                        label=paste("Select colour for",c[i]),
                        value=baseColours[i]
                    )
                },c)
            }
        })
        
        output$areaProfile <- renderPlot({
            areaPlots$areaProfile
        })
        
        output$exportAreaPDF <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("area_plot_",tt,".pdf", sep='')
            },
            content=function(con) {
                ggsave(filename=con,plot=areaPlots$areaProfile,
                    width=14,height=7)
            }
        )
        
        output$exportAreaPNG <- downloadHandler(
            filename=function() {
                tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                paste("area_plot_",tt,".png", sep='')
            },
            content=function(con) {
                ggsave(filename=con,plot=areaPlots$areaProfile,
                    width=14,height=7)
            }
        )
        
        output$areaExplorerMessages <- renderUI({
            lapply(areaExplorerMessages$messages,function(x) {
                switch(x$type,
                    INFO = {
                        div(class="info-box",x$msg)
                    },
                    SUCCESS = {
                        div(class="success-box",x$msg)
                    },
                    WARNING = {
                        div(class="warn-box",x$msg)
                    },
                    ERROR = {
                        div(class="error-box",x$msg)
                    }
                )
            })
        })
        
        isolate({
            g <- input$areaGeneName
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"areaGeneName",
                    choices=geneNames[i],
                    selected=input$areaGeneName,
                    server=TRUE)
            }
        })
        
        observe({
            updateTextInput(session,"customChrInArea",label="",
                value=input$customChrA)
        })
        
        observe({
            if (length(customRegionsInArea$name)==0)
                shinyjs::disable("removeCustomRegionInArea")
            else
                shinyjs::enable("removeCustomRegionInArea")
        })
        
        observe({
            customRegionsInArea <- addTheCustomRegionInArea()
            customRegionsInArea <- removeTheCustomRegionInArea()
        })
        
        # Trimming textbox
        observe({
            if (input$areaSumStatType!="trimmed") {
                shinyjs::disable("areaTrimPct")
                shinyjs::enable("createAreaProfile")
            }
            else {
                shinyjs::enable("areaTrimPct")
                trimp <- as.numeric(input$areaTrimPct)
                if (trimp<0 || trimp>0.5) {
                    output$areaExplorerError <- renderUI({
                        div(class="error-message",paste("The trimming ",
                            "must be between 0 and 0.5!",sep=""))
                    })
                    shinyjs::disable("createAreaProfile")
                }
                else {
                    output$areaExplorerError <- renderUI({div()})
                    if (isEmpty(input$areaGeneName) 
                        && length(customArea$name)==0)
                        shinyjs::disable("createAreaProfile")
                    else
                        shinyjs::enable("createAreaProfile")
                    currentOpts$trim <- trimp
                }
            }
        })
        
        observe({
            # Engage button status
            customArea$chromosome <- input$customChrA
            customArea$start <- input$customStartA
            customArea$end <- input$customEndA
            if (any(is.null(customArea$chromosome) || customArea$chromosome=="",
                is.null(customArea$start) || customArea$start=="",
                is.null(customArea$end) || customArea$end=="",
                customArea$start>=customArea$end)
                && input$areaTypeRadio=="area")
                shinyjs::disable("createAreaProfile")
            else
                shinyjs::enable("createAreaProfile")
        })
        
        observe({
            tryCatch({
                shinyjs::disable("createAreaProfile")
                createAreaProfile()
            },error=function(e) {
                areaPlots$rendered <- FALSE
                return()
            },
            finally={
                if (isEmpty(input$areaGeneName) 
                    && length(customArea$name)==0)
                    shinyjs::disable("createAreaProfile")
                else
                    shinyjs::enable("createAreaProfile")
            })
        })
        
        observe({
            updateFlanksA()
            updateCurrentArea()
            updateAreaColours()
            updateAreaExplorerMessages()
        })
        
        ########################################################################
        
        # Fifth main panel A - Expression explorer
        
        currentTables <- reactiveValues()
        
        createCurrentTables <- reactive({
            for (c in currentMetadata$class)
                currentTables[[c]] <- NULL
        })
        
        countDataTables <- reactive({
            dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
            M <- currentMetadata$final
            s <- unique(as.character(M$source))
            d <- unique(as.character(M$dataset))
            c <- unique(as.character(M$class))
            lapply(c,function(x,M,D) {
                samples <- as.character(M[which(as.character(M$class)==x),
                    "sample_id"])
                switch(input$rnaExpressionGeneList,
                    all = {
                        switch(input$rnaExpressionMeasureRadio,
                            raw = {
                                tab <- D$counts[,samples,drop=FALSE]
                                average <- round(apply(tab,1,
                                    input$rnaExpressionAverageRadio))
                                deviation <- round(apply(tab,1,
                                    input$rnaExpressionDeviationRadio))
                            },
                            norm = {
                                tab <- D$norm[,samples,drop=FALSE]
                                average <- round(apply(tab,1,
                                    input$rnaExpressionAverageRadio))
                                deviation <- round(apply(tab,1,
                                    input$rnaExpressionDeviationRadio))
                            },
                            rpkm = {
                                tab <- round(edgeR::rpkm(
                                    D$counts[,samples,drop=FALSE],
                                    gene.length=D$length,
                                    lib.size=unlist(D$libsize[samples]),
                                    ),digits=6)
                                average <- round(apply(tab,1,
                                    input$rnaExpressionAverageRadio),digits=6)
                                deviation <- round(apply(tab,1,
                                    input$rnaExpressionDeviationRadio),digits=6)
                            },
                            rpgm = {
                                tab <- round(D$norm[,samples,
                                    drop=FALSE]/D$length,digits=6)
                                average <- round(apply(tab,1,
                                    input$rnaExpressionAverageRadio),digits=6)
                                deviation <- round(apply(tab,1,
                                    input$rnaExpressionDeviationRadio),digits=6)
                            }
                        )
                    },
                    select = {
                        g <- input$selectExpressionGeneName
                        if (isEmpty(g))
                            tab <- data.frame(
                                name=character(0),
                                average=numeric(0),
                                deviation=numeric(0),
                                sample=character(0)
                            )
                        else {
                            switch(input$rnaExpressionMeasureRadio,
                                raw = {
                                    tab <- D$counts[g,samples,drop=FALSE]
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio))
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio))
                                },
                                norm = {
                                    tab <- D$norm[g,samples,drop=FALSE]
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio))
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio))
                                },
                                rpkm = {
                                    tab <- round(edgeR::rpkm(
                                        D$counts[g,samples,drop=FALSE],
                                        gene.length=D$length[g],
                                        lib.size=unlist(D$libsize[samples])
                                        ),digits=6)
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio),
                                            digits=6)
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio),
                                            digits=6)
                                },
                                rpgm = {
                                    tab <- round(D$norm[g,
                                        samples,drop=FALSE]/D$length[g],
                                            digits=6)
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio),
                                            digits=6)
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio),
                                            digits=6)
                                }
                            )
                        }
                    },
                    custom = {
                        g <- input$rnaExpressionCustomList
                        if (isEmpty(g))
                            tab <- data.frame(
                                name=character(0),
                                average=numeric(0),
                                deviation=numeric(0),
                                sample=character(0)
                            )
                        else {
                            g <- strsplit(g,split="\n")[[1]]
                            m <- match(g,as.character(dbGene$gene_name))
                            na <- which(is.na(m))
                            if (length(na)>0)
                                m <- m[-na]
                            gene_id <- names(dbGene)[m]
                            switch(input$rnaExpressionMeasureRadio,
                                raw = {
                                    tab <- D$counts[gene_id,samples,drop=FALSE]
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio))
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio))
                                },
                                norm = {
                                    tab <- D$norm[gene_id,samples,drop=FALSE]
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio))
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio))
                                },
                                rpkm = {
                                    tab <- round(edgeR::rpkm(
                                        D$counts[gene_id,samples,drop=FALSE],
                                        gene.length=D$length[gene_id],
                                        lib.size=unlist(D$libsize[samples])
                                        ),digits=6)
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio),digits=6)
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio),digits=6)
                                },
                                rpgm = {
                                    tab <- round(D$norm[gene_id,
                                        samples,drop=FALSE]/D$length[gene_id],
                                        digits=6)
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio),digits=6)
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio),digits=6)
                                }
                            )
                        }
                    }
                )
                if (!isEmpty(g)) {
                    switch(input$rnaExpressionScaleRadio,
                        natural = {
                            tab <- tab
                            average <- average
                            deviation <- deviation
                        },
                        log2 = {
                            tab <- round(log2(tab+1),digits=6)
                            average <- round(apply(tab,1,
                                input$rnaExpressionAverageRadio),digits=6)
                            deviation <- round(apply(tab,1,
                                input$rnaExpressionDeviationRadio),digits=6)
                        }
                    )
                    
                    # Use alternative ids if present
                    if (!is.null(M$alt_id))
                        colnames(tab) <- as.character(M[which(
                            as.character(M$class)==x),"alt_id"])
                    # Attach some annotation
                    meta <- data.frame(
                        #chromosome=as.character(seqnames(dbGene[rownames(tab)])),
                        #start=start(dbGene[rownames(tab)]),
                        #end=end(dbGene[rownames(tab)]),
                        name=as.character(dbGene[rownames(tab)]$gene_name),
                        average=average,
                        deviation=deviation
                    )
                    currentTables[[x]] <- cbind(meta,tab)
                }
                else 
                    meta <- NULL
                output[[paste("countTable",x,sep="_")]] <- 
                    DT::renderDataTable(
                        cbind(meta,tab),
                        class="display compact",
                        rownames=FALSE,
                        options=list(
                            searchHighlight=TRUE,
                            pageLength=10,
                            lengthMenu=c(10,20,50,100)
                        )
                    )
            },M,loadedData[[s]][[d]])
        })
        
        handleExpressionSampleSelection <- reactive({
            c <- currentMetadata$class
            lapply(c,function(x) {
                observeEvent(input[[paste("clearCountSelection_",x,sep="")]],{
                    proxy <- dataTableProxy(paste("countTable_",x,sep=""))
                    selectRows(proxy,NULL)
                })
            })
            lapply(c,function(x) {
                observeEvent(input[[paste("invertCountSelection_",x,sep="")]],{
                    N <- input[[paste("countTable_",x,"_rows_all",
                        sep="")]]
                    sel <- input[[paste("countTable_",x,"_rows_selected",
                        sep="")]]
                    if (length(sel)>0) {
                        N <- N[-sel]
                        proxy <- dataTableProxy(paste("countTable_",x,sep=""))
                        selectRows(proxy,N)
                    }
                })
            })
        })
        
        handleExpressionDownloadSelection <- reactive({
            c <- currentMetadata$class
            lapply(c,function(x) {
                output[[paste("exportCountSelection_",x,sep="")]] <- 
                    downloadHandler(
                        filename=function() {
                            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                            paste(x,"_gene_expression_",tt,".txt", sep='')
                        },
                        content=function(con) {
                            sel <- input[[paste("countTable_",x,
                                "_rows_selected",sep="")]]
                            if (length(sel)>0)
                                write.table(currentTables[[x]][sel,],file=con,
                                    sep="\t",quote=FALSE,col.names=NA)
                        }
                    )
                output[[paste("exportCountAll_",x,sep="")]] <- downloadHandler(
                    filename=function() {
                        tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                        paste(x,"_gene_expression_",tt,".txt", sep='')
                    },
                    content=function(con) {
                        write.table(currentTables[[x]],file=con,sep="\t",
                                quote=FALSE,col.names=NA)
                    }
                )
            })
        })
        
        output$rnaExpressionTables <- renderUI({
            M <- currentMetadata$final
            if (is.null(M))
                div(
                    style="display:inline-block; margin:5px;",
                    h4(paste("Please create a dataset first using the",
                        "'Data selector' tab."))
                )
            else {
                s <- unique(as.character(M$source))
                d <- unique(as.character(M$dataset))
                c <- unique(as.character(M$class))
                countDataTables()
                lapply(c,function(x,s,d) {
                    list(
                        div(
                            class="small table-container",
                            h4(paste("Gene expression for ",x)),
                            DT::dataTableOutput(paste("countTable",x,sep="_")),
                            div(
                                style="display:inline-block; margin:5px;",
                                actionButton(
                                    paste("clearCountSelection_",x,sep=""),
                                    paste("Clear selection"),
                                    class="btn-xs"
                                )
                            ),
                            div(
                                style="display:inline-block; margin:5px;",
                                actionButton(
                                    paste("invertCountSelection_",x,sep=""),
                                    paste("Invert selection"),
                                    class="btn-xs"
                                )
                            ),
                            div(
                                style="display:inline-block; margin:5px;",
                                downloadButton(
                                    paste("exportCountSelection_",x,sep=""),
                                    paste("Export selection"),
                                    class="btn-xs"
                                )
                            ),
                            div(
                                style="display:inline-block; margin:5px;",
                                downloadButton(
                                    paste("exportCountAll_",x,sep=""),
                                    paste("Export all"),
                                    class="btn-xs"
                                )
                            )
                        ),
                        hr()
                    )
                },s,d)
            }
        })
        
        isolate({
            g <- input$selectExpressionGeneName
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"selectExpressionGeneName",
                    choices=geneNames[i],
                    selected=input$selectExpressionGeneName,
                    server=TRUE)
            }
        })
        
        observe({
            createCurrentTables()
            handleExpressionSampleSelection()
            handleExpressionDownloadSelection()
        })
        
        # Fifth main panel B - Custom RNA calculator
        
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
        
        createCurrentCustomRnaTables <- reactive({
            for (c in currentMetadata$class)
                currentCustomRnaTables$tables[[c]] <- NULL
        })
        
        addTheCustomRnaRegion <- eventReactive(input$addRnaCustomRegion,{
            regionValidate <- c(
                name=FALSE,
                strand=FALSE,
                chrom=FALSE,
                start=FALSE,
                end=FALSE,
                startBeforeEnd=FALSE,
                regionExistsName=FALSE,
                regionExistsLocus=FALSE
            )
            errMsg <- c(
                name="",
                strand="",
                chrom="",
                start="",
                end="",
                startBeforeEnd="",
                regionExistsName="",
                regionExistsLocus=""
            )
            if (input$customNameExpr=="") {
                regionValidate["name"] <- TRUE
                errMsg["name"] <- "A region name is required!"
            
            }
            if (input$customStrandExpr=="") {
                regionValidate["strand"] <- TRUE
                errMsg["strand"] <- "A region strand is required!"
            }
            if (input$customChrExpr=="") {
                regionValidate["chrom"] <- TRUE
                errMsg["chrom"] <- "A region chromosome is required!"
            
            }
            if (input$customStartExpr=="" 
                || as.integer(input$customStartExpr)<=0
                || is.na(suppressWarnings(
                    as.numeric(input$customStartExpr)))) {
                regionValidate["start"] <- TRUE
                errMsg["start"] <- 
                    "Region start base must be a valid positive integer!"
            }
            if (input$customEndExpr=="" 
                || as.numeric(input$customEndExpr)<=0
                || is.na(suppressWarnings(
                    as.integer(input$customEndExpr)))) {
                regionValidate["end"] <- TRUE
                errMsg["end"] <- 
                    "Region end base must be a valid positive integer!"
            }
            if (!(regionValidate["start"] || regionValidate["end"])) {
                # Check start before end
                s <- as.integer(input$customStartExpr)
                e <- as.integer(input$customEndExpr)
                if (e-s<=0) {
                    regionValidate["startBeforeEnd"] <- TRUE
                    errMsg["startBeforeEnd"] <- paste("Region start must be ",
                        "smaller than region end and the region length must ",
                        "be greater than 1!",sep="")
                }
                
                # Check the new region does not already exists
                cc <- customRnaRegions$chromosome
                ss <- customRnaRegions$start
                ee <- customRnaRegions$end
                tt <- customRnaRegions$strand
                if (input$customNameExpr %in% customRnaRegions$name) {
                    regionValidate["regionExistsName"] <- TRUE
                    errMsg["regionExistsName"] <- paste("Custom regions must ",
                        "have each a unique name! ",input$customNameExpr," ",
                        "already exists.",sep="")
                }
                if (length(ss)>0) {
                    for (i in 1:length(ss)) {
                        if (input$customChrExpr==cc[i] 
                            && input$customStartExpr==ss[i] 
                            && input$customEndExpr==ee[i] 
                            && input$customStrandExpr==tt[i])
                        regionValidate["regionExistsLocus"] <- TRUE
                        errMsg["regionExistsLocus"] <- 
                            "Custom regions must have unique coordinates!"
                    }
                }
            }
            
            if (any(regionValidate)) {
                output$customRegionExprError <- renderUI({
                    div(class="error-message",paste(errMsg,collapse="\n"))
                })
                return("")
            }
            else {
                output$customRegionExprError <- renderUI({div()})
                customRnaRegions$chromosome <- c(customRnaRegions$chromosome,
                    input$customChrExpr)
                customRnaRegions$start <- c(customRnaRegions$start,
                    input$customStartExpr)
                customRnaRegions$end <- c(customRnaRegions$end,
                    input$customEndExpr)
                customRnaRegions$strand <- c(customRnaRegions$strand,
                    input$customStrandExpr)
                customRnaRegions$name <- c(customRnaRegions$name,
                    input$customNameExpr)
            }
        })
        
        removeTheCustomRnaRegion <- eventReactive(
            input$removeRnaCustomRegion,{
            output$customRegionExprError <- renderUI({div()})
            j <- as.integer(input$customRegionListExpr_rows_selected)
            if (length(j)>0) {
                ii <- customRnaRegions$name[j]
                customRnaRegions$chromosome <- 
                    customRnaRegions$chromosome[-j]
                customRnaRegions$start <- customRnaRegions$start[-j]
                customRnaRegions$end <- customRnaRegions$end[-j]
                customRnaRegions$strand <- customRnaRegions$strand[-j]
                customRnaRegions$name <- customRnaRegions$name[-j]
            }
        })
        
        calcCustomRnaCounts <- eventReactive(input$calculateCustomRegionRna,{
            if (is.null(currentMetadata$final)) {
                output$customRnaCalcError <- renderUI({
                    div(class="error-message",paste("You must create a ",
                        "dataset first!",sep=""))
                })
            }
            else {
                output$customRnaCalcError <- renderUI({div()})
                if (input$customExpressionFromGeneExplorer=="fromge")
                    customRnaRegions <- customRegions
                customInputRegions <- makeGRangesFromDataFrame(
                    df=data.frame(
                        chromosome=as.character(customRnaRegions$chromosome),
                        start=as.numeric(customRnaRegions$start),
                        end=as.numeric(customRnaRegions$end),
                        strand=as.character(customRnaRegions$strand),
                        gene_id=as.character(customRnaRegions$name)
                    ),
                    keep.extra.columns=TRUE
                )
                currentCustomRnaTables$lengths <- width(customInputRegions)
                names(currentCustomRnaTables$lengths) <- 
                    names(customInputRegions)
                #customInputRegions <- as.list(customInputRegions)
                names(customInputRegions) <- as.character(customRnaRegions$name)
                progress <- shiny::Progress$new()
                progress$initialize(
                    session,
                    min=0,
                    max=length(customInputRegions)
                )
                progress$set(message="",value=0)
                on.exit(progress$close())
                
                updateProgress <- function(value=NULL,detail=NULL) {
                    if (is.null(value)) {
                        value <- progress$getValue()
                        value <- value + 1
                    }
                    progress$set(value=value,detail=detail)
                }
                
                customCounts <- tryCatch({
                    getCustomCounts(
                        coords=customInputRegions,
                        samples=as.character(currentMetadata$final$sample_id),
                        config=currentMetadata$final,
                        progressFun=updateProgress,
                        rc=RC
                    )},
                    error=function(e) {
                    },finally=""
                )
                
                for (c in unique(as.character(currentMetadata$final$class))) {
                    sc <- as.character(currentMetadata$final[
                        currentMetadata$final$class==c,"sample_id"])
                    currentCustomRnaTables$tables[[c]] <- 
                        customCounts[,sc,drop=FALSE]
                }
            }
        })
        
        customCountRnaDataTables <- reactive({
            M <- currentMetadata$final
            s <- unique(as.character(M$source))
            d <- unique(as.character(M$dataset))
            c <- unique(as.character(M$class))
            le <- currentCustomRnaTables$lengths
            lapply(c,function(x,M,le,D) {
                tab <- currentCustomRnaTables$tables[[x]]
                if (is.null(tab)) {
                    tab <- data.frame(
                        name=character(0),
                        average=numeric(0),
                        deviation=numeric(0),
                        sample=character(0)
                    )
                    meta <- NULL
                }
                else {
                    tab <- as.matrix(tab)
                    switch(input$rnaCustomMeasureRadio,
                        raw = {
                            average <- round(apply(tab,1,
                                input$rnaCustomAverageRadio))
                            deviation <- round(apply(tab,1,
                                input$rnaCustomDeviationRadio))
                        },
                        #norm = {
                        #   average <- round(apply(tab,1,
                        #       input$rnaCustomAverageRadio))
                        #   deviation <- round(apply(tab,1,
                        #       input$rnaCustomDeviationRadio))
                        #},
                        rpkm = {
                            tab <- round(edgeR::rpkm(tab,gene.length=le,
                                lib.size=unlist(D$libsize[colnames(tab)]),),
                                digits=6)
                            average <- round(apply(tab,1,
                                input$rnaCustomAverageRadio),digits=6)
                            deviation <- round(apply(tab,1,
                                input$rnaCustomDeviationRadio),digits=6)
                        },
                        rpgm = {
                            for (ss in colnames(tab))
                                tab[,ss] <- round(tab[,ss]/le,digits=6)
                            average <- round(apply(tab,1,
                                input$rnaCustomAverageRadio),digits=6)
                            deviation <- round(apply(tab,1,
                                input$rnaCustomDeviationRadio),digits=6)
                        }
                    )
                    if (!isEmpty(tab)) {
                        switch(input$rnaCustomScaleRadio,
                            natural = {
                                tab <- tab
                                average <- average
                                deviation <- deviation
                            },
                            log2 = {
                                tab <- round(log2(tab+1),digits=6)
                                average <- round(apply(tab,1,
                                    input$rnaCustomAverageRadio),digits=6)
                                deviation <- round(apply(tab,1,
                                    input$rnaCustomDeviationRadio),digits=6)
                            }
                        )
                        # Use alternative ids if present
                        if (!is.null(M$alt_id))
                            colnames(tab) <- as.character(M[which(
                                as.character(M$class)==x),"alt_id"])
                        # Attach some annotation
                        meta <- data.frame(
                            name=rownames(tab),
                            average=average,
                            deviation=deviation
                        )
                    }
                    else 
                        meta <- NULL
                }
                output[[paste("customCountRnaTable",x,sep="_")]] <- 
                    DT::renderDataTable(
                        cbind(meta,tab),
                        class="display compact",
                        rownames=FALSE,
                        options=list(
                            searchHighlight=TRUE,
                            pageLength=10,
                            lengthMenu=c(10,20,50,100)
                        )
                    )
            },M,le,loadedData[[s]][[d]])
        })
        
        handleCustomRnaSampleSelection <- reactive({
            c <- currentMetadata$class
            lapply(c,function(x) {
                observeEvent(input[[paste("clearCustomRnaCountSelection_",x,
                    sep="")]],{
                    proxy <- dataTableProxy(paste("customCountRnaTable_",x,
                        sep=""))
                    selectRows(proxy,NULL)
                })
            })
            lapply(c,function(x) {
                observeEvent(input[[paste("invertCustomRnaCountSelection_",x,
                    sep="")]],{
                    N <- input[[paste("customCountRnaTable_",x,"_rows_all",
                        sep="")]]
                    sel <- input[[paste("customCountRnaTable_",x,
                        "_rows_selected",sep="")]]
                    if (length(sel)>0) {
                        N <- N[-sel]
                        proxy <- dataTableProxy(paste("customCountRnaTable_",x,
                            sep=""))
                        selectRows(proxy,N)
                    }
                })
            })
        })
        
        handleCustomRnaDownloadSelection <- reactive({
            c <- currentMetadata$class
            lapply(c,function(x) {
                output[[paste("exportCustomRnaCountSelection_",x,sep="")]] <- 
                    downloadHandler(
                        filename=function() {
                            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                            paste(x,"_custom_expression_",tt,".txt", sep='')
                        },
                        content=function(con) {
                            sel <- input[[paste("customCountRnaTable_",x,
                                "_rows_selected",sep="")]]
                            if (length(sel)>0)
                                write.table(currentCustomRnaTables$tables[[
                                    x]][sel,],file=con,sep="\t",quote=FALSE,
                                    col.names=NA)
                        }
                    )
                output[[paste("exportCustomRnaCountAll_",x,sep="")]] <- 
                    downloadHandler(
                        filename=function() {
                            tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                            paste(x,"_custom_expression_",tt,".txt", sep='')
                        },
                        content=function(con) {
                            write.table(currentCustomRnaTables$tables[[x]],
                                file=con,sep="\t",quote=FALSE,col.names=NA)
                        }
                    )
            })
        })
        
        output$customRegionListExpr <- DT::renderDataTable(
            data.frame(
                chromosome=customRnaRegions$chromosome,
                start=customRnaRegions$start,
                end=customRnaRegions$end,
                strand=customRnaRegions$strand,
                name=customRnaRegions$name
            ),
            class="display compact",
            rownames=FALSE,
            options=list(
                dom='t',
                paging=FALSE
            )
        )
        
        output$rnaCustomTables <- renderUI({
            M <- currentMetadata$final
            if (is.null(M))
                div(
                    style="display:inline-block; margin:5px;",
                    h4(paste("Please create a dataset first using the",
                        "'Data selector' tab."))
                )
            else {
                s <- unique(as.character(M$source))
                d <- unique(as.character(M$dataset))
                c <- unique(as.character(M$class))
                customCountRnaDataTables()
                lapply(c,function(x,s,d) {
                    list(
                        div(
                            class="small table-container",
                            h4(paste("Expression for ",x)),
                            DT::dataTableOutput(paste("customCountRnaTable",x,
                                sep="_")),
                            div(
                                style="display:inline-block; margin:5px;",
                                actionButton(
                                    paste("clearCustomRnaCountSelection_",x,
                                        sep=""),
                                    paste("Clear selection"),
                                    class="btn-xs"
                                )
                            ),
                            div(
                                style="display:inline-block; margin:5px;",
                                actionButton(
                                    paste("invertCustomRnaCountSelection_",x,
                                        sep=""),
                                    paste("Invert selection"),
                                    class="btn-xs"
                                )
                            ),
                            div(
                                style="display:inline-block; margin:5px;",
                                downloadButton(
                                    paste("exportCustomRnaCountSelection_",x,
                                        sep=""),
                                    paste("Export selection"),
                                    class="btn-xs"
                                )
                            ),
                            div(
                                style="display:inline-block; margin:5px;",
                                downloadButton(
                                    paste("exportCustomRnaCountAll_",x,sep=""),
                                    paste("Export all"),
                                    class="btn-xs"
                                )
                            )
                        ),
                        hr()
                    )
                },s,d)
            }
        })
        
        output$setChrsExpr <- renderUI({
            selectizeInput(
                inputId="customChrExpr",
                label="", 
                choices=c("",
                    getValidChromosomes("hg19")
                ),
                options=list(
                    placeholder="Chrom"
                )
            )
        })
        
        observe({
            customRnaRegions <- addTheCustomRnaRegion()
            customRnaRegions <- removeTheCustomRnaRegion()
        })
        
        observe({
            customCountRnaDataTables()
            handleCustomRnaSampleSelection()
            handleCustomRnaDownloadSelection()
        })
        
        observe({
            if (length(customRnaRegions$name)!=0 
                || (input$customExpressionFromGeneExplorer=="fromge"
                    && length(customRegions$name)!=0))
                shinyjs::enable("calculateCustomRegionRna")
            else
                shinyjs::disable("calculateCustomRegionRna")
        })
        
        observe({
            if (length(customRnaRegions$name)==0)
                shinyjs::disable("removeRnaCustomRegion")
            else
                shinyjs::enable("removeRnaCustomRegion")
        })
        
        observe({
            tryCatch({
                shinyjs::disable("calculateCustomRegionRna")
                calcCustomRnaCounts()
            },error=function(e) {
                return()
            },
            finally={
                if (length(customRnaRegions$name)==0)
                    shinyjs::disable("calculateCustomRegionRna")
                else
                    shinyjs::enable("calculateCustomRegionRna")
            })
        })
        
        ########################################################################
        
        # Sixth main panel - Genome browser
        
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

loadJBrowse <- function(source,dataset,config,org="hg19") {
    urlBase <- "http://epigenomics.fleming.gr/bigseqcbrowse/index.html?"
    
    tracksBase <- paste("http://epigenomics.fleming.gr/bigseqcvis_tracks",
        org,sep="/")
    
    ind <- which(as.character(config$source)==source 
        & as.character(config$dataset)==dataset)
    subconf <- config[ind,]
    
    if (!is.null(subconf$alt_id))
        initTracks <- paste(as.character(subconf$sample_id),
            as.character(subconf$alt_id),"xy",sep="_")
    else
        initTracks <- paste(as.character(subconf$sample_id),"xy",sep="_")
    
    initTracks <- c(initTracks,"ucsc_human_hg19","ensGene","refGene")
    initTracks <- paste(initTracks,collapse=",")
    
    query <- paste("data=",tracksBase,"&tracks=",initTracks,sep="");
    
    return(paste(urlBase,query,sep=""))
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
