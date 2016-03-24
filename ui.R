# ui.R

require(DT)
require(shinyjs)

source("config/init.R")

shinyUI(fluidPage(
    shinyjs::useShinyjs(),
    tags$head(
        tags$link(
            rel="stylesheet",
            type="text/css",
            href="bigseqcvis.css"
        ),
        tags$link(
            rel="stylesheet",
            type="text/css",
            href="pace.css"
        ),
        tags$script(
            src="pace.js"
        )
    ),
    navbarPage(
        title="BigSeqCVis (beta)",
        tabPanel("Data selector",
            fluidRow(column(6,
                fluidRow(column(12,
                    wellPanel(
                        h4("Select input data"),
                        fluidRow(column(5,
                            htmlOutput("dataSource")
                        ),column(5,
                            htmlOutput("dataDataset")
                        ),column(2,
                            htmlOutput("dataGenome")
                        )),
                        fluidRow(column(12,
                            htmlOutput("dataSelectHint")
                        )),
                        fluidRow(column(12,
                            htmlOutput("dataCustomSamples")
                        )),
                        fluidRow(column(6,
                            div(
                                class="pull-left",
                                actionButton(
                                    inputId="loadDataset",
                                    label="Load selected data",
                                    icon=icon("truck")
                                )
                            )
                        ),column(3,
                            div(
                                class="pull-right",
                                actionButton(
                                    inputId="clearDataset",
                                    label="Clear dataset",
                                    icon=icon("trash")
                                )
                            )
                        ),column(3,
                            div(
                                class="pull-right",
                                disabled(actionButton(
                                    inputId="createDataset",
                                    label="Create dataset",
                                    icon=icon("table"),
                                    class="btn-sample-select"
                                ))
                            )
                        ))
                    )
                ))
            ),column(6,
                fluidRow(column(12,
                    wellPanel(
                        h4("Current dataset"),
                        fluidRow(column(12,
                            div(
                                class="small table-contaner",
                                DT::dataTableOutput("currentDatasetTable")
                            )
                        ))
                    )
                )),
                fluidRow(column(12,
                    wellPanel(
                        h4("Mesages"),
                        fluidRow(column(12,
                            div(
                                class="message-box",
                                htmlOutput("dataSelectorMessages")
                            )
                        ))
                    )
                ))
            ))
        ),
        tabPanel("Gene explorer",
            fluidRow(column(3,
                fluidRow(column(12,
                    wellPanel(
                        h4("Plot data"),
                        tabsetPanel(
                            id="geneType",
                            tabPanel(
                                title="Genes",
                                fluidRow(column(12,
                                    selectizeInput(
                                        inputId="geneGeneName",
                                        label="",
                                        multiple=TRUE,
                                        choices=NULL,
                                        options=list(
                                            placeholder="Select genes",
                                            selectOnTab=TRUE#,
                                            #create=FALSE#,
                                            #valueField="symbol",
                                            #labelField="symbol",
                                            #searchField="symbol",
                                            #load=I("function(query,callback) {
                                            #   if (!query.length) return callback();
                                            #   $.ajax({
                                            #       url: 'http://rest.genenames.org/search/symbol/'+encodeURIComponent(query+'*'),
                                            #       type: 'GET',
                                            #       dataType: 'json',
                                            #       error: function() {
                                            #           callback();
                                            #       },
                                            #       success: function(res) {
                                            #           callback(res.response.docs);
                                            #       }
                                            #   });
                                            #}")
                                        )
                                    )
                                )),
                                fluidRow(column(12,
                                    div(class="small table-container",
                                        DT::dataTableOutput("knownGeneList")
                                    )
                                ))
                            ),
                            tabPanel(
                                title="Custom",
                                fluidRow(column(8,
                                    textInput(
                                        inputId="customName", 
                                        label="", 
                                        value="",
                                        placeholder="Name"
                                    )
                                ),column(4,
                                    selectizeInput(
                                        inputId="customStrand",
                                        label="",
                                        choices=c("","+","-"),
                                        options=list(placeholder="Strand")
                                    )
                                )),
                                fluidRow(column(4,
                                    htmlOutput("setChrs")
                                ),column(4,
                                    textInput(
                                        inputId="customStart", 
                                        label="", 
                                        value="",
                                        placeholder="Start"
                                    )
                                ),column(4,
                                    textInput(
                                        inputId="customEnd",
                                        label="",
                                        value="",
                                        placeholder="End"
                                    )
                                )),
                                fluidRow(column(12,
                                    div(class="small table-container",
                                        DT::dataTableOutput("customRegionList")
                                    )
                                )),
                                fluidRow(column(12,
                                    htmlOutput("customRegionError")
                                )),
                                fluidRow(br()),
                                fluidRow(column(4,""
                                ),column(4,
                                    div(
                                        class="pull-right",
                                        style="display:inline-block",
                                        actionButton("removeCustomRegion",
                                            "Remove",icon=icon("minus"))
                                    )
                                ),column(4,
                                    div(
                                        class="pull-right",
                                        style="display:inline-block",
                                        actionButton("addCustomRegion",
                                            "Add",icon=icon("plus"))
                                    )
                                ))
                            )
                        )
                    )
                )),
                fluidRow(column(12,
                    wellPanel(
                        h4("Plot options"),
                        h5("Flanking",style="font-size:1.2em;margin-top:20px;"),
                        fluidRow(column(6,
                            textInput(
                                inputId="upstreamFlank", 
                                label="Upstream", 
                                value=2000
                            )
                        ),
                        column(6,
                            textInput(
                                inputId="downstreamFlank",
                                label="Downstream",
                                value=2000
                            )
                        )),
                        fluidRow(column(12,
                            htmlOutput("regionFlankError")
                        )),
                        fluidRow(br()),
                        fluidRow(column(12,
                            htmlOutput("geneExplorerColours")
                        )),
                        fluidRow(br()),
                        fluidRow(column(8,
                            radioButtons(
                                inputId="geneSumStatType",
                                label="Gene profile averaging",
                                choices=list(
                                    "Mean"="mean",
                                    "Median"="median",
                                    "Trimmed mean"="trimmed"
                                )
                            )
                        ),column(4,
                            textInput(
                                inputId="geneTrimPct", 
                                label="Trim fraction", 
                                value="0.1"
                            )
                        )),
                        fluidRow(br()),
                        fluidRow(column(8,
                            htmlOutput("geneExplorerError")
                        ),column(4,
                             div(
                                 class="pull-right",
                                 style="display:inline-block",
                                 actionButton(
                                    inputId="createGeneProfile",
                                    label="Engage!",
                                    icon=icon("rocket")
                                )
                             )
                        ))
                    )
                ))
            ),column(9,
                fluidRow(column(12,
                    plotOutput("geneProfile",height="640px")
                )),
                fluidRow(column(8,
                    div("")
                ),
                column(2,
                    downloadButton(
                        outputId="exportGenePNG",
                        label="Export PNG",
                        #icon=icon("file-image-o"),
                        class="pull-right"
                    )
                ),column(2,
                    downloadButton(
                        outputId="exportGenePDF",
                        label="Export PDF",
                        #icon=icon("file-pdf-o"),
                        class="pull-right"
                    )
                )),
                fluidRow(br()),
                fluidRow(column(12,
                    wellPanel(
                        h4("Messages"),
                        div(
                            class="message-box",
                            htmlOutput("geneExplorerMessages")
                        )
                    )
                ))
            ))
        ),
        tabPanel("Area explorer",
            fluidRow(column(3,
                fluidRow(column(12,
                    wellPanel(
                        h4("Plot data"),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="areaTypeRadio",
                                label="Select type of genomic area",
                                #inline=TRUE,
                                choices=list(
                                    "Custom genomic area"="area",
                                    "Around gene"="gene"
                                )
                            )
                        )),
                        conditionalPanel(
                            condition="input.areaTypeRadio=='area'",
                            fluidRow(column(4,
                                htmlOutput("setChrsA")
                            ),column(4,
                                textInput(
                                    inputId="customStartA", 
                                    label="", 
                                    value="",
                                    placeholder="Start"
                                )
                            ),column(4,
                                textInput(
                                    inputId="customEndA",
                                    label="",
                                    value="",
                                    placeholder="End"
                                )
                            )),
                            div(class="small",helpText(
                                paste("Known genes overlaping the area above ",
                                    "will be automatically detected. Check ",
                                    "the following box if you want to add ",
                                    "custom areas to the automatic detection.")
                            )),
                            checkboxInput(
                                inputId="addCustomRegionsInAreaCheck", 
                                label="Include custom transcribed regions", 
                                value=FALSE
                            ),
                            conditionalPanel(
                                condition="input.addCustomRegionsInAreaCheck",
                                radioButtons(
                                    inputId="customRegionsFromGeneExplorer",
                                    label="Custom regions to include",
                                    choices=list(
                                        "From gene explorer (if any, same chromosome only)"="fromge",
                                        "Other custom regions"="other"
                                    )
                                ),
                                conditionalPanel(
                                    condition="input.customRegionsFromGeneExplorer=='other'",
                                    helpText("Please describe custom regions"),
                                    fluidRow(column(8,
                                        textInput(
                                            inputId="customNameInArea",
                                            label="", 
                                            value="",
                                            placeholder="Name"
                                        )
                                    ),column(4,
                                        selectizeInput(
                                            inputId="customStrandInArea",
                                            label="",
                                            choices=c("","+","-"),
                                            options=list(placeholder="Strand")
                                        )
                                    )),
                                    fluidRow(column(4,
                                        disabled(textInput(
                                            inputId="customChrInArea",
                                            label="", 
                                            value="",
                                            placeholder="Chrom"
                                        ))
                                    ),column(4,
                                        textInput(
                                            inputId="customStartInArea", 
                                            label="", 
                                            value="",
                                            placeholder="Start"
                                        )
                                    ),column(4,
                                        textInput(
                                            inputId="customEndInArea",
                                            label="",
                                            value="",
                                            placeholder="End"
                                        )
                                    )),
                                    fluidRow(column(12,
                                        div(class="small table-container",
                                            DT::dataTableOutput(
                                                "customRegionListInArea"
                                            )
                                        )
                                    )),
                                    fluidRow(column(12,
                                        htmlOutput("customRegionErrorInArea")
                                    )),
                                    fluidRow(br()),
                                    fluidRow(column(4,""
                                    ),column(4,
                                        div(
                                            class="pull-right",
                                            style="display:inline-block",
                                            actionButton("removeCustomRegionInArea",
                                                "Remove",icon=icon("minus"))
                                        )
                                    ),column(4,
                                        div(
                                            class="pull-right",
                                            style="display:inline-block",
                                            actionButton("addCustomRegionInArea",
                                                "Add",icon=icon("plus"))
                                        )
                                    ))
                                )
                            )
                        ),
                        conditionalPanel(
                            condition="input.areaTypeRadio=='gene'",
                            fluidRow(column(12,
                                selectizeInput(
                                    inputId="areaGeneName",
                                    label="",
                                    multiple=FALSE,
                                    choices="",
                                    options=list(
                                        placeholder="Select a gene",
                                        selectOnTab=TRUE
                                    )
                                )
                            )),
                            h5("Flanking",style="font-size:1.2em;"),
                            fluidRow(column(6,
                                textInput(
                                    inputId="upstreamFlankA", 
                                    label="Upstream", 
                                    value=2000
                                )
                            ),column(6,
                                textInput(
                                    inputId="downstreamFlankA",
                                    label="Downstream",
                                    value=2000
                                )
                            ))
                        ),
                        fluidRow(column(12,
                            htmlOutput("areaFlankError")
                        ))
                    )
                )),
                fluidRow(column(12,
                    wellPanel(
                        h4("Plot options"),
                        fluidRow(column(12,
                            htmlOutput("areaExplorerColours")
                        )),
                        fluidRow(br()),
                        fluidRow(column(8,
                            radioButtons(
                                inputId="areaSumStatType",
                                label="Area profile averaging",
                                choices=list(
                                    "Mean"="mean",
                                    "Median"="median",
                                    "Trimmed mean"="trimmed"
                                )
                            )
                        ),column(4,
                            textInput(
                                inputId="areaTrimPct", 
                                label="Trim fraction", 
                                value="0.1"
                            )
                        )),
                        fluidRow(column(8,
                            htmlOutput("areaExplorerError")
                        ),column(4,
                            div(
                                class="pull-right",
                                style="display:inline-block",
                                actionButton(
                                    inputId="createAreaProfile",
                                    label="Engage!",
                                    icon=icon("rocket")
                                )
                            )
                        ))
                    )
                ))
            ),column(9,
                fluidRow(column(12,
                    plotOutput("areaProfile",height="640px")
                )),
                fluidRow(column(8,
                    div("")
                ),
                column(2,
                    downloadButton(
                        outputId="exportAreaPNG",
                        label="Export PNG",
                        #icon=icon("file-image-o"),
                        class="pull-right"
                    )
                ),column(2,
                    downloadButton(
                        outputId="exportAreaPDF",
                        label="Export PDF",
                        #icon=icon("file-pdf-o"),
                        class="pull-right"
                    )
                )),
                fluidRow(br()),
                fluidRow(column(12,
                    wellPanel(
                        h4("Messages"),
                        div(
                            class="message-box",
                            htmlOutput("areaExplorerMessages")
                        )
                    )
                ))
            ))
        ),
        navbarMenu("Expression",
            tabPanel("Explorer",
                fluidRow(column(3,
                    fluidRow(column(12,
                        wellPanel(
                            h4("Gene settings"),
                            fluidRow(column(12,
                                radioButtons(
                                    inputId="rnaExpressionGeneList",
                                    label="Select genes",
                                    choices=list(
                                        "Select from list"="select",
                                        "Custom list"="custom",
                                        "All genes"="all"
                                    )
                                ),
                                conditionalPanel(
                                    condition=
                                        "input.rnaExpressionGeneList=='select'",
                                    selectizeInput(
                                        inputId="selectExpressionGeneName",
                                        label="",
                                        multiple=TRUE,
                                        choices=NULL,
                                        options=list(
                                            placeholder="Select genes",
                                            selectOnTab=TRUE
                                        )
                                    )
                                ),
                                conditionalPanel(
                                    condition=
                                        "input.rnaExpressionGeneList=='custom'",
                                    tags$textarea(
                                        id="rnaExpressionCustomList",
                                        class="bsv-textarea",
                                        rows=5,cols=40,
                                        placeholder=paste("Paste gene names ",
                                            "separated by newlines",sep="")
                                    )
                                )
                            ))
                        )
                    )),
                    fluidRow(column(12,
                        wellPanel(
                            h4("Expression settings"),
                            fluidRow(column(12,
                                radioButtons(
                                    inputId="rnaExpressionMeasureRadio",
                                    label="Select RNA-Seq expression measure",
                                    choices=list(
                                        "Raw counts"="raw",
                                        "DESeq Normalized counts"="norm",
                                        "RPKM"="rpkm",
                                        "RPGM"="rpgm"
                                    )
                                )
                            )),
                            fluidRow(column(12,
                                radioButtons(
                                    inputId="rnaExpressionScaleRadio",
                                    label="Select expression measure scale",
                                    choices=list(
                                        "Natural"="natural",
                                        "log2"="log2"
                                    )
                                )
                            )),
                            fluidRow(column(12,
                                radioButtons(
                                    inputId="rnaExpressionAverageRadio",
                                    label="Select expression measure averaging",
                                    choices=list(
                                        "Mean"="mean",
                                        "Median"="median"
                                    )
                                )
                            )),
                            fluidRow(column(12,
                                radioButtons(
                                    inputId="rnaExpressionDeviationRadio",
                                    label="Select expression deviation measure",
                                    choices=list(
                                        "Standard deviation"="sd",
                                        "Median Absolute Deviation"="mad",
                                        "Interquartile Range"="IQR"
                                    )
                                )
                            ))
                        )
                    ))
                ),column(9,
                    fluidRow(column(12,
                        wellPanel(
                            htmlOutput("rnaExpressionTables")
                        )
                    ))
                ))
            ),
            tabPanel("Calculator",
                fluidRow(column(3,
                    fluidRow(column(12,
                        wellPanel(
                            h4("Region settings"),
                            radioButtons(
                                inputId="customExpressionFromGeneExplorer",
                                label="Custom region settings",
                                choices=list(
                                    "From gene explorer"="fromge",
                                    "Custom regions"="custom"
                                )
                            ),
                            conditionalPanel(condition=
                                "input.customExpressionFromGeneExplorer=='custom'",
                                fluidRow(column(8,
                                    textInput(
                                        inputId="customNameExpr", 
                                        label="", 
                                        value="",
                                        placeholder="Name"
                                    )
                                ),column(4,
                                    selectizeInput(
                                        inputId="customStrandExpr",
                                        label="",
                                        choices=c("","+","-","*"),
                                        options=list(placeholder="Strand")
                                    )
                                )),
                                fluidRow(column(4,
                                    htmlOutput("setChrsExpr")
                                ),column(4,
                                    textInput(
                                        inputId="customStartExpr", 
                                        label="", 
                                        value="",
                                        placeholder="Start"
                                    )
                                ),column(4,
                                    textInput(
                                        inputId="customEndExpr",
                                        label="",
                                        value="",
                                        placeholder="End"
                                    )
                                )),
                                fluidRow(column(12,
                                    div(class="small",
                                        DT::dataTableOutput(
                                            "customRegionListExpr")
                                    )
                                )),
                                fluidRow(column(12,
                                    htmlOutput("customRegionExprError")
                                )),
                                fluidRow(br()),
                                fluidRow(column(4,""
                                ),column(4,
                                    div(
                                        class="pull-right",
                                        style="display:inline-block",
                                        actionButton("removeRnaCustomRegion",
                                            "Remove",icon=icon("minus"))
                                    )
                                ),column(4,
                                    div(
                                        class="pull-right",
                                        style="display:inline-block",
                                        actionButton("addRnaCustomRegion",
                                            "Add",icon=icon("plus"))
                                    )
                                ))
                            ),
							fluidRow(br()),
							fluidRow(column(8,
								htmlOutput("customRnaCalcError")
							),column(4,
								div(
									class="pull-right",
									style="display:inline-block",
									actionButton(
										inputId="calculateCustomRegionRna",
										label="Engage!",
										icon=icon("rocket")
									)
								)
							))
                        )
                    )),
                    fluidRow(column(12,
                        wellPanel(
                            h4("Expression settings"),
                            fluidRow(column(12,
                                radioButtons(
                                    inputId="rnaCustomMeasureRadio",
                                    label="Select RNA-Seq expression measure",
                                    choices=list(
                                        "Raw counts"="raw",
                                        #"DESeq Normalized counts"="norm",
                                        "RPKM"="rpkm",
                                        "RPGM"="rpgm"
                                    )
                                )
                            )),
                            fluidRow(column(12,
                                radioButtons(
                                    inputId="rnaCustomScaleRadio",
                                    label="Select expression measure scale",
                                    choices=list(
                                        "Natural"="natural",
                                        "log2"="log2"
                                    )
                                )
                            )),
                            fluidRow(column(12,
                                radioButtons(
                                    inputId="rnaCustomAverageRadio",
                                    label="Select expression measure averaging",
                                    choices=list(
                                        "Mean"="mean",
                                        "Median"="median"
                                    )
                                )
                            )),
                            fluidRow(column(12,
                                radioButtons(
                                    inputId="rnaCustomDeviationRadio",
                                    label="Select expression deviation measure",
                                    choices=list(
                                        "Standard deviation"="sd",
                                        "Median Absolute Deviation"="mad",
                                        "Interquartile Range"="IQR"
                                    )
                                )
                            ))
                        )
                    ))
                ),column(9,
                    fluidRow(column(12,
                        wellPanel(
                            htmlOutput("rnaCustomTables")
                        )
                    ))
                ))
            ),
            tabPanel("Differential",
                fluidRow(column(3,
                    fluidRow(column(12,
                        div("Coming soon!")
                    ))
                ),column(9,
                    div("Coming soon!")
                ))
            )
        ),
        tabPanel("Genome browser",
            fluidRow(column(12,
                htmlOutput("genomeBrowser")
            ))
        ),
        tabPanel("About",
            fluidRow(column(12,includeHTML("www/about.html")))
        )
    )
))

