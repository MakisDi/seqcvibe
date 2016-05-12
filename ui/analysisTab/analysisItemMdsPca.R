mdsPcaTabPanel <- function() {
    fluidRow(column(3,
        wellPanel(
            h4("Setings"),
            tabsetPanel(
                id="mdsPcaSettings",
                tabPanel(
                    title="Genes",
                    fluidRow(column(12,
                        h4("Gene settings"),
                        radioButtons(
                            inputId="rnaMdsPcaGeneList",
                            label="Select genes",
                            choices=list(
                                "All genes with non-zero counts"="all",
                                "All expressed* genes"="expr",
                                "Custom list"="custom",
                                "Select from list"="select"
                            )
                        ),
                        conditionalPanel(
                            condition=
                                "input.rnaMdsPcaGeneList=='select'",
                            selectizeInput(
                                inputId="selectMdsPcaGeneName",
                                label="Select genes for analysis" ,
                                multiple=TRUE,
                                choices=NULL,
                                options=list(
                                    placeholder="Type some gene names",
                                    selectOnTab=TRUE
                                )
                            )
                        ),
                        conditionalPanel(
                            condition=
                                "input.rnaMdsPcaGeneList=='custom'",
                            tags$textarea(
                                id="rnaMdsPcaCustomList",
                                class="bsv-textarea",
                                rows=5,cols=30,
                                placeholder=paste("Paste gene names ",
                                    "separated by newlines",sep="")
                            )
                        ),
                        div(
                            class="small",
                            helpText(paste("*Passing the filters in the ",
                                "differential expression analysis tab (an ",
                                "analysis must have run first).",sep=""))
                        )
                    )),
                    fluidRow(column(12,
                        h4("Expression measurement"),
                        radioButtons(
                            inputId="rnaMdsPcaMeasureRadio",
                            label="Select expression measure",
                            choices=list(
                                "Normalized counts"="counts",
                                "RPKM"="rpkm",
                                "RPGM"="rpgm"
                            )
                        ),
                        radioButtons(
                            inputId="rnaMdsPcaScaleRadio",
                            label="Select expression scale",
                            choices=list(
                                "Natural"="natural",
                                "log2"="log2"
                            )
                        )
                    )),
                    fluidRow(column(12,
                        div(
                            style="font-weight:bold;margin-top:10px;",
                            "Custom regions"
                        ),
                        checkboxInput(
                            inputId="includeCustomRegionsInDimRed",
                            label=paste("Include regions from expression ",
                                "calculator"),
                            value=FALSE
                        )
                    ))
                ),
                tabPanel(
                    title="MDS/PCA",
                    fluidRow(column(12,
                        h4("MDS/PCA settings"),
                        radioButtons(
                            inputId="rnaDimRedMethod",
                            label="Select variance projection method",
                            choices=list(
                                "MultiDimensional Scaling (MDS)"="mds",
                                "Principal Component Analysis (PCA)"="pca"
                            )
                        ),
                        conditionalPanel(
                            condition="input.rnaDimRedMethod=='mds'",
                            selectizeInput(
                                inputId="rnaMdsDistMethod",
                                label="Select (dis)similarity metric",
                                choices=c(
                                    "Euclidean"="euclidean",
                                    "Maximum"="maximum",
                                    "Manhattan"="manhattan",
                                    #"Canberra"="canberra",
                                    "Minkowski"="minkowski",
                                    "Pearson correlation"="pearson",
                                    "Spearman correlation"="spearman",
                                    "Cosine"="cosine"
                                )
                            ),
                            selectizeInput(
                                inputId="rnaMdsKDim",
                                label="Select number of dimensions",
                                choices=NULL
                            )
                        ),
                        conditionalPanel(
                            condition="input.rnaDimRedMethod=='pca'",
                            div(style="font-weight:bold","PCA settings"),
                            fluidRow(column(6,
                                checkboxInput(
                                    inputId="rnaPcaDoCenter",
                                    label="Center values",
                                    value=TRUE
                                )
                            ),column(6,
                                checkboxInput(
                                    inputId="rnaPcaDoScale",
                                    label="Scale values",
                                    value=FALSE
                                )
                            ))
                        )
                    ))
                )
            ),
            fluidRow(br()),
            fluidRow(column(8,
                htmlOutput("rnaMdsPcaSettingsError")
            ),column(4,
                 div(
                     class="pull-right",
                     style="display:inline-block",
                     actionButton(
                        inputId="performRnaMdsPca",
                        label="Engage!",
                        icon=icon("rocket")
                    )
                 )
            ))
        ),
        wellPanel(
            fluidRow(column(12,
                h4("Coordinate plotting"),
                selectizeInput(
                    inputId="rnaDimRedXAxis",
                    label="Horizontal (x) axis projection",
                    choices=NULL
                ),
                selectizeInput(
                    inputId="rnaDimRedYAxis",
                    label="Verical (y) axis projection",
                    choices=NULL
                ),
                div(
                    class="small",
                    helpText(paste("When the variance projection method is ",
                        "MDS, PC refers to 'Principal Coordinate'. When the ",
                        "former is PCA, the latter refers to 'Principal ",
                        "Component'.",sep=""))
                ),
                div(
                    style="font-weight:bold;margin-top:10px;",
                    "Custom regions"
                ),
                checkboxInput(
                    inputId="rnaMdsPcaTogglePointNames",
                    label="Toggle sample names",
                    value=FALSE
                )
            )),
            fluidRow(column(12,
                htmlOutput("rnaMdsPcaPlotColours")
            ))
        )
    ),column(9,
        fluidRow(column(12,
            conditionalPanel(
                condition="input.rnaDimRedMethod=='mds'",
                tabsetPanel(
                    id="rnaDimRedMdsPlots",
                    tabPanel(
                        title="Eigenvectors",
                        plotOutput(
                            outputId="rnaMdsPlot",
                            click="rnaMdsPlotClick",
                            #dblclick="rnaMdsPlotDblClick",
                            brush=brushOpts(
                                id="rnaMdsPlotBrush",
                                resetOnNew=TRUE
                            ),
                            height="640px"
                        )
                    )
                )
            ),
            conditionalPanel(
                condition="input.rnaDimRedMethod=='pca'",
                tabsetPanel(
                    id="rnaDimRedPcaPlots",
                    tabPanel(
                        title="Scree plot",
                        plotOutput("rnaPcaScreePlot",height="640px")
                    ),
                    tabPanel(
                        title="Scores",
                        plotOutput(
                            outputId="rnaPcaScoresPlot",
                            click="rnaPcaScoresPlotClick",
                            #dblclick="rnaPcaScoresPlotDblClick",
                            brush=brushOpts(
                                id="rnaPcaScoresPlotBrush",
                                resetOnNew=TRUE
                            ),
                            height="640px"
                        )
                    ),
                    tabPanel(
                        title="Loadings 2D",
                        plotOutput(
                            outputId="rnaPcaLoadingsPlot",
                            click="rnaPcaLoadingsPlotClick",
                            #dblclick="rnaPcaLoadingsPlotDblClick",
                            brush=brushOpts(
                                id="rnaPcaLoadingsPlotBrush",
                                resetOnNew=TRUE
                            ),
                            height="640px"
                        )
                    ),
                    tabPanel(
                        title="Loadings 1D",
                        plotOutput(
                            outputId="rnaPcaRankedLoadingsPlot",
                            click="rnaPcaRankedLoadingsPlotClick",
                            #dblclick="rnaPcaScoresPlotDblClick",
                            brush=brushOpts(
                                id="rnaPcaRankedLoadingsPlotBrush",
                                resetOnNew=TRUE
                            ),
                            height="640px"
                        )
                    ),
                    tabPanel(
                        title="Biplot",
                        plotOutput("rnaPcaBiplotPlot",height="640px")
                    )
                )
            )
        )),
        fluidRow(br()),
        fluidRow(column(12,
            conditionalPanel(
                condition="input.rnaDimRedMethod=='mds'",
                verbatimTextOutput("rnaMdsPcaDisplay")
            )
        )),
        wellPanel(
            fluidRow(column(12,
                htmlOutput("rnaMdsPcaGenesSamples")
            ))
        )
    ))
}
