correlationTabPanel <- function() {
    fluidRow(column(3,
        wellPanel(
            h4("Gene settings"),
            fluidRow(column(12,
                h4("Gene settings"),
                radioButtons(
                    inputId="rnaCorrelationGeneList",
                    label="Select genes",
                    choices=list(
                        "All genes with non-zero counts"="all",
                        "Custom list"="custom",
                        "Select from list"="select"
                    )
                ),
                conditionalPanel(
                    condition=
                        "input.rnaCorrelationGeneList=='select'",
                    selectizeInput(
                        inputId="selectCorrelationGeneName",
                        label="Select genes to correlate",
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
                        "input.rnaCorrelationGeneList=='custom'",
                    tags$textarea(
                        id="rnaCorrelationCustomList",
                        class="bsv-textarea",
                        rows=5,cols=30,
                        placeholder=paste("Paste gene names ",
                            "separated by newlines",sep="")
                    )
                )
            )),
            fluidRow(column(12,
                h4("Expression measurement"),
                radioButtons(
                    inputId="rnaCorrelationMeasureRadio",
                    label="Select expression measure",
                    choices=list(
                        "Normalized counts"="counts",
                        "RPKM"="rpkm",
                        "RPGM"="rpgm"
                    )
                ),
                radioButtons(
                    inputId="rnaCorrelationScaleRadio",
                    label="Select expression scale",
                    choices=list(
                        "Natural"="natural",
                        "log2"="log2"
                    )
                )
            )),


            fluidRow(column(12,
                h4("Correlation settings"),
                radioButtons(
                    inputId="rnaCorrelateWhat",
                    label="Select type of expression correlation",
                    choices=list(
                        "Sample-wise (all dataset samples)"="samples",
                        "Gene-wise (all selected genes)"="allgenes",
                        "Gene-wise (reference gene against selected)"=
                            "refgene"
                    )
                ),
                conditionalPanel(
                    condition="input.rnaCorrelateWhat=='refgene'",
                    selectizeInput(
                        inputId="rnaCorrelationRefGene",
                        label="Select a reference gene",
                        multiple=FALSE,
                        choices="",
                        options=list(
                            placeholder="Type a gene name",
                            selectOnTab=TRUE
                        )
                    )
                )
            )),
            fluidRow(column(12,
                selectizeInput(
                    inputId="rnaCorrelationMethod",
                    label="Select correlation method",
                    multiple=FALSE,
                    choices=c(
                        Pearson="pearson",
                        Spearman="spearman"
                    )
                )
            )),
            

            fluidRow(br()),
            fluidRow(column(8,
                htmlOutput("rnaCorrelationSettingsError")
            ),column(4,
                 div(
                     class="pull-right",
                     style="display:inline-block",
                     actionButton(
                        inputId="performRnaCorrelation",
                        label="Engage!",
                        icon=icon("rocket")
                    )
                 )
            ))
        )
    ),column(9,
        fluidRow(column(9,
            htmlOutput("correlationOutput")
        ),column(3,
            wellPanel(
                h4("Heatmap controls"),
                fluidRow(column(12,
                    colourInput(
                        inputId="corrColourHigh",
                        label="Colour for high correlation",
                        value="#0000FF"
                    ),
                    colourInput(
                        inputId="corrColourNo",
                        label="Colour for no correlation",
                        value="#BEBEBE"
                    ),
                    colourInput(
                        inputId="corrColourLow",
                        label="Colour for low correlation",
                        value="#FFFF00"
                    )
                )),
                fluidRow(column(12,
                    checkboxInput(
                        inputId="checkRnaSymmCorHeatmap",
                        label="Symmetric colors in heatmap",
                        value=FALSE
                    )
                )),
                fluidRow(br()),
                h4("Export controls"),
                fluidRow(column(12,
                    downloadButton(
                        outputId="exportCorrelationPNG",
                        label="Export PNG",
                        #icon=icon("file-image-o"),
                        class="btn-with-margin"
                    )
                )),
                fluidRow(column(12,
                    downloadButton(
                        outputId="exportCorrelationPDF",
                        label="Export PDF",
                        #icon=icon("file-pdf-o"),
                        class="btn-with-margin"
                    )
                )),
                div(style="display:block;height:160px;","")
            )
        )),
        fluidRow(br()),
        fluidRow(column(5,
            wellPanel(
                h4("Correlation matrix"),
                htmlOutput("rnaCorrelationMatrix")
            )
        ),column(7,
            wellPanel(
                h4("Data matrix"),
                div(
                    class="small table-container",
                    DT::dataTableOutput("rnaCorDataMatrix")
                )
            )
        ))
    ))
}
