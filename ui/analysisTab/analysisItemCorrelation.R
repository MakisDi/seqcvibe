correlationTabPanel <- function() {
    fluidRow(column(3,
        wellPanel(
            h4("Correlation settings"),
            tabsetPanel(
                id="correlationSettings",
                tabPanel(
                    title="Genes",
                    fluidRow(br()),
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
                    ))
                ),
                tabPanel(
                    title="Correlation",
                    fluidRow(br()),
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
                    ))
                )            
            ),
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
        fluidRow(column(12,
            htmlOutput("correlationOutput")
        )),
        fluidRow(br()),
        fluidRow(column(8,
            div("")
        ),column(2,
            downloadButton(
                outputId="exportCorrelationPNG",
                label="Export PNG",
                #icon=icon("file-image-o"),
                class="pull-right"
            )
        ),column(2,
            downloadButton(
                outputId="exportCorrelationPDF",
                label="Export PDF",
                #icon=icon("file-pdf-o"),
                class="pull-right"
            )
        ))
    ))
}
