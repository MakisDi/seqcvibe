differentialExpressionTabPanel <- function() {
    fluidRow(column(3,
        fluidRow(column(12,
            wellPanel(
                h4("Pipeline settings"),
                tabsetPanel(
                    id="pipSettingsTabset",
                    tabPanel(
                        title="General",
                        fluidRow(br()),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDePipeline",
                                label="Select pipeline",
                                choices=list(
                                    "DESeq"="deseq",
                                    "edgeR"="edger",
                                    "voom"="limma"
                                )
                            )
                        )),
                        fluidRow(column(12,
                            htmlOutput("rnaDePipelineControl")
                        )),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDeNormalizeWhen",
                                label="Apply filters when",
                                choices=list(
                                    "Before normalization"="prenorm",
                                    "After normalization"="postnorm"
                                ),
                                selected="postnorm"
                            )
                        )),
                        fluidRow(column(12,
                            selectInput(
                                inputId="rnaDeMTC",
                                label="Multiple testing correction",
                                choices=list(
                                    "False Discovery Rate (FDR)"=c(
                                        "Benjamini-Hochberg"="BH",
                                        "Benjamini-Yekutieli"="BY"
                                    ),
                                    "Family Wise Error Rate (FWER)"=c(
                                        "Holm"="holm",
                                        "Hommel"="hommel",
                                        "Bonferroni"="bonferroni"
                                    )
                                )
                            )
                        ))
                    ),
                    tabPanel(
                        title="Filtering",
                        fluidRow(br()),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDeGeneFilter",
                                label="Basic gene filter",
                                choices=list(
                                    "Median expression"="median",
                                    "Mean expression"="mean",
                                    "Quantile"="quantile",
                                    "Known genes"="known"
                                )
                            ),
                            conditionalPanel(
                                condition="input.rnaDeGeneFilter=='quantile'",
                                textInput(
                                    inputId="rnaDeQuantileFilter", 
                                    label="Quantile", 
                                    value="0.25"
                                )
                            ),
                            conditionalPanel(
                                condition="input.rnaDeGeneFilter=='known'",
                                selectizeInput(
                                    inputId="rnaDeKnownFilter",
                                    label="List of genes for filtering",
                                    multiple=TRUE,
                                    choices=NULL,
                                    options=list(
                                        placeholder="Select genes",
                                        selectOnTab=TRUE
                                    )
                                )
                            )
                        )),
                        fluidRow(column(12,
                            div(style="font-weight:bold","Gene length filter"),
                            checkboxInput(
                                inputId="rnaDeGeneLengthFilter",
                                label="Gene length filtering",
                                value=FALSE
                            ),
                            conditionalPanel(
                                condition="input.rnaDeGeneLengthFilter",
                                textInput(
                                    inputId="rnaDeGeneLengthFilterValue",
                                    label="Gene length in bp",
                                    value="500"
                                )
                            )
                        )),
                        fluidRow(column(12,
                            div(style="font-weight:bold","Gene biotype filters"),
                            checkboxInput(
                                inputId="rnaDeBiotypeFilter",
                                label="Biotype expression filtering",
                                value=FALSE
                            ),
                            conditionalPanel(
                                condition="input.rnaDeBiotypeFilter",
                                htmlOutput("checkboxBiotypeListRna")
                            )
                        ))
                    )
                ),
                fluidRow(br()),
                fluidRow(column(8,
                    htmlOutput("rnaDeSettingsError")
                ),column(4,
                     div(
                         class="pull-right",
                         style="display:inline-block",
                         actionButton(
                            inputId="performDeAnalysis",
                            label="Engage!",
                            icon=icon("rocket")
                        )
                     )
                ))
            ),
            wellPanel(
                h4("Results table settings"),
                tabsetPanel(
                    id="rnaDeTableSettings",
                    tabPanel(
                        title="Score",
                        fluidRow(br()),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="statThresholdType",
                                label="Statistical score threshold",
                                inline=TRUE,
                                choices=c(
                                    "p-value"="pvalue",
                                    "FDR"="fdr"
                                )
                            ),
                            conditionalPanel(
                                condition="input.statThresholdType=='pvalue'",
                                sliderInput(
                                    inputId="pvalue",
                                    label="p-value",
                                    min=0,
                                    max=1,
                                    value=0.05,
                                    step=0.01
                                )
                            ),
                            conditionalPanel(
                                condition="input.statThresholdType=='fdr'",
                                sliderInput(
                                    inputId="fdr",
                                    label="FDR",
                                    min=0,
                                    max=1,
                                    value=0.05,
                                    step=0.01
                                )
                            )
                        )),
                        fluidRow(column(12,
                            div(
                                style="font-weight:bold",
                                "Fold change threshold"
                            ),
                            div(
                                class="small",
                                helpText(paste("Select the fold change ",
                                    "display scale (natural or log2) from the ",
                                    "'Scale' tab above."))
                            ),
                            conditionalPanel(
                                condition=
                                    "input.rnaDeValueScaleRadio=='natural'",
                                sliderInput(
                                    inputId="fcNatural",
                                    label="Fold change (natural)",
                                    min=0,
                                    max=10,
                                    value=c(0.5,2),
                                    step=0.5
                                )
                            ),
                            conditionalPanel(
                                condition="input.rnaDeValueScaleRadio=='log2'",
                                sliderInput(
                                    inputId="fcLog",
                                    label="Fold change (log2)",
                                    min=-5,
                                    max=5,
                                    value=c(-1,1),
                                    step=0.5
                                )
                            )
                        ))
                    ),
                    tabPanel(
                        title="Filtering",
                        fluidRow(br()),
                        fluidRow(column(12,
                            htmlOutput("setDeChrs")
                        )),
                        fluidRow(column(12,
                            selectizeInput(
                                inputId="rnaDeShowSpecificGenes",
                                label="Show only selected genes",
                                multiple=TRUE,
                                choices=NULL,
                                options=list(
                                    placeholder="Select genes",
                                    selectOnTab=TRUE
                                )
                            )
                        )),
                        fluidRow(column(12,
                            div(style=
                                "font-weight:bold","Filter by gene biotype"),
                            checkboxInput(
                                inputId="rnaDeAnalyzedBiotypeFilter",
                                label="Biotype expression filtering",
                                value=FALSE
                            ),
                            conditionalPanel(
                                condition="input.rnaDeAnalyzedBiotypeFilter",
                                htmlOutput("checkboxBiotypeListAnalyzedRna")
                            )
                        ))
                    ),
                    tabPanel(
                        title="Scale",
                        fluidRow(br()),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDeValueCompRadio",
                                label="Select summary value components",
                                choices=list(
                                    "Counts"="counts",
                                    "RPKM"="rpkm",
                                    "RPGM"="rpgm"
                                )
                            )
                        )),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDeValueScaleRadio",
                                label="Select summary value scale",
                                choices=list(
                                    "Natural"="natural",
                                    "log2"="log2"
                                )
                            )
                        )),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDeValueAverageRadio",
                                label="Select summary value averaging",
                                choices=list(
                                    "Mean"="mean",
                                    "Median"="median"
                                )
                            )
                        )),
                        fluidRow(column(12,
                            radioButtons(
                                inputId="rnaDeValueDeviationRadio",
                                label="Select summary value deviation",
                                choices=list(
                                    "Standard deviation"="sd",
                                    "Median Absolute Deviation"="mad",
                                    "Interquartile Range"="IQR"
                                )
                            )
                        ))
                    )
                )
            )
        ))
    ),column(9,
        fluidRow(column(12,
            plotOutput("rnaDeMAPlot",height="640px")
        )),
        wellPanel(
            tabsetPanel(
                id="rnaDeAnalysisResults",
                tabPanel(
                    title="Summary",
                    fluidRow(br()),
                    fluidRow(column(12,
                        htmlOutput("rnaDeAnalysisSummary")
                    ))
                ),
                tabPanel(
                    title="Annotation",
                    fluidRow(br()),
                    fluidRow(column(12,
                        htmlOutput("rnaDeAnalysisAnnotation")
                    ))
                ),
                tabPanel(
                    title="Flags",
                    fluidRow(br()),
                    fluidRow(column(12,
                        htmlOutput("rnaDeAnalysisFlags")
                    ))
                ),
                tabPanel(
                    title="All",
                    fluidRow(br()),
                    fluidRow(column(12,
                        htmlOutput("rnaDeAnalysisAll")
                    ))
                )
            )
        )
    ))
}

mdsPcaTabPanel <- function() {
    fluidRow(column(3,
        wellPanel(
            fluidRow(column(12,
                radioButtons(
                    inputId="rnaDimRedMethod",
                    label="Select variance projection method",
                    choices=list(
                        "MultiDimensional Scaling (MDS)"="mds",
                        "Principal Component Analysis (PCA)"="mean"
                    )
                ),
                conditionalPanel(
                    condition="input.rnaDimRedMethod=='pca'",
                    helpText("Some text")
                ),
                conditionalPanel(
                    condition="input.rnaDimRedMethod=='mds'",
                    radioButtons(
                        inputId="rnaMdsCorrMethod",
                        label="Select correlation metric",
                        choices=list(
                            "Pearson"="pearson",
                            "Spearman"="spearman"
                        )
                    )
                )
            ))
        )
    ),column(9,
        wellPanel(
            fluidRow(column(12,
                h4("Coming soon!")
            ))
        )
    ))
}

correlationTabPanel <- function() {
    fluidRow(column(3,
        wellPanel(
            fluidRow(column(12,
                h4("Coming soon!")
            ))
        )
    ),column(9,
        wellPanel(
            fluidRow(column(12,
                h4("Coming soon!")
            ))
        )
    ))
}

pathwayTabPanel <- function() {
    fluidRow(column(3,
        wellPanel(
            fluidRow(column(12,
                h4("Coming soon!")
            ))
        )
    ),column(9,
        wellPanel(
            fluidRow(column(12,
                h4("Coming soon!")
            ))
        )
    ))
}

goTabPanel <- function() {
    fluidRow(column(3,
        wellPanel(
            fluidRow(column(12,
                h4("Coming soon!")
            ))
        )
    ),column(9,
        wellPanel(
            fluidRow(column(12,
                h4("Coming soon!")
            ))
        )
    ))
}
