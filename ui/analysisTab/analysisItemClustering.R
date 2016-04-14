clusteringTabPanel <- function() {
    fluidRow(column(3,
        wellPanel(
            fluidRow(column(12,
                h4("Gene settings"),
                fluidRow(column(12,
                    radioButtons(
                        inputId="rnaClusterGeneList",
                        label="Select genes for clustering",
                        choices=list(
                            "Selected genes"="select",
                            "Differentially expressed genes"="degenes"
                        )
                    ),
                    conditionalPanel(
                        condition=
                            "input.rnaClusterGeneList=='select'",
                        selectizeInput(
                            inputId="selectClusteringGeneName",
                            label="",
                            multiple=TRUE,
                            choices=NULL,
                            options=list(
                                placeholder="Select genes to cluster",
                                selectOnTab=TRUE
                            )
                        )
                    ),
                    conditionalPanel(
                        condition=
                            "input.rnaClusterGeneList=='degenes'",
                        div(class="small",
                            helpText(paste("Genes currently present in the ",
                                "differentially expressed genes table in the ",
                                "'Differential expression' panel will be ",
                                "clustered. Please be careful with the number ",
                                "of genes to be clustered (e.g. p-value=1). ",
                                "Very long computations will be interrupted ",
                                "for application safety.",sep=""))
                        )
                    )
                ))
            ))
        ),
        wellPanel(
            fluidRow(column(12,
                h4("Clustering settings"),
                selectizeInput(
                    inputId="selectClusteringDistance",
                    label="Select distance metric",
                    choices=c(
                        "Euclidean"="euclidean",
                        "Maximum"="maximum",
                        "Manhattan"="manhattan",
                        "Canberra"="canberra",
                        "Minkowski"="minkowski",
                        "Pearson correlation"="pearson",
                        "Spearman correlation"="spearman",
                        "Cosine"="cosine"
                    )
                ),
                selectizeInput(
                    inputId="selectClusteringLinkage",
                    label="Select linkage function (dendrogram)",
                    choices=c(
                        "Average"="average",
                        "Complete"="complete",
                        "Single"="single",
                        "McQuitty"="mcquitty",
                        "Median"="median",
                        "Centroid"="centroid",
                        "Ward v1"="ward1",
                        "Ward v2"="ward2"
                    )
                ),
                radioButtons(
                    inputId="rnaClusterWhat",
                    label="Select clustering dimensions",
                    choices=list(
                        "Cluster both rows and columns"="both",
                        "Cluster rows"="row",
                        "Cluster columns"="column",
                        "No clustering (only show image)"="none"
                    )
                )
            )),
            fluidRow(column(6,
                textInput(
                    inputId="kRowInput",
                    label="# of row clusters",
                    value="1"
                )
            ),column(6,
                textInput(
                    inputId="kColInput",
                    label="# of column clusters",
                    value="1"
                )
            )),
            fluidRow(br()),
            fluidRow(column(8,
                htmlOutput("rnaClusteringSettingsError")
            ),column(4,
                 div(
                     class="pull-right",
                     style="display:inline-block",
                     actionButton(
                        inputId="performRnaClustering",
                        label="Engage!",
                        icon=icon("rocket")
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

