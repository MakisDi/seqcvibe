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
