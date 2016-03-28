differentialExpressionTabPanel <- function() {
	fluidRow(column(3,
		fluidRow(column(12,
			wellPanel(
				h4("Pipeline settings"),
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
					radioButtons(
						inputId="rnaDeGeneFilter",
						label="Select gene filter",
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
				)),
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
					radioButtons(
						inputId="foldThresholdType",
						label="Fold threshold",
						inline=TRUE,
						choices=c(
							"Natural scale"="natural",
							"log2 scale"="log2"
						)
					),
					conditionalPanel(
						condition="input.foldThresholdType=='natural'",
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
						condition="input.foldThresholdType=='log2'",
						sliderInput(
							inputId="fcLog",
							label="Fold change (log2)",
							min=-5,
							max=5,
							value=c(-1,1),
							step=0.5
						)
					)
				)),
				fluidRow(column(12,
					htmlOutput("setDeChrs")
				)),
				fluidRow(column(12,
					div(style="font-weight:bold","Filter by gene biotype"),
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
			)
		))
	),column(9,
		fluidRow(column(12,
			plotOutput("rnaDeMAPlot",height="640px")
		))
	))
}

mdsPcaTabPanel <- function() {
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
