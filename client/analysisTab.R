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
			)
		))
	),column(9,
		fluidRow(column(12,
			wellPanel(
				htmlOutput("rnaDeMAPlot")
			)
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
