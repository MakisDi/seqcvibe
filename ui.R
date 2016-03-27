# ui.R

require(DT)
require(shinyjs)

source("config/init.R")
source("client/dataSelectorTab.R")
source("client/signalViewerTab.R")
source("client/expressionViewerTab.R")
source("client/analysisTab.R")
source("client/helpTab.R")

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
			dataSelectorTabPanel()
        ),
        navbarMenu("Signal viewer",
			tabPanel("Gene signal",
				geneSignalTabPanel()
			),
			tabPanel("Area signal",
				areaSignalTabPanel()
			)
		),
        navbarMenu("Expression viewer",
            tabPanel("Known genes",
                expressionExplorerTabPanel()
            ),
            tabPanel("Calculator",
				expressionCalculatorTabPanel()
            )
        ),
        navbarMenu("Analysis",
            tabPanel("Differential expression",
				differentialExpressionTabPanel()
            ),
            tabPanel("MDS/PCA analysis",
				mdsPcaTabPanel()
            ),
            tabPanel("Correlation analysis",
				correlationTabPanel()
            ),
            tabPanel("Pathway analysis",
				pathwayTabPanel()
            ),
            tabPanel("GO analysis",
				goTabPanel()
            )
        ),
        tabPanel("Genome browser",
            fluidRow(column(12,
                htmlOutput("genomeBrowser")
            ))
        ),
        navbarMenu("Help",
            tabPanel("Documentation",
				docTabPanel()
            ),
            tabPanel("FAQ",
				faqTabPanel()
            )
        ),
        tabPanel("About",
            fluidRow(column(12,includeHTML("www/about.html")))
        )
    )
))
