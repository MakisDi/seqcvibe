# ui.R

require(DT)
require(shinyjs)

source("config/init.R")
source("ui/dataSelectorTab/dataSelectorItem.R")
source("ui/signalViewerTab/signalViewerItemGeneSignal.R")
source("ui/signalViewerTab/signalViewerItemAreaSignal.R")
source("ui/expressionViewerTab/expressionViewerItemKnownGenes.R")
source("ui/expressionViewerTab/expressionViewerItemCalculator.R")
source("ui/analysisTab/analysisItemDiffExpr.R")
source("ui/analysisTab/analysisItemMdsPca.R")
source("ui/analysisTab/analysisItemCorrelation.R")
source("ui/analysisTab/analysisItemClustering.R")
source("ui/analysisTab/analysisItemGoPathway.R")
source("ui/genomeBrowserTab/genomeBrowserItem.R")
source("ui/helpTab/helpItemDoc.R")
source("ui/helpTab/helpItemFaq.R")

shinyUI(fluidPage(
    shinyjs::useShinyjs(),
    tags$head(
        tags$link(
            rel="stylesheet",
            type="text/css",
            href="seqcvibe.css"
        ),
        tags$link(
            rel="stylesheet",
            type="text/css",
            href="pace.css"
        ),
        tags$script(
            src="pace.js"
        ),
        tags$script(HTML(
            "$(function() {
                setTimeout(function() {
                var vals = [0,1e-16,1e-8,1e-4,0.001,0.005,0.01,0.02,0.03,0.04,
                    0.05,0.1,0.2,0.3,0.4,0.5,1];
                $('#pvalue').data('ionRangeSlider').update({
                    'values': vals,
                    'from': 10
                });
                $('#fdr').data('ionRangeSlider').update({
                    'values': vals,
                    'from': 10
                });
                var fNatVals = [0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,
                    0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,
                    6,7,8,9,10]
                $('#fcNatural').data('ionRangeSlider').update({
                    'values': fNatVals,
                    'from': 10,
                    'to': 20
                });
            },5)})"
        ))
    ),
    #conditionalPanel(
	#	condition="input.universeLoading",
	#	div(
	#		class="fullscreen",div(
	#			class="splash","Starting SeqCVIBE..."
	#		)
	#	)
    #),
    navbarPage(
        title="SeqCVIBE (beta)",
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
            tabPanel("Clustering analysis",
                clusteringTabPanel()
            ),
            tabPanel("Correlation analysis",
                correlationTabPanel()
            ),
            tabPanel("MDS/PCA analysis",
                mdsPcaTabPanel()
            )#,
            #tabPanel("GO/Pathway analysis",
            #    goPathwayTabPanel()
            #)
        ),
        tabPanel("Genome browser",
            genomeBrowserTabPanel()
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
