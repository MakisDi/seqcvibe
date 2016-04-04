# ui.R

require(DT)
require(shinyjs)

source("config/init.R")
source("client/dataSelectorTab.R")
source("client/signalViewerTab.R")
source("client/expressionViewerTab.R")
source("client/analysisTab.R")
source("client/genomeBrowserTab.R")
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
            },5)})"
        ))
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
