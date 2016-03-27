geneSignalTabPanel <- function() {
	fluidRow(column(3,
		fluidRow(column(12,
			wellPanel(
				h4("Plot data"),
				tabsetPanel(
					id="geneType",
					tabPanel(
						title="Genes",
						fluidRow(column(12,
							selectizeInput(
								inputId="geneGeneName",
								label="",
								multiple=TRUE,
								choices=NULL,
								options=list(
									placeholder="Select genes",
									selectOnTab=TRUE#,
									#create=FALSE#,
									#valueField="symbol",
									#labelField="symbol",
									#searchField="symbol",
									#load=I("function(query,callback) {
									#   if (!query.length) return callback();
									#   $.ajax({
									#       url: 'http://rest.genenames.org/search/symbol/'+encodeURIComponent(query+'*'),
									#       type: 'GET',
									#       dataType: 'json',
									#       error: function() {
									#           callback();
									#       },
									#       success: function(res) {
									#           callback(res.response.docs);
									#       }
									#   });
									#}")
								)
							)
						)),
						fluidRow(column(12,
							div(class="small table-container",
								DT::dataTableOutput("knownGeneList")
							)
						))
					),
					tabPanel(
						title="Custom",
						fluidRow(column(8,
							textInput(
								inputId="customName", 
								label="", 
								value="",
								placeholder="Name"
							)
						),column(4,
							selectizeInput(
								inputId="customStrand",
								label="",
								choices=c("","+","-"),
								options=list(placeholder="Strand")
							)
						)),
						fluidRow(column(4,
							htmlOutput("setChrs")
						),column(4,
							textInput(
								inputId="customStart", 
								label="", 
								value="",
								placeholder="Start"
							)
						),column(4,
							textInput(
								inputId="customEnd",
								label="",
								value="",
								placeholder="End"
							)
						)),
						fluidRow(column(12,
							div(class="small table-container",
								DT::dataTableOutput("customRegionList")
							)
						)),
						fluidRow(column(12,
							htmlOutput("customRegionError")
						)),
						fluidRow(br()),
						fluidRow(column(4,""
						),column(4,
							div(
								class="pull-right",
								style="display:inline-block",
								actionButton("removeCustomRegion",
									"Remove",icon=icon("minus"))
							)
						),column(4,
							div(
								class="pull-right",
								style="display:inline-block",
								actionButton("addCustomRegion",
									"Add",icon=icon("plus"))
							)
						))
					)
				)
			)
		)),
		fluidRow(column(12,
			wellPanel(
				h4("Plot options"),
				h5("Flanking",style="font-size:1.2em;margin-top:20px;"),
				fluidRow(column(6,
					textInput(
						inputId="upstreamFlank", 
						label="Upstream", 
						value=2000
					)
				),
				column(6,
					textInput(
						inputId="downstreamFlank",
						label="Downstream",
						value=2000
					)
				)),
				fluidRow(column(12,
					htmlOutput("regionFlankError")
				)),
				fluidRow(br()),
				fluidRow(column(12,
					htmlOutput("geneExplorerColours")
				)),
				fluidRow(br()),
				fluidRow(column(8,
					radioButtons(
						inputId="geneSumStatType",
						label="Gene profile averaging",
						choices=list(
							"Mean"="mean",
							"Median"="median",
							"Trimmed mean"="trimmed"
						)
					)
				),column(4,
					textInput(
						inputId="geneTrimPct", 
						label="Trim fraction", 
						value="0.1"
					)
				)),
				fluidRow(br()),
				fluidRow(column(8,
					htmlOutput("geneExplorerError")
				),column(4,
					 div(
						 class="pull-right",
						 style="display:inline-block",
						 actionButton(
							inputId="createGeneProfile",
							label="Engage!",
							icon=icon("rocket")
						)
					 )
				))
			)
		))
	),column(9,
		fluidRow(column(12,
			plotOutput("geneProfile",height="640px")
		)),
		fluidRow(column(8,
			div("")
		),
		column(2,
			downloadButton(
				outputId="exportGenePNG",
				label="Export PNG",
				#icon=icon("file-image-o"),
				class="pull-right"
			)
		),column(2,
			downloadButton(
				outputId="exportGenePDF",
				label="Export PDF",
				#icon=icon("file-pdf-o"),
				class="pull-right"
			)
		)),
		fluidRow(br()),
		fluidRow(column(12,
			wellPanel(
				h4("Messages"),
				div(
					class="message-box",
					htmlOutput("geneExplorerMessages")
				)
			)
		))
	))
}

areaSignalTabPanel <- function() {
	fluidRow(column(3,
		fluidRow(column(12,
			wellPanel(
				h4("Plot data"),
				fluidRow(column(12,
					radioButtons(
						inputId="areaTypeRadio",
						label="Select type of genomic area",
						#inline=TRUE,
						choices=list(
							"Custom genomic area"="area",
							"Around gene"="gene"
						)
					)
				)),
				conditionalPanel(
					condition="input.areaTypeRadio=='area'",
					fluidRow(column(4,
						htmlOutput("setChrsA")
					),column(4,
						textInput(
							inputId="customStartA", 
							label="", 
							value="",
							placeholder="Start"
						)
					),column(4,
						textInput(
							inputId="customEndA",
							label="",
							value="",
							placeholder="End"
						)
					)),
					div(class="small",helpText(
						paste("Known genes overlaping the area above will be ",
							"automatically detected. Check the following box ",
							"if you want to add custom areas to the automatic ",
							"detection.")
					)),
					checkboxInput(
						inputId="addCustomRegionsInAreaCheck", 
						label="Include custom transcribed regions", 
						value=FALSE
					),
					conditionalPanel(
						condition="input.addCustomRegionsInAreaCheck",
						radioButtons(
							inputId="customRegionsFromGeneExplorer",
							label="Custom regions to include",
							choices=list(
								"From gene explorer (if any, same chromosome only)"="fromge",
								"Other custom regions"="other"
							)
						),
						conditionalPanel(
							condition=
								"input.customRegionsFromGeneExplorer=='other'",
							helpText("Please describe custom regions"),
							fluidRow(column(8,
								textInput(
									inputId="customNameInArea",
									label="", 
									value="",
									placeholder="Name"
								)
							),column(4,
								selectizeInput(
									inputId="customStrandInArea",
									label="",
									choices=c("","+","-"),
									options=list(placeholder="Strand")
								)
							)),
							fluidRow(column(4,
								disabled(textInput(
									inputId="customChrInArea",
									label="", 
									value="",
									placeholder="Chrom"
								))
							),column(4,
								textInput(
									inputId="customStartInArea", 
									label="", 
									value="",
									placeholder="Start"
								)
							),column(4,
								textInput(
									inputId="customEndInArea",
									label="",
									value="",
									placeholder="End"
								)
							)),
							fluidRow(column(12,
								div(class="small table-container",
									DT::dataTableOutput(
										"customRegionListInArea"
									)
								)
							)),
							fluidRow(column(12,
								htmlOutput("customRegionErrorInArea")
							)),
							fluidRow(br()),
							fluidRow(column(4,""
							),column(4,
								div(
									class="pull-right",
									style="display:inline-block",
									actionButton("removeCustomRegionInArea",
										"Remove",icon=icon("minus"))
								)
							),column(4,
								div(
									class="pull-right",
									style="display:inline-block",
									actionButton("addCustomRegionInArea",
										"Add",icon=icon("plus"))
								)
							))
						)
					)
				),
				conditionalPanel(
					condition="input.areaTypeRadio=='gene'",
					fluidRow(column(12,
						selectizeInput(
							inputId="areaGeneName",
							label="",
							multiple=FALSE,
							choices="",
							options=list(
								placeholder="Select a gene",
								selectOnTab=TRUE
							)
						)
					)),
					h5("Flanking",style="font-size:1.2em;"),
					fluidRow(column(6,
						textInput(
							inputId="upstreamFlankA", 
							label="Upstream", 
							value=2000
						)
					),column(6,
						textInput(
							inputId="downstreamFlankA",
							label="Downstream",
							value=2000
						)
					))
				),
				fluidRow(column(12,
					htmlOutput("areaFlankError")
				))
			)
		)),
		fluidRow(column(12,
			wellPanel(
				h4("Plot options"),
				fluidRow(column(12,
					htmlOutput("areaExplorerColours")
				)),
				fluidRow(br()),
				fluidRow(column(8,
					radioButtons(
						inputId="areaSumStatType",
						label="Area profile averaging",
						choices=list(
							"Mean"="mean",
							"Median"="median",
							"Trimmed mean"="trimmed"
						)
					)
				),column(4,
					textInput(
						inputId="areaTrimPct", 
						label="Trim fraction", 
						value="0.1"
					)
				)),
				fluidRow(column(8,
					htmlOutput("areaExplorerError")
				),column(4,
					div(
						class="pull-right",
						style="display:inline-block",
						actionButton(
							inputId="createAreaProfile",
							label="Engage!",
							icon=icon("rocket")
						)
					)
				))
			)
		))
	),column(9,
		fluidRow(column(12,
			plotOutput("areaProfile",height="640px")
		)),
		fluidRow(column(8,
			div("")
		),
		column(2,
			downloadButton(
				outputId="exportAreaPNG",
				label="Export PNG",
				#icon=icon("file-image-o"),
				class="pull-right"
			)
		),column(2,
			downloadButton(
				outputId="exportAreaPDF",
				label="Export PDF",
				#icon=icon("file-pdf-o"),
				class="pull-right"
			)
		)),
		fluidRow(br()),
		fluidRow(column(12,
			wellPanel(
				h4("Messages"),
				div(
					class="message-box",
					htmlOutput("areaExplorerMessages")
				)
			)
		))
	))
}
