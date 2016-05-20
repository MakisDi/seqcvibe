# Function initializing SeqCVIBE universe (packages, persistent variables, etc.)
initUniverse <- function(session) {
	# Initial page loading indicator, until all content is loaded
	ftProgress <- shiny::Progress$new()
	ftProgress$initialize(session,min=0,max=10)
	ftProgress$set(message="Starting:",value=0)
	on.exit(ftProgress$close())
	# Progress update function
	updateFtProgress <- function(value=NULL,detail=NULL) {
		if (is.null(value)) {
			value <- ftProgress$getValue()
			value <- value + 1
		}
		ftProgress$set(value=value,detail=detail)
	}
	
	# Load packages
	updateFtProgress(value=1,detail="Loading DT")
	require(DT)
	updateFtProgress(value=2,detail="Loading ggplot2")
	require(ggplot2)
	updateFtProgress(value=3,detail="Loading GenomicRanges")
	require(GenomicRanges)
	updateFtProgress(value=4,detail="Loading GenomicAlignments")
	require(GenomicAlignments)
	updateFtProgress(value=5,detail="Loading rtracklayer")
	require(rtracklayer)
	updateFtProgress(value=6,detail="Loading ggbio")
	require(ggbio)
	updateFtProgress(value=7,detail="Loading d3heatmap")
	require(d3heatmap)
	#require(plotly)

	# Load additional functions
	updateFtProgress(value=8,detail="Loading SeqCVIBE lib")
	source("lib/control.R")
	source("lib/util.R")

	# Load metadata
	updateFtProgress(value=9,detail="Loading metadata")
	metadata <- read.delim("config/metadata.txt")
	rownames(metadata) <- as.character(metadata$sample_id)

	# Intialize metadata reactive content
	sources <- unique(as.character(metadata$source))
	datasets <- unique(as.character(metadata$dataset[
		which(as.character(metadata$source)==sources[1])]))
	classes <- unique(as.character(metadata$class[
		which(as.character(metadata$source)==sources[1] 
			& as.character(metadata$dataset)==datasets[1])]))
	genomes <- unique(as.character(metadata$genome))
	genome <- genomes[1]

	# Load data file hash
	source("config/data_files.R")
	allClasses <- unique(as.character(metadata$class))
	baseColours <- c("#B40000","#00B400","#0000B4","#B45200","#9B59B6","#21BCBF",
		"#BC4800","#135C34","#838F00","#4900B5")
	baseColours <- rep(baseColours,length.out=length(allClasses))
	names(baseColours) <- allClasses

	# Keep track of loaded annotations and load the first
	updateFtProgress(value=10,detail="Loading initial genome")
	load(file.path("genome",genome,"gene.rda"))
	loadedGenomes <- vector("list",length(genomes))
	names(loadedGenomes) <- genomes
	for (gen in genomes) {
		loadedGenomes[[gen]] <- list(
			geneNames=NULL,
			dbGene=NULL,
			dbExon=NULL
		)
	}
	geneNames <- names(gene)
	names(geneNames) <- as.character(elementMetadata(gene)$gene_name)
	loadedGenomes[[genome]] <- list(
		geneNames=geneNames,
		dbGene=gene,
		dbExon=NULL
	)

	# Keep track of loaded data and load the first
	loadedData <- vector("list",length(sources))
	names(loadedData) <- sources
	for (s in sources) {
		dd <- unique(as.character(metadata$dataset[
			which(as.character(metadata$source)==s)]))
		loadedData[[s]] <- vector("list",length(dd))
		names(loadedData[[s]]) <- dd
	}

	# Restrict the number of cores dedicated to SeqCVIBE
	#RC <- 0.25
	RC <- NULL

	# Hack for p-value, FDR display in sliders and natural fold change
	statScoreValues <- c(0,1e-16,1e-8,1e-4,0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.1,
		0.2,0.3,0.4,0.5,1)
	fcNatValues <- c(0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,
		0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10)
		
	return(list(
		metadata=metadata,
		sources=sources,
		datasets=datasets,
		classes=classes,
		genomes=genomes,
		baseColours=baseColours,
		loadedGenomes=loadedGenomes,
		loadedData=loadedData,
		statScoreValues=statScoreValues,
		fcNatValues=fcNatValues,
		geneNames=geneNames,
		RC=RC
	))
}
