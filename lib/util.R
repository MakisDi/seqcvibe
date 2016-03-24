updateMessages <- function(messageContainer,type,msg,clear=FALSE) {
    if (clear)
        messageContainer$messages <- list(
            type=type,
            msg=msg
        )
    else {
        n <- length(messageContainer$messages)
        messageContainer$messages[[n+1]] <- list(type=type,msg=msg)
    }
    return(messageContainer)
}

bsvMessage <- function(...,messageContainer=NULL) {
    if (is.null(messageContainer))
        message(...)
    else if (is.character(messageContainer)) {
        #print(messageContainer)
        write(paste(...,sep=""),file=messageContainer,append=TRUE)
    }
    else {              
        messageContainer <- updateMessages(
            messageContainer,
            type="INFO",
            msg=paste(getTime("INFO"),...,sep="")
        )
    }
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        m <- FALSE
    else {
        m <- TRUE
        ncores <- parallel::detectCores()
        if (ncores==1) 
            m <- FALSE
        else {
            if (!missing(rc) && !is.null(rc))
                ncores <- ceiling(rc*ncores)
        }
    }
    if (m)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}

ssCI <- function(fit) {
    res <- (fit$yin - fit$y)/(1-fit$lev)
    sigma <- sqrt(var(res))
    upper <- fit$y + 3.0*sigma*sqrt(fit$lev)
    lower <- fit$y - 3.0*sigma*sqrt(fit$lev)
    return(list(lower=lower,upper=upper))
}

splitBySeqname <- function(gr,verbose=FALSE,rc=NULL) {
    if (verbose)
        message("Splitting input regions by seqname...")
    gr.list <- cmclapply(levels(seqnames(gr)),function(x,lib) {
        if (verbose)
            message("  ",x)
        tmp <- lib[seqnames(lib)==x]
        if (length(tmp)>0) return(tmp) else return(NULL)
    },gr,rc=rc)
    names(gr.list) <- levels(seqnames(gr))
    null <- which(sapply(gr.list,is.null))
    if (length(null)>0)
        gr.list <- gr.list[-null]
    return(gr.list)
}

splitVector <- function(x,n,seed=42) {
    bin.size <- floor(length(x)/n)
    dif <- length(x) - bin.size*n 
    bin.fac <- rep(bin.size,n)
    # Random bin increase size to avoid problems
    set.seed(seed)
    add <- sample(1:n,dif)
    bin.fac[add] <- bin.fac[add]+1
    f <- factor(rep(1:n,bin.fac))
    return(split(x,f))
}

checkConfig <- function(config) {
    if (!all(c("sample_id","sample_dir") %in% names(config)))
        stop("Config must contain at least sample ids and dirs")
    n <- nrow(config)
    if (!("source") %in% names(config))
        config$source <- rep("Unknown source",n)
    if (!("dataset") %in% names(config))
        config$dataset <- rep("Unknown dataset",n)
    if (!("class") %in% names(config))
        config$class <- rep("Unknown class",n)
    if (!("alt_id") %in% names(config))
        config$alt_id <- rep("N/A",n)
    if (!("norm_factor") %in% names(config))
        config$norm_factor <- rep(1,n)
    if (!("library_strategy") %in% names(config))
        config$library_strategy <- rep("N/A",n)
    if (!("quality") %in% names(config))
        config$quality <- rep("N/A",n)
    return(config)
}

getValidChromosomes <- function(org) {
    switch(org,
        hg18 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg19 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg38 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm9 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        rn5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX"
            ))
        },
        dm3 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        danRer7 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        panTro4 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        susScr3 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY"
            ))
        },
        tair10 = {
            return(c(
                "chr1","chr2","chr3","chr4","chr5"
            ))
        }
    )
}

getFriendlyName <- function(org) {
    switch(org,
        hg18 = { return("Human") },
        hg19 = { return("Human") },
        hg38 = { return("Human") },
        mm9 = { return("Mouse") },
        mm10 = { return("Mouse") },
        rn5 = { return("Rat") },
        dm3 = { return("Fruitfly") },
        danrer7 = { return("Zebrafish") },
        pantro4 = { return("Chimpanzee") },
        susscr3 = { return("Pig") },
        tair10 = { return("Arabidopsis") }
    )
}

getOrganismName <- function(org) {
    switch(org,
        hg18 = { return("Homo sapiens") },
        hg19 = { return("Homo sapiens") },
        hg38 = { return("Homo sapiens") },
        mm9 = { return("Mus musculus") },
        mm10 = { return("Mus musculus") },
        rn5 = { return("Rattus norvegicus") },
        dm3 = { return("Drosophila melanogaster") },
        danrer7 = { return("Danio rerio") },
        pantro4 = { return("Pan troglodytes") },
        susscr3 = { return("Sus scrofa") },
        tair10 = { return("Arabidopsis thaliana") }
    )
}

getTime <- function(type) {
    tt <- format(Sys.time(),format="%Y-%m-%d %H:%M:%S")
    return(paste(type," ",tt,": ",sep=""))
}

isEmpty <- function(x) {
    return(is.null(x) || x=="")
}
