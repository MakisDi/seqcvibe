makeFoldChange <- function(contrast,sample.list,data.matrix,scale="natural",
    log.offset=1) {
    conds <- strsplit(contrast,"_vs_")[[1]]
    fold.mat <- matrix(0,nrow(data.matrix),length(conds)-1)
    for (i in 1:(length(conds)-1)) { # Last condition is ALWAYS reference
        samples.nom <- sample.list[[conds[i]]]
        samples.denom <- sample.list[[conds[length(conds)]]]
        nom <- data.matrix[,match(samples.nom,colnames(data.matrix)),drop=FALSE]
        denom <- data.matrix[,match(samples.denom,colnames(data.matrix)),
            drop=FALSE]
        if (!is.matrix(nom)) 
            nom <- as.matrix(nom) # Cover the case with no replicates...
        if (!is.matrix(denom)) 
            denom <- as.matrix(denom)
        mean.nom <- apply(nom,1,mean)
        mean.denom <- apply(denom,1,mean)
        if (any(mean.nom==0)) 
            mean.nom <- mean.nom + log.offset
        if (any(mean.denom==0)) 
            mean.denom <- mean.denom + log.offset
        if (scale=="natural")
            fold.mat[,i] <- mean.nom/mean.denom
        else if (scale=="log2")
            fold.mat[,i] <- mean.nom - mean.denom
    }
    rownames(fold.mat) <- rownames(data.matrix)
    colnames(fold.mat) <- paste(conds[1:(length(conds)-1)],"_vs_",
        conds[length(conds)],sep="")
    return(fold.mat)
}

makeA <- function(contrast,sample.list,data.matrix,log.offset=1) {
    conds <- strsplit(contrast,"_vs_")[[1]]
    a.mat <- matrix(0,nrow(data.matrix),length(conds)-1)
    for (i in 1:(length(conds)-1)) { # Last condition is ALWAYS reference
        samples.trt <- sample.list[[conds[i]]]
        samples.cnt <- sample.list[[conds[length(conds)]]]
        trt <- data.matrix[,match(samples.trt,colnames(data.matrix)),drop=FALSE]
        cnt <- data.matrix[,match(samples.cnt,colnames(data.matrix)),drop=FALSE]
        if (!is.matrix(trt)) 
            trt <- as.matrix(trt) # Cover the case with no replicates...
        if (!is.matrix(cnt)) 
            cnt <- as.matrix(cnt)
        mean.trt <- apply(trt,1,mean)
        mean.cnt <- apply(cnt,1,mean)
        if (any(mean.trt==0)) 
            mean.trt <- mean.trt + log.offset
        if (any(mean.cnt==0)) 
            mean.cnt <- mean.cnt + log.offset
        a.mat[,i] <- 0.5*(log2(mean.trt)+log2(mean.cnt))
    }
    rownames(a.mat) <- rownames(data.matrix)
    colnames(a.mat) <- paste(conds[1:(length(conds)-1)],"_vs_",
        conds[length(conds)],sep="")
    return(a.mat)
}

makeStat <- function(samples,data.mat,stat,value) {    
    stat.data <- data.mat[,match(samples,colnames(data.mat))]
    if (!is.matrix(stat.data)) 
        stat.data <- as.matrix(stat.data)
    if (value=="counts") {
        switch(stat,
            mean = {
                return(apply(stat.data,1,function(x) {
                    return(as.integer(round(mean(x))))
                }))
            },
            median = {
                stat.calc <- apply(stat.data,1,function(x) {
                    return(as.integer(round(median(x))))
                })
            },
            sd = {
                return(apply(stat.data,1,function(x) {
                    return(as.integer(round(sd(x))))
                }))
            },
            mad = {
                return(apply(stat.data,1,function(x) {
                    return(as.integer(round(mad(x))))
                }))
            },
            cv = {
                return(apply(stat.data,1,function(x,s) {
                    return(round(sd(x))/round(mean(x)))
                }))
            },
            rcv = {
                return(apply(stat.data,1,function(x,s) {
                    return(round(mad(x))/round(median(x))) 
                }))
            }
        )
    }
    else {
        switch(stat,
            mean = {
                return(apply(stat.data,1,function(x) {
                    return(round(mean(x),6))
                }))
            },
            median = {
                stat.calc <- apply(stat.data,1,function(x) {
                    return(round(median(x),6))
                })
            },
            sd = {
                return(apply(stat.data,1,function(x) {
                    return(round(sd(x),6))
                }))
            },
            mad = {
                return(apply(stat.data,1,function(x) {
                    return(round(mad(x),6))
                }))
            },
            cv = {
                return(apply(stat.data,1,function(x,s) {
                    return(round(sd(x)/mean(x),6))
                }))
            },
            rcv = {
                return(apply(stat.data,1,function(x,s) {
                    return(round(mad(x)/median(x),6))
                }))
            }
        )
    }
}

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

getBiotypes <- function(org) {
    if (!(org %in% c("hg18","hg19","hg38","mm9","mm10","rn5","dm3","danrer7",
        "pantro4","susscr3")))
        return(NULL)
    switch(org,
        hg18 = {
            return(c("unprocessed_pseudogene","pseudogene","miRNA",
                "retrotransposed","protein_coding","processed_pseudogene",
                "snRNA","snRNA_pseudogene","Mt_tRNA_pseudogene",
                "miRNA_pseudogene","misc_RNA","tRNA_pseudogene","snoRNA",
                "scRNA_pseudogene","rRNA_pseudogene","snoRNA_pseudogene","rRNA","misc_RNA_pseudogene","IG_V_gene","IG_D_gene","IG_J_gene",
                "IG_C_gene","IG_pseudogene","scRNA"))
        },
        hg19 = {
            return(c("pseudogene","lincRNA","protein_coding","antisense",
                "processed_transcript","snRNA","sense_intronic","miRNA",
                "misc_RNA","snoRNA","rRNA","polymorphic_pseudogene",
                "sense_overlapping","3prime_overlapping_ncrna","TR_V_gene",
                "TR_V_pseudogene","TR_D_gene","TR_J_gene","TR_C_gene",
                "TR_J_pseudogene","IG_C_gene","IG_C_pseudogene","IG_J_gene",
                "IG_J_pseudogene","IG_D_gene","IG_V_gene","IG_V_pseudogene"))
        },
        hg38 = {
            return(c("protein_coding","polymorphic_pseudogene","lincRNA",
                "unprocessed_pseudogene","processed_pseudogene","antisense",
                "processed_transcript","transcribed_unprocessed_pseudogene",
                "sense_intronic","unitary_pseudogene","IG_V_gene",
                "IG_V_pseudogene","TR_V_gene","sense_overlapping",
                "transcribed_processed_pseudogene","miRNA","snRNA","misc_RNA",
                "rRNA","snoRNA","IG_J_pseudogene","IG_J_gene","IG_D_gene",
                "3prime_overlapping_ncrna","IG_C_gene","IG_C_pseudogene",
                "pseudogene","TR_V_pseudogene","Mt_tRNA","Mt_rRNA",
                "translated_processed_pseudogene","TR_J_gene","TR_C_gene",
                "TR_D_gene","TR_J_pseudogene","LRG_gene"))
        },
        mm9 = {
            return(c("pseudogene","snRNA","protein_coding","antisense","miRNA",
                "lincRNA","snoRNA","processed_transcript","misc_RNA","rRNA",
                "sense_overlapping","sense_intronic","polymorphic_pseudogene",
                "non_coding","3prime_overlapping_ncrna","IG_C_gene",
                "IG_J_gene","IG_D_gene","IG_V_gene","ncrna_host"))
        },
        mm10 = {
            return(c("pseudogene","snRNA","protein_coding","antisense","miRNA",
                "snoRNA","lincRNA","processed_transcript","misc_RNA","rRNA",
                "sense_intronic","sense_overlapping","polymorphic_pseudogene",
                "IG_C_gene","IG_J_gene","IG_D_gene","IG_LV_gene","IG_V_gene",
                "IG_V_pseudogene","TR_V_gene","TR_V_pseudogene",
                "3prime_overlapping_ncrna"))
        },
        dm3 = {
            return(c("protein_coding","ncRNA","snoRNA","pre_miRNA","pseudogene",
                "snRNA","tRNA","rRNA"))
        },
        rn5 = {
            return(c("protein_coding","pseudogene","processed_pseudogene",
                "miRNA","rRNA","misc_RNA"))
        },
        danrer7 = {
            return(c("antisense","protein_coding","miRNA","snoRNA","rRNA",
                "lincRNA","processed_transcript","snRNA","pseudogene",
                "sense_intronic","misc_RNA","polymorphic_pseudogene",
                "IG_V_pseudogene","IG_C_pseudogene","IG_J_pseudogene",
                "non_coding","sense_overlapping"
            ))
        },
        pantro4 = {
            return(c("protein_coding","pseudogene","processed_pseudogene",
                "miRNA","rRNA","snRNA","snoRNA","misc_RNA"))
        },
        susscr3 = {
            return(c("antisense","protein_coding","lincRNA","pseudogene",
                "processed_transcript","miRNA","rRNA","snRNA","snoRNA",
                "misc_RNA","non_coding","IG_C_gene","IG_J_gene",
                "IG_V_gene","IG_V_pseudogene"))
        },
        tair10 = {
            return(c("miRNA","ncRNA","protein_coding","pseudogene","rRNA",
                "snoRNA","snRNA","transposable_element","tRNA"))
        }
    )
}

getTime <- function(type) {
    tt <- format(Sys.time(),format="%Y-%m-%d %H:%M:%S")
    return(paste(type," ",tt,": ",sep=""))
}

isEmpty <- function(x) {
    return(is.null(x) || x=="")
}
