# This R script is used to download the data in the matrix format from GEO and then map probesets to genes based
# on the platform annotation. The major functions in this R script are refactored from Nasim's original harmonize.R

library(GEOquery)
library(GEOmetadb)
library(limma)
library(logger)
library(stringr)
library(ggplot2)
library(autoplotly)

# Have to make sure this sqlite database exist to collect data.
if ( !file.exists( "GEOmetadb.sqlite" ) ) {
    geometadbfile <- getSQLiteFile(destdir = ".") # download full DB
} else {
    geometadbfile <- "GEOmetadb.sqlite"
}

set.seed(1234)
options(stringsAsFactors = F)
log_threshold(DEBUG)

# Process the downloaded and export the processed data into a dest directory.
process.gse <- function(gse.id, 
                        dest.dir,
                        need.plot = TRUE, # Check if we need to have output
                        gsm.ids = NA) {
    # Download the data
    # Create a temp folder to download the data if needed
    temp.dir <- paste(dest.dir, "temp", sep = "/")
    if (!dir.exists(temp.dir)) {
        dir.create(temp.dir)
    }
    gse.data <- getGEO(gse.id, destdir = temp.dir, AnnotGPL = TRUE)
    # The results are organized as GPL
    genes.expression <- data.frame()
    # Collect gpl information for batch correction pca analysis
    # This basically is a map from gsm to gpl
    gsms.gpls <- data.frame()
    for (i in 1 : length(gse.data)) {
        gpl.gse.data <- gse.data[[i]]
        gpl.gse.gene.expression <- process.gpl.data(gpl.gse.data,
                                                    gsm.ids)
        if (is.null(gpl.gse.gene.expression)) {
            log_info(paste("No GSM left aftering filter for", gpl.gse.data@annotation))
            next
        }
        gsms <- colnames(gpl.gse.gene.expression)
        gsms.gpls <- rbind(gsms.gpls, data.frame(gsm=gsms, gpl=gpl.gse.data@annotation))
        # Rows are genes and columns are gsms
        genes.expression <- merge(genes.expression, gpl.gse.gene.expression, by = 0, all = T)
        # After the above call, there is a new row.names added and the row.names is gone.
        # Need to do something about it
        rownames(genes.expression) <- genes.expression$Row.names
        genes.expression$Row.names <- NULL # Remove this row
    }
    if (dim(genes.expression)[1] == 0) {
        log_info(paste("An empty gene expression data: ", gse.id, sep = ""))
        return(NA)
    }
    # Dump the data 
    file.name <- paste(dest.dir, paste(gse.id, ".txt", sep = ""), sep = "/")
    write.table(genes.expression,
                file = file.name,
                quote = FALSE,
                sep = "\t")
    if (need.plot) {
        file.name <- paste(dest.dir, paste(gse.id, ".pdf", sep=""), sep="/")
        eda(genes.expression, gsms.gpls, file.name)
    }
    genes.expression
}

# Process a GPL specific GSE data set.
process.gpl.data <- function(gpl.gse.data, 
                             gsm.ids = NA,
                             gene = NA) {
    # This function is used to do mapping from genes to probesets
    # expressions should be a dataframe
    expressions <- gpl.gse.data@assayData$exprs
    # Check if we need to do a filtering
    if (!is.null(gsm.ids)) { # Use is.null for this list
        which.cols <- colnames(expressions) %in% gsm.ids
        old.size <- dim(expressions)[2]
        expressions <- expressions[, which.cols]
        new.size <- dim(expressions)[2]
        diff <- old.size - new.size;
        if (diff > 0) {
            log_info(paste("Some GSM ids are removed: ", diff, sep = ""))
        }
        # If there is nothing left, just return na
        if (dim(expressions)[2] == 0) {
            return(NULL)
        }
    }
    # Make sure gpl.gse.data is ExpressionSet
    if (!is(gpl.gse.data, "ExpressionSet")) {
        stop("The passed argument should be an ExpressionSet object.")
    }
    # Ignore probsets that cannot be mapped to genes
    features <- gpl.gse.data@featureData@data
    # Upper case the col names for easy handling
    colnames(features) <- toupper(colnames(features))
    if (is.na(gene)) {
        # Need to determine what is the gene symbol col
        if ("GENE SYMBOL" %in% colnames(features)) {
            symbol.col.name <- "GENE SYMBOL" # Affy
        }else if ("SYMBOL" %in% colnames(features)) {
            symbol.col.name <- "SYMBOL" # Illumina
        }else {
            map.results <- map.to.genes(features)
            if (is.null(map.results))
                stop("Cannot find the column for gene symbols in the features!")
            symbol.col.name <- map.results$'gene.col.name'
            features <- map.results$'features' # This hould be modified with a new column
        }
        genes <- features[, symbol.col.name]
        # Remove rows that don't have genes
        genes <- genes[which(genes != "")]
        # Some probes can be mapped to multiple genes, esp. in Affy
        # Many of genes there have their own unique probes. In this version,
        # these genes will be ingored. For example, IL28A /// IL28B
        # IL28A is mapped to another probe. However, IL28B is not. We may need to 
        # tolerate the information loss a little bit here. Should be less than 3%!
        which <- grep("///", genes) # e.g. IL28A /// IL28B. 
        if (length(which) > 0) {
            genes <- genes[-which]
        }
        # Remove duplicated genes
        genes <- unique(genes)
        # genes <- sample(genes, 100)
    }else {
        genes <- c(gene)
    }
    map.probe.to.genes <- lapply(genes, function(gene){
        probes <- features$ID[features[, symbol.col.name] %in% gene]
        which.probes <- rownames(expressions) %in% probes
        if (sum(which.probes) > 0) {
            gene.expression <- expressions[which.probes, ]
            if (sum(which.probes) > 1) {
                row <- apply(gene.expression, 2, median)
            }else {
                row <- gene.expression
            }
            row.df <- data.frame(t(row))
            rownames(row.df) <- c(gene)
            return(row.df)
        }
    })
    # Create a new dataframe
    genes.expressions <- data.frame(do.call(rbind, map.probe.to.genes))
    genes.expressions
}

# If features don't have a gene symbol column, we will do our best to map to
# genes using external file
map.to.genes <- function(features) {
    # First check if UniGene column exists
    gene.col.name <- 'GENE SYMBOL'
    features.col.names <- colnames(features)
    # Two colnames for mapping
    feature.col.name <- NA
    map.col.name <- NA
    if ("UNIGENE" %in% features.col.names) {
        log_info("Using UniGene to map gene symbols")
        feature.col.name <- "UNIGENE"
        map.col.name <- "UNIGENE"
        
    }else if("ENTREZ_GENE_ID" %in% features.col.names) {
        log_info("Using ENTREZ_GENE_ID to map gene symbols")
        feature.col.name <- "ENTREZ_GENE_ID"
        map.col.name <- "GENE_ID"
    }else if ("LOCUSLINK" %in% features.col.names) {
        log_info("Usiing LOCUSLINK to map gene symbols")
        feature.col.name <- "LOCUSLINK"
        map.col.name <- "GENE_ID"
    }
    if (is.na(feature.col.name) || is.na(map.col.name)) return(NULL)
    # Have to make sure this file exists
    unigene.2.gene <- read.delim(unigene.gene.map.file, header = T, check.names = F)
    colnames(unigene.2.gene) <- toupper(colnames(unigene.2.gene))
    indices <- match(features[, feature.col.name], unigene.2.gene[, map.col.name])
    log_info(paste("Total mapped genes:", sum(!is.na(indices))))
    features[, gene.col.name] <- unigene.2.gene[, gene.col.name][indices]
    rtn <- list("gene.col.name" = gene.col.name, "features" = features)
    # View(features)
    return(rtn)
}

# Perform some simple EDA and save the results into a folder.
eda <- function(expressions, 
                gsms.gpls, # Provide a map from gsms to gpls
                file.name) {
    # Make sure file.name ends with .pdf
    if (!str_ends(file.name, ".pdf")) {
        file.name = paste(file.name, ".pdf", sep = "")
    }
    # pca.plot <- ggplot(plotdata, aes(x=PC1, y=PC2)) + geom_point()
    # print(pca.plot)
    pdf(file = file.name,  wi = 12, he = 9)
    # Use genes having all values
    clean.expressions <- expressions[complete.cases(expressions), ]
    # Need this df for the following plots
    pca(clean.expressions)
    # To see if we need to do a batch correction
    pca.with.batch.correction(clean.expressions, gsms.gpls)
    plot.mad.mean(expressions)
    dev.off()
}

plot.mad.mean <- function(expressions) {
    df <- expressions[sample(rownames(expressions), 1000), ]
    # The following should be gene-wise information
    # Need to have a better control
    mad.values <- c()
    mean.values <- c()
    for (i in 1 : dim(df)[1]) {
        row <- df[i, ]
        row <- row[!is.na(row)]
        if (length(row) == 0) next
        mad.values <- c(mad.values, mad(row))
        mean.values <- c(mean.values, mean(row))
    }
    plot(sort(mad.values))
    plot(sort(mean.values))
    smoothScatter(x = mean.values, y = mad.values)
}

pca <- function(expressions,
                title = "",
                need.no.scale = TRUE) {
    df <- prepare.pca.df(expressions)
    if (is.null(df)) {
        log_info("Cannot generate a DataFrame for PCA.")
        return(NULL)
    }
    pca <- prcomp(df, scale.=T, center = T) # Use whatever data we have
    # Plot first 2 PCs
    plotdata <- data.frame(SampleID=rownames(pca$x),
                           PC1=pca$x[,1],
                           PC2=pca$x[,2])
    pca.plot <- ggplot(plotdata, aes(x=PC1, y=PC2)) + geom_point() +
        labs(title = paste(title, "PCA with Scale and Center", sep=""))
    print(pca.plot)
    if (!need.no.scale) {
        return(pca)
    }
    pca <- prcomp(df, scale.=F, center = F) # Use whatever data we have
    # Plot first 2 PCs
    plotdata <- data.frame(SampleID=rownames(pca$x),
                           PC1=pca$x[,1],
                           PC2=pca$x[,2])
    pca.plot <- ggplot(plotdata, aes(x=PC1, y=PC2)) + geom_point() + 
        labs(title = paste(title, "PCA without Scale and Center", sep = ""))
    print(pca.plot)
    pca
}

# Generate a dataframe for PCA
prepare.pca.df <- function(expressions) {
    if (dim(expressions)[1] == 0) {
        log_info("No genes are in the expressions dataframe for PCA analysis.")
        return(NULL)
    }
    df <- expressions[sample(rownames(expressions), 1000), ]
    # Make sure genes are columns and samples are rows
    df <- t(df)
    # To avoid error in the following prcomp, pick up genes having enough variance
    vars <- apply(df, 2, var)
    # This is arbirary number
    which <- vars > 1.0E-5
    df <- df[, which]
    # Remove columns having NA: 
    # https://stackoverflow.com/questions/2643939/remove-columns-from-dataframe-where-all-values-are-na
    df <- df[, colSums(is.na(df)) < nrow(df)]
    df
}

# Use autoplotly to plot an interactive PCA for manual investigation. The code is modified
# from this web page: https://terrytangyuan.github.io/2018/02/12/autoplotly-intro/
interactive.pca <- function(expressions,
                            meta.file) {
    df <- prepare.pca.df(expressions)
    pca <- prcomp(df, scale.=T, center = T) # Use whatever data we have
    # Want to get some meta information to color the plot
    meta <- read.delim(meta.file, sep = "\t", header = T, check.names = T)
    which <- meta[, 'Expsample.Repository.Accession'] %in% rownames(df)
    meta.rows <- meta[which, ]
    tmp <- meta.rows
    # Sort based on df
    indices <- match(rownames(df), meta.rows[, 'Expsample.Repository.Accession'])
    meta.rows <- meta.rows[indices, ]
    rownames(meta.rows) <- meta.rows[, 'Expsample.Repository.Accession']
    # Use the call from autoplotly
    plotdata <- autoplotly(pca, frame = TRUE, data = meta.rows, 
                           label = TRUE, label.size = 2,
                           originalData = FALSE,
                           colour = "Expsample.Repository.Accession")
    print(plotdata)
    pca
}

pca.with.batch.correction <- function(expressions, gsms.gpls) {
    if (dim(expressions)[1] == 0) {
        log_info("No genes are in the expressions dataframe for batch analysis.")
        return(NULL)
    }
    gpls <- unique(gsms.gpls$gpl)
    if (length(gpls) == 1) {
        return() # There is no need to do
    }
    log_info("Doing GPL correction...")
    # Usually this should not be needed. Just to make sure
    gsms.gpls <- gsms.gpls[match(colnames(expressions), gsms.gpls$gsm), ]
    corrected.expressions <- removeBatchEffect(expressions, gsms.gpls$gpl)
    pca(corrected.expressions, title = "GPL Corrected ")
}

# Process a set of GSE accessions from the passed file
process.file <- function(file.name, 
                         dest.dir,
                         gse.id = NA,
                         gse.col.name = 'GSE Accession',
                         gsm.col.name = 'Expsample Repository Accession') {
    immport.data <- read.delim(file.name, header = T, check.names = F)
    gse.accessions <- unique(immport.data[, gse.col.name])
    log_debug(paste("Total GSEs: ", length(gse.accessions), sep = ""))
    gpls <- unique(immport.data[, "GPL Accession"])
    log_debug(paste("Total GPL:", length(gpls)))
    log_debug()
    for (i in 1:length(gse.accessions)) {
        gse <- gse.accessions[i]
        if (!is.na(gse.id) && gse != gse.id) { # For a test
            next
        }
        log_debug(paste("Processing ", gse, sep = ""))
        which <- immport.data[, gse.col.name] == gse
        gsms <- immport.data[, gsm.col.name][which]
        log_debug(paste("Total GSMs: ", length(gsms), sep = ""))
        process.gse(gse.id = gse,
                    dest.dir = dest.dir,
                    need.plot = T,
                    gsm.ids = gsms)
    }
}

# The following variables are used as a configuration and make sure they exist
unigene.gene.map.file <- "UniGeneMapper_071321.txt"
if (!file.exists(unigene.gene.map.file)) {
    stop(paste(unigene.gene.map.file, "doesn't exist!"))
}
# geneacc.gene.mal.file <- "/Volumes/ssd/datasets/NCBI/UniGene/UniGeneMapperWithAcc_071321.txt"
# if (!file.exists(geneacc.gene.mal.file)) {
#     stop(paste(geneacc.gene.mal.file, "doesn't exist!"))
# }
# This is for test
dest.dir <- "output"
file.name <- "ImmuneExposureGeneExpression.txt"
gse.id <- "GSE13485" # Test GPL7567 using UniGene map
gse.id <- "GSE22768" # Test GPL10647 using ENTREZ_GENE_ID map
# gse.id <- "GSE22121" # Test GPL9700 and GPL10465 using LOCUSLINK: There are only 
#                     #  about 7,000 columns in the expression matrix file though much more rows in the GPL
gse.id <- NA # To run all GSEs
# gse.id <- "GSE47353" # There are two files related to the platform used by this GSE: GPL6244.
gse.id <- "GSE52005" # Error related if diff > 0
gse.id <- "GSE30101" # Same platform but three batches
gse.id <- "GSE59714"
log_info("Starting processing the file...")
# process.file(file.name, dest.dir, gse.id)
log_info("Done!")