# -----------------------------------------------------
# Read absolute path inputs for expression, phenotype, and user-selections files
# Specify output absolute path and file name as the last arg to have generated results to be saved in
# -----------------------------------------------------
# This script is modified from differential_genes.R

# -----------------------------------------------------
usePackage("pacman")
# -----------------------------------------------------
# load libraries 
# -----------------------------------------------------
p_load("jsonlite")
p_load("tidyverse")
p_load("dplyr")
p_load("doParallel")
p_load("limma")

# The file paths for pre-generated files are hard-coded here. This may be enhanced so that these file names
# can be configured inside Java directly

dir <- args[2]
pheno.path <- paste(dir, "biosample_metadata2.csv", sep = "/")
all.exp.path <- paste(dir, "all_expr_df_final.csv", sep = "/")
all.exp.adjusted.path <- paste(dir, "all_expr_df_final_adjusted.csv", sep = "/")

# ------ Load data and then cached them -----------------------------------------------------
set.seed(1234)
options(stringsAsFactors = F)
# registerDoParallel(detectCores() - 2) # optional for back-end speed 
# -----------------------------------------------------
# read in data for pre-processed biosample expression and associated phenotypes 
# -----------------------------------------------------
all.pheno.dat <- read.csv(pheno.path) %>% as.data.frame()
# We will do dynamic loading for these two big data files
all.exp.data <- NULL
all.exp.adjusted.data <- NULL

load_exp_data <- function(file.name) {
  all.expr.dat <- read_csv(all.exp.path) %>% as.data.frame()
  rownames(all.expr.dat) <- all.expr.dat[, 1]
  all.expr.dat <- all.expr.dat[, -1]
}


# Perform the differential expression analysis via limma
do_diff_exp_analysis <- function(selection.json.text) {
  
  user.select <- fromJSON(selection.json.text, simplifyVector = TRUE)
  
  # -----------------------------------------------------
  # Filter biosamples by users study design selections &
  # make sure biosamples exist & have matching orders 
  # -----------------------------------------------------
  pheno.dat <- all.pheno.dat[which(all.pheno.dat$gsm %in% user.select$GSMids), ]
  # Based on the user's choice
  expr.dat <- NULL
  if (user.select$platformCorrection) {
    if (is.null(all.exp.adjusted.data)) {
      all.exp.adjusted.data <<- load_exp_data(all.exp.adjusted.path)
    }
    expr.dat <- all.exp.adjusted.data[ ,which(colnames(all.exp.adjusted.data) %in% pheno.dat$gsm)]
  }
  else {
    if (is.null(all.exp.data)) { # Use the expression data directly
      all.exp.data <<- load_exp_data(all.exp.path)
    }
    expr.dat <- all.exp.data[ ,which(colnames(all.exp.data) %in% pheno.dat$gsm)]
  }
  # Make sure it has the same order as in the sample meta data. Not sure why (GW)?
  expr.dat <- expr.dat[ ,match(pheno.dat$gsm, colnames(expr.dat))]
  
  # -----------------------------------------------------
  # make factors for groups selected 
  # -----------------------------------------------------
  # one of the study variables selection lengths has to be greater than 2 
  # treat groups as factors 
  # -----------------------------------------------------
  if (length(user.select$studyCohort) == 1){
    pheno.dat[,which(colnames(pheno.dat) %in% user.select$studyCohort)] <- as.factor(pheno.dat[,which(colnames(pheno.dat) %in% user.select$studyCohort)])
  } else {
    for (s in user.select$studyCohort){
      pheno.dat[,which(colnames(pheno.dat) %in% s)] <- as.factor(pheno.dat[,which(colnames(pheno.dat) %in% s)])
    }
  }
  
  # -----------------------------------------------------
  # Filter variable genes for stronger analysis results 
  # based on users request 
  # Note: 0.6 seems to be a good threshold for this slice of data 
  # -----------------------------------------------------
  if (user.select$variableGenes){
    mdf.mad <- apply(expr.dat, 1, mad)
    variable.genes <- names(mdf.mad[mdf.mad >  0.6])
    
    expr.dat <- expr.dat[which(rownames(expr.dat) %in% variable.genes), ]
  }
  # -----------------------------------------------------
  # conduct differential expression analysis
  # add time interactions 
  # -----------------------------------------------------
  if (user.select$modelTime){
    # ex. ~ age_group + immport_vaccination_time + age_group:immport_vaccination_time
    synergy.terms <- lapply(user.select$studyCohort, function(s){
      paste(s, "immport_vaccination_time", sep =':')
    })
    synergy.terms <- unlist(synergy.terms)
    
    de.formula <- as.formula(paste("~" , paste(c(user.select$studyCohort, "immport_vaccination_time", synergy.terms), collapse = ' + ')))
    
  } else {
    # ex. ~ age_group + immport_vaccination_time_groups + age_group:immport_vaccination_time_groups
    pheno.dat$immport_vaccination_time_groups <- NA
    pheno.dat$immport_vaccination_time_groups[pheno.dat$immport_vaccination_time %in% user.select$analysisGroups$group1] <- "A"
    pheno.dat$immport_vaccination_time_groups[pheno.dat$immport_vaccination_time %in% user.select$analysisGroups$group2] <- "B"
    pheno.dat$immport_vaccination_time_groups <- as.factor(pheno.dat$immport_vaccination_time_groups)
    
    synergy.terms <- lapply(user.select$studyCohort, function(s){
      paste(s, "immport_vaccination_time_groups", sep =':')
    })
    synergy.terms <- unlist(synergy.terms)
    
    de.formula <- as.formula(paste("~" , paste(c(user.select$studyCohort, "immport_vaccination_time_groups", synergy.terms), collapse = ' + ')))
  }
  
  fit <- limma::eBayes(limma::lmFit(expr.dat, model.matrix(de.formula, data = pheno.dat)))
  top.table <- limma::topTable(fit, 2000, coef=2)
  
  # -----------------------------------------------------
  # leave filtering to frontend 
  # write table as json file 
  # -----------------------------------------------------
  toJSON(as_tibble(top.table, rownames = "gene_name"))
}

# ----------------------------------
# Very simple RESTful API via plumber
# ----------------------------------

#* Echo back the input. This is more like a test code to check the service.
#* @param msg The message to echo
#* @get /echo
function(msg="") {
  list(msg = paste0("The message is: '", msg, "'"))
}

#* Perform the differential expression analysis.
#* @param selection The JSON selection from the web front-end
#* @post /doDiffExpAnalysis
function(selection.json) {
  if (is.na(selection.json)) {
    return(list(msg = "No selection.json is provided!"))
  }
  do_diff_exp_analysis(selection.json)
}
