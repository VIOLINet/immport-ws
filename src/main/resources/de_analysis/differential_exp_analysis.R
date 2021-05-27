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
pheno.path <- paste(dir, "biosample_metadata.csv", sep = "/")
all.exp.path <- paste(dir, "gene_expression_matrix.csv", sep = "/")

# ------ Load data and then cached them -----------------------------------------------------
set.seed(1234)
options(stringsAsFactors = F)
# registerDoParallel(detectCores() - 2) # optional for back-end speed 
# -----------------------------------------------------
# read in data for pre-processed biosample expression and associated phenotypes 
# -----------------------------------------------------
all.pheno.dat <- read.csv(pheno.path) %>% as.data.frame()
# # We will do dynamic loading for the expression data
all.exp.data <- NULL

load_exp_data <- function(file.name) {
  all.expr.dat <- read_csv(all.exp.path) %>% as.data.frame()
  rownames(all.expr.dat) <- all.expr.dat[, 1]
  all.expr.dat <- all.expr.dat[, -1]
}

# Make sure there are at least two values for the selected variable
validate_variables <- function(variables, pheno.dat) {
  validate.variables <- c() # To be returned
  for (var in variables) {
    values <- pheno.dat[, var]
    unique.values <- unique(values)
    if ((length(unique.values)) > 1) {
      validate.variables <- c(validate.variables, var)
    }
  }
  validate.variables
}

# Perform the differential expression analysis via limma
do_diff_exp_analysis <- function(selection.json.text) {
  user.select <- fromJSON(selection.json.text, simplifyVector = TRUE)
  # -----------------------------------------------------
  # Filter biosamples by users study design selections &
  # make sure biosamples exist & have matching orders 
  # -----------------------------------------------------
  pheno.dat <- NA
  if (!is.null(user.select$GSMids)) {
    pheno.dat <- all.pheno.dat[which(all.pheno.dat$gsm %in% user.select$GSMids), ]
  }else {
    pheno.dat <- all.pheno.dat
    for (s in names(user.select$studyCohort)) {
      print(paste(s, user.select$studyCohort[s][[1]], sep = ": "))
      # print(pheno.dat[, s])
      selected <- which(pheno.dat[, s] %in% user.select$studyCohort[s][[1]])
      pheno.dat <- pheno.dat[selected, ]
    }
    # View(pheno.dat)
  }
  # View(pheno.dat)
  if (is.null(all.exp.data)) { # Use the expression data directly
    all.exp.data <<- load_exp_data(all.exp.path)
  }
  # Filter the expression to the user's selected samples
  expr.dat <- all.exp.data[ , which(colnames(all.exp.data) %in% pheno.dat$gsm)]
  # Make sure it has the same order as in the sample meta data.
  expr.dat <- expr.dat[ ,match(pheno.dat$gsm, colnames(expr.dat))]
  
  # For test
  # which <- rownames(expr.dat) %in% c("IGLL3")
  # selected.view <- expr.dat[which, ]
  # View(selected.view)
  # -----------------------------------------------------
  # make factors for groups selected 
  # -----------------------------------------------------
  # one of the study variables selection lengths has to be greater than 2 
  # treat groups as factors 
  # TODO: This converting may be removed. Need to check. Note by G.W.
  # -----------------------------------------------------
  if (length(user.select$studyVariables) == 1){
    pheno.dat[,which(colnames(pheno.dat) %in% user.select$studyVariables)] <- as.factor(pheno.dat[,which(colnames(pheno.dat) %in% user.select$studyVariables)])
  } else {
    for (s in user.select$studyVariables){
      pheno.dat[,which(colnames(pheno.dat) %in% s)] <- as.factor(pheno.dat[,which(colnames(pheno.dat) %in% s)])
    }
  }
  # -----------------------------------------------------
  # Filter variable genes for stronger analysis results 
  # based on users request 
  # Note: 0.6 seems to be a good threshold for this slice of data 
  # This should never be done when using limma:ebays. Note by G.W.
  # -----------------------------------------------------
  # if (user.select$variableGenes){
  #   mdf.mad <- apply(expr.dat, 1, mad)
  #   variable.genes <- names(mdf.mad[mdf.mad >  0.6])
  #   
  #   expr.dat <- expr.dat[which(rownames(expr.dat) %in% variable.genes), ]
  # }
  # -----------------------------------------------------
  # conduct differential expression analysis
  # -----------------------------------------------------
  
  # Collect all variables to be corrected
  # Other variables that should be corrected
  other.vars <- c()
  if (user.select$platformCorrection) {
    other.vars <- c(other.vars, "gpl") # add the platform as a new variable for correction
  }
  if (!is.null(user.select$usePairedData) && user.select$usePairedData) {
    other.vars <- c(other.vars, "immport_subject_accession") # paired analysis
  }
  
  # For test. It should be controlled by the user interface
  other.vars <- c(other.vars, "type_subtype")
  
  total.vars <- other.vars
  if (!is.null(user.select$studyVariables)) {
    total.vars <- c(total.vars, user.select$studyVariables)
  }
  # Make sure variable has at least two values
  total.vars <- validate_variables(total.vars, pheno.dat)
  
  coef.name <- NULL
  if (user.select$modelTime){
    # # ex. ~ age_group + immport_vaccination_time + age_group:immport_vaccination_time
    # Don't consider synergy.terms for the time being.
    # synergy.terms <- lapply(user.select$studyVariables, function(s){
    #   paste(s, "immport_vaccination_time", sep =':')
    # })
    # synergy.terms <- unlist(synergy.terms)
    coef.name <- "immport_vaccination_time"
  } else {
    # ex. ~ age_group + immport_vaccination_time_groups + age_group:immport_vaccination_time_groups
    pheno.dat$immport_vaccination_time_groups <- NA
    pheno.dat$immport_vaccination_time_groups[pheno.dat$immport_vaccination_time %in% user.select$analysisGroups$group1] <- "A"
    pheno.dat$immport_vaccination_time_groups[pheno.dat$immport_vaccination_time %in% user.select$analysisGroups$group2] <- "B"
    pheno.dat$immport_vaccination_time_groups <- as.factor(pheno.dat$immport_vaccination_time_groups)
    coef.name <- "immport_vaccination_time_groups"
    # # Don't consider interaction term for the time being. It is difficult to give then a biological explanation.
    # synergy.terms <- lapply(total.vars, function(s){
    #   paste(s, coef.name, sep =':')
    # })
    # synergy.terms <- unlist(synergy.terms)
    # total.vars <- c(total.vars, synergy.terms)
  }
  print(total.vars)
  # Make sure the coef.name is the first parametmer so that we can use coef = 2 in top.table
  de.formula <- as.formula(paste("~" , paste(c(coef.name, total.vars), collapse = ' + ')))
  print(de.formula)
  design <- model.matrix(de.formula, data = pheno.dat)
  # View(design)
  print(dim(design))
  print(dim(expr.dat))
  fit <- limma::lmFit(expr.dat, design)
  fit <- limma::eBayes(fit)
  # View(fit$coefficients)
  # View(fit$p.value)
  # Try to get all genes
  top.table <- limma::topTable(fit, coef=2, number = Inf)
  # View(top.table)
  # write.table(top.table, "top_table.csv", sep = ",", quote = F, col.names = T)
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
