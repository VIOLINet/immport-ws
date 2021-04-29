#!/usr/bin/env Rscript
# -----------------------------------------------------
# Read absolute path inputs for expression, phenotype, and user-selections files
# Specify output absolute path and file name as the last arg to have generated results to be saved in
# -----------------------------------------------------
args = commandArgs(trailingOnly=T)

if(length(args) < 4) {
  print("Missing one of the absolute paths to generate analysis!")
} 

pheno.path <- args[1]
exp.path <- args[2]
selections.path <- args[3]
output.path <- args[4]
# -----------------------------------------------------
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("pacman")
# -----------------------------------------------------
# load libraries 
# -----------------------------------------------------
p_load("jsonlite")
p_load("tidyverse")
p_load("dplyr")
p_load("doParallel")
p_load("limma")

# -----------------------------------------------------
set.seed(1234)
options(stringsAsFactors = F)
# registerDoParallel(detectCores() - 2) # optional for back-end speed 
# -----------------------------------------------------
# read in data for pre-processed biosample expression and associated phenotypes 
# -----------------------------------------------------
all.pheno.dat <- read.csv(pheno.path) %>% as.data.frame()
all.expr.dat <- read_csv(exp.path) %>% as.data.frame()
rownames(all.expr.dat) <- all.expr.dat[, 1]
all.expr.dat <- all.expr.dat[, -1]

user.select <- read_json(selections.path, simplifyVector = TRUE)

# -----------------------------------------------------
# Filter biosamples by users study design selections &
# make sure biosamples exist & have matching orders 
# -----------------------------------------------------
pheno.dat <- all.pheno.dat[which(all.pheno.dat$gsm %in% user.select$GSMids), ]

expr.dat <- all.expr.dat[ ,which(colnames(all.expr.dat) %in% pheno.dat$gsm)]
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
# Note: changed to fetch upper quantile variable genes based on mad of expressions
# -----------------------------------------------------
if (user.select$variableGenes == T){
  mdf.mad <- apply(expr.dat, 1, mad)
  variable.genes <- names(mdf.mad[mdf.mad >  summary(mdf.mad)[[5]]])
  
  expr.dat <- expr.dat[which(rownames(expr.dat) %in% variable.genes), ]
}
# -----------------------------------------------------
# conduct differential expression analysis
# add time interactions 
# -----------------------------------------------------
if (user.select$modelTime == F){
  # ex. ~ age_group + immport_vaccination_time + age_group:immport_vaccination_time
  synergy.terms <- lapply(user.select$studyCohort, function(s){
    paste(s, "immport_vaccination_time", sep =':')
  })
  synergy.terms <- unlist(synergy.terms)
  pheno.dat$immport_subject_accession <- as.factor(pheno.dat$immport_subject_accession)
  
  de.formula <- as.formula(paste("~" , paste(c("immport_vaccination_time", "immport_subject_accession", user.select$studyCohort, synergy.terms), collapse = ' + ')))
  
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
  
  de.formula <- as.formula(paste("~" , paste(c("immport_vaccination_time_groups", "immport_subject_accession", user.select$studyCohort, synergy.terms), collapse = ' + ')))
}

fit <- limma::eBayes(limma::lmFit(expr.dat, model.matrix(de.formula, data = pheno.dat)))
top.table <- limma::topTable(fit, 30000, coef=2)

# -----------------------------------------------------
# leave filtering to frontend 
# write table as json file 
# -----------------------------------------------------
write_json(as_tibble(top.table, rownames = "gene_name"), output.path)
