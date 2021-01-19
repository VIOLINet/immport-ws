#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)

# Two parameters are needed: port and the data directory
if(length(args) < 2) {
  print("No parameters are assigned for the script. Use the default parameters.")
  args <- c(8087, "../data")
}

# The bootstrap code to start the RESTful API in R for differential expression analysis
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

usePackage("plumber")

pr("differential_exp_analysis.R")  %>% pr_run(port=as.integer(args[1]))