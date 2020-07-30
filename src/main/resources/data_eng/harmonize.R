usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("pacman")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# install on first run -----------
# BiocManager::install("GEOquery")
# BiocManager::install("GEOmetadb")
library(GEOquery)
library(GEOmetadb)
# ----------------------------------------------------
p_load(jsonlite)
p_load(dplyr)
p_load(tidyr)
p_load(lubridate)
p_load(purrr)
p_load(tidyverse)
p_load(ggplot)


set.seed(1234)
options(stringsAsFactors = F)
# ---------------------------------------------------------------
# NOTE: This db doesn't uptodate frequently and may be missing GSM/GSE/GPL mappings.
# However, it's a faster resource for the samples that have information to fetch.
if ( !file.exists( "GEOmetadb.sqlite" ) ) {
  geometadbfile <- getSQLiteFile() # download full DB
} else {
  geometadbfile <- "GEOmetadb.sqlite"
}
# ---------------------------------------------------------------
# TODO: these files can be replaced with an immport sql query result arg param.
# they may have been fetched and assembled from different resources 
# here we use pdata.immport since it's the closest to immport data structure 
query.result <- read_json("query_studyid.json", simplifyVector = TRUE)
query.result <- query.result[grepl("GSM", query.result$gsm), ]
query.result <- query.result[!duplicated(query.result[c("gsm", "exposureMaterialId", "studyAccession")]), ]
# segment time into int and units if not seperated. 
if ("time" %in% names(query.result)){
  query.result <- query.result %>% separate(time, c("time_num", "time_unit"), " ")
  query.result$time_num <- as.integer(query.result$time_num)
}


pdata.immport <- readxl::read_excel("Immport_flu_vaccine_VO_2_geo_v3.xlsx")
# ---------------------------------------------------------------
# EDA: Counts and distribution of unique data features and indices 
summary(pdata.immport)

n_distinct(pdata.immport$immport_geo_accession) 
n_distinct(pdata.immport$immport_immune_exposure_materia_id) 
n_distinct(pdata.immport$immport_study_accession)

unique(query.result$time_unit) # "Hours" "Days" 
# ---------------------------------------------------------------
# find all the GEO entities of one type associated with another GEO entity. gsm -> (gpl and gse)
meta.ids <- geoConvert(pdata.immport$immport_geo_accession, out_type = c("gpl", "gse"))
meta.ids <- merge(x = meta.ids$gpl, y = meta.ids$gse, by = "from_acc", all.x = TRUE)
names(meta.ids) <- c("gsm", "gpl", "gse")

# EDA: Counts
n_distinct(meta.ids$gsm) 
n_distinct(meta.ids$gpl)
n_distinct(meta.ids$gse) 

sum(is.na(meta.ids$gsm))
sum(is.na(meta.ids$gpl))
sum(is.na(meta.ids$gse)) 
# ---------------------------------------------------------------
# Fetch missing GSM's inf from geo and EDA: Counts
missing.gse <- lapply(seq_along(meta.ids$gsm[is.na(meta.ids$gse)]), function(i){
  
  tryCatch({temp <- getGEO(meta.ids$gsm[is.na(meta.ids$gse)][i])
  return(temp@header$series_id)}, 
  error = function(e) paste(meta.ids$gsm[is.na(meta.ids$gse)][i], sep = "***"))
})
missing.gse.df <- tibble( gsm = meta.ids$gsm[is.na(meta.ids$gse)], gse = unlist(missing.gse))
meta.ids$gse[is.na(meta.ids$gse)] <- missing.gse.df$gse

dim(meta.ids)
dim(meta.ids[!duplicated(meta.ids[c("gsm", "gpl", "gse")]), ])
n_distinct(meta.ids$gsm) 
n_distinct(meta.ids$gpl) 
n_distinct(meta.ids$gse) 
# ---------------------------------------------------------------
# Fetch platform information 
get.platform.info <- function(platform){
  info <- getGEO(platform)
  dat <- tibble(gpl = platform, organism = info@header$organism, platform_name = info@header$manufacturer, platform_desc = info@header$title)
  return(dat)
}

platform <- meta.ids$gpl %>% unique()
platform.list <- lapply(platform, function(p){
  return(get.platform.info(p))
})
platform.info <- do.call(rbind, platform.list)


meta.ids <- merge(x = meta.ids, y = unique(platform.info[c("gpl", "platform_name", "platform_desc")]), by = "gpl", all.x = TRUE)
meta.ids$platform_name <- str_replace(meta.ids$platform_name, "Inc.", "") %>% str_replace(., ",", "") %>% str_squish()
dim(meta.ids) 
# ---------------------------------------------------------------
# TODO: change 12 with standard feature name of the column in future. 
names(pdata.immport)[12] <- "gsm"
final.meta <- merge(x = meta.ids, y = pdata.immport, by = "gsm", all.x = TRUE)
dim(final.meta)  
# write.csv(final.meta, "final_meta_edison.csv")
# ---------------------------------------------------------------
# Summary df/matrix of immport studies to aid data mapping and integration 
# for downstream analysis 
# ---------------------------------------------------------------
summary.df.study <- final.meta %>% 
  group_by(immport_study_accession) %>% 
  summarise(vaccine = paste(unique(immport_immune_exposure_materia_id), collapse=", "), 
            platform_count = n_distinct(platform_name), 
            platform_name = paste(unique(platform_name), collapse=", "), 
            platform_desc = paste(unique(platform_desc), collapse=", "),
            gse_count = n_distinct(gse),
            gpl_count = n_distinct(gpl),
            gsm_count = n_distinct(gsm),
            gpl = paste(unique(gpl), collapse=", "), 
            gse = paste(unique(gse), collapse=", ")) 

sum(summary.df.study$gsm_count)
# ---------------------------------------------------------------
# Fetch upstream normalization method for each study get GEO
# ---------------------------------------------------------------
norm.list <- lapply(seq_along(summary.df.study$immport_study_accession), function(s){
  gsm.vec <- unlist(strsplit(summary.df.study[s, "gse"][[1]][1], split=", "))
  gsm.vec <- gsm.vec[which(!gsm.vec %in% "NA")]
  
  l <- c()
  norm.method  <- c()
  
  if (length(gsm.vec) > 0){
    
    for (i in 1:length(gsm.vec)) {
      gse <- getGEO(gsm.vec[i])
      for (j in 1:length(gse)){
        norm.method[j] <- unique(gse[[j]]@phenoData@data$data_processing)
        l[i] <- norm.method
      }
    }
    
  }
  return(paste(unique(l), collapse=", "))
})

summary.df.study$data_processing <- unlist(norm.list)
# write.csv(summary.df.study, "summary_metadata.csv")
# ------------------------------------------------------------
# Assess sample overlap between studies 
# ------------------------------------------------------------
final.meta.unique <- final.meta[!duplicated(final.meta$gsm), ] # duplicated doesn't overlap either
l <- list() 
v <- c()
m = 1
for (i in summary.df.study$immport_study_accession){
  print(i)
  n = 1
  v <- c()
  for (j in summary.df.study$immport_study_accession){
    v[n] <- sum(final.meta.unique$gsm[which(final.meta.unique$immport_study_accession %in% i)] %in% 
                  final.meta.unique$gsm[which(final.meta.unique$immport_study_accession %in% j)])
    n = n+1
  }
  print(v)
  l[[m]] <- v
  m = m+1
}
sample.overlap.count <- data.frame(do.call(cbind, l))
names(sample.overlap.count) <- summary.df.study$immport_study_accession
rownames(sample.overlap.count) <- summary.df.study$immport_study_accession
# write.csv(sample.overlap.count, "sample_overlap_count_study.csv")

# ------------------------------------------------------------
# given a study ID, assess GSE and GPL mappings, choose a one-one GSE for this study that matches 
# samples counts for this study and doesn't overlap with other studies, 
# then fetch GEO pheno and feature metadata along upstream normalized expr() probes X GSM 
# ------------------------------------------------------------
sid <- "SDY212"
sid.gse <- strsplit(summary.df.study$gse[summary.df.study$immport_study_accession %in% sid], ", ")[[1]]

length(sid.gse) # check how many GSE - assess which to choose from summary.df.study and sample overlaps in GSE

gse.level <- 1
sid.gse[gse.level]
sid.dat <- getGEO(sid.gse[gse.level])

length(sid.dat) # how many GPL 
sid.dat[[1]]@annotation # GPL

sid.fdata <- sid.dat[[1]]@featureData@data
sid.pdata <- sid.dat[[1]]@phenoData@data
sid.exprs <- sid.dat[[1]]@assayData$exprs

dim(sid.exprs)
sid.exprs <- sid.exprs[, which(colnames(sid.exprs) %in% pdata.immport$gsm)]
sid.pdata <- sid.pdata[which(rownames(sid.pdata) %in% pdata.immport$gsm), ]
dim(sid.exprs)

# ------------------------------------------------------------
# collapse expressions into unique gene X GSM df/matrix 
# workflow and data structs differ between illumina and affy 
# ------------------------------------------------------------
# TODO: if statement to handle the cases 
##########################################
# ----- illumina bead array case -----
##########################################
sid.fdata$ILMN_Gene[grep("[.]", sid.fdata$ILMN_Gene)] <- gsub("\\..*", "", sid.fdata$ILMN_Gene[grep("[.]", sid.fdata$ILMN_Gene)])
sid.dup.genes <- unique(sid.fdata$ILMN_Gene[duplicated(sid.fdata$ILMN_Gene)])

collapsed.probes.list <- lapply(sid.dup.genes, function(gene){
  print(gene)
  sid.probes <- rownames(sid.fdata)[sid.fdata$ILMN_Gene %in% gene]
  print(sid.probes)
  return(apply(sid.exprs[which(rownames(sid.exprs) %in% sid.probes), ], 2, median))
})
collapsed.probes <- data.frame(do.call(rbind, collapsed.probes.list))
rownames(collapsed.probes) <- sid.dup.genes


sid.unique.genes <- sid.fdata$ILMN_Gene[which(!sid.fdata$ILMN_Gene %in% sid.dup.genes)]
sid.unique.genes.probes <- sid.fdata[which(sid.fdata$ILMN_Gene %in% sid.unique.genes), c("ID", "ILMN_Gene")]
unique.probes <- sid.exprs[which(rownames(sid.exprs) %in% sid.unique.genes.probes$ID), ]


sid.unique.genes.probes <- sid.unique.genes.probes[match(rownames(unique.probes), sid.unique.genes.probes$ID), ]
rownames(unique.probes) <- sid.unique.genes.probes$ILMN_Gene

sid.final.exprs <- data.frame(rbind(collapsed.probes, unique.probes))

# check gene names, log transform, and expression dimention 
rownames(sid.final.exprs)[grep("[.]", rownames(sid.final.exprs))]
summary(sid.final.exprs)
summary(log2(sid.final.exprs))
dim(sid.final.exprs)

##########################################
# ----- Affymetrix bead array case -----
##########################################
sid.fdata$`Gene Symbol`[sid.fdata$`Gene Symbol`==""] <- NA

length(sid.fdata$`Gene Symbol`)
length(unlist(strsplit(sid.fdata$`Gene Symbol`, split = " /// ")))

gs <- strsplit(sid.fdata$`Gene Symbol`, split = " /// ")
gs.expanded <- data.frame(ID = rep(sid.fdata$ID, sapply(gs, length)), gene_symbol = unlist(gs))
sid.fdata2 <- merge(x = gs.expanded, y = sid.fdata, by = "ID", all.x = TRUE)


sid.fdata2$gene_symbol[grep("[.]", sid.fdata2$gene_symbol)] <- gsub("\\..*", "", sid.fdata2$gene_symbol[grep("[.]", sid.fdata2$gene_symbol)])

sid.dup.genes <- unique(sid.fdata2$gene_symbol[duplicated(sid.fdata2$gene_symbol)])
sid.dup.genes <- sid.dup.genes[!is.na(sid.dup.genes)]

collapsed.probes.list <- lapply(sid.dup.genes, function(gene){
  print(gene)
  sid.probes <- sid.fdata2$ID[sid.fdata2$gene_symbol %in% gene]
  print(sid.probes)
  
  if(is.null(dim(sid.exprs[which(rownames(sid.exprs) %in% sid.probes), ]))) {
    return(sid.exprs[which(rownames(sid.exprs) %in% sid.probes), ])
  }else{
    return(apply(sid.exprs[which(rownames(sid.exprs) %in% sid.probes), ], 2, median))
  }
  
})

collapsed.probes <- data.frame(do.call(rbind, collapsed.probes.list))
rownames(collapsed.probes) <- sid.dup.genes

sid.unique.genes <- sid.fdata2$gene_symbol[which(!sid.fdata2$gene_symbol %in% sid.dup.genes)]
sid.unique.genes <- sid.unique.genes[!is.na(sid.unique.genes)]
sid.unique.genes.probes <- unique(sid.fdata2[which(sid.fdata2$gene_symbol %in% sid.unique.genes), c("ID", "gene_symbol")])
unique.probes <- sid.exprs[which(rownames(sid.exprs) %in% sid.unique.genes.probes$ID), ]


sid.unique.genes.probes <- sid.unique.genes.probes[match(rownames(unique.probes), sid.unique.genes.probes$ID), ]
rownames(unique.probes) <- sid.unique.genes.probes$gene_symbol

sid.final.exprs <- data.frame(rbind(collapsed.probes, unique.probes))
rownames(sid.final.exprs)[duplicated(rownames(sid.final.exprs))]

# check gene names, log transform, and expression dimention 
summary(sid.final.exprs)
summary(log2(sid.final.exprs))
rownames(sid.final.exprs)[grep("[.]", rownames(sid.final.exprs))]
# --------------------------------------------------
# log transform if needed 
# --------------------------------------------------
sid.final.exprs.test <- log2(sid.final.exprs) 
summary(sid.final.exprs.test)
sid.final.exprs.test <- do.call(data.frame,lapply(sid.final.exprs.test, function(x) replace(x, is.infinite(x), 0)))
sid.final.exprs.test <- do.call(data.frame,lapply(sid.final.exprs.test, function(x) replace(x, is.na(x), 0)))

summary(sid.final.exprs.test)
sid.final.exprs <- sid.final.exprs.test
# --------------------------------------------------
# EDA: look into final expression matrix data structure 
# PCA, mean, median, mad 
# NOTE: sample outliers may be removed via PCA 
# here the final PCA merged via choosen GSE data didn't show any sample outliers
df <- sid.final.exprs[sample(rownames(sid.final.exprs), 1000), ]

PC <- prcomp(df, scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(SampleID=rownames(PC$rotation),
                       PC1=PC$rotation[,1],
                       PC2=PC$rotation[,2])


ggplot(plotdata, aes(x=PC1, y=PC2)) + geom_point()

mdf.mad <- apply(df, 1, mad)
mdf.mean <- apply(df, 1, mean)
smoothScatter(x = mdf.mean, y = mdf.mad)
plot(sort(mdf.mean))
plot(sort(mdf.mad))
# ----------------------------------------------------
# save results for downstream merging and integration 
# ----------------------------------------------------
dir.create(paste0("all_data/", sid.gse[gse.level]))


write.csv(sid.final.exprs, paste0(paste0(paste0("all_data/", sid.gse[gse.level]), "/"), "final_exprs.csv"))
write.csv(sid.fdata, paste0(paste0(paste0("all_data/", sid.gse[gse.level]), "/"), "fdata.csv"))
write.csv(sid.pdata, paste0(paste0(paste0("all_data/", sid.gse[gse.level]), "/"), "pdata.csv"))
write.csv(sid.exprs, paste0(paste0(paste0("all_data/", sid.gse[gse.level]), "/"), "exprs.csv"))


pdf(file = paste0(paste0("all_data/",   paste0(sid.gse[gse.level] ,"/")), paste0(sid.gse[gse.level], ".pdf")),  wi = 12, he = 9)
ggplot(plotdata, aes(x=PC1, y=PC2)) + geom_point()
plot(sort(mdf.mean))
plot(sort(mdf.mad))
smoothScatter(x = mdf.mean, y = mdf.mad)
dev.off()
# ----------------------------------------------------
