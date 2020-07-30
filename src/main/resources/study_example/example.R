usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("pacman")


# ----------------------------------------------------
p_load(dplyr)
p_load(tidyr)
p_load(purrr)
p_load(tidyverse)
p_load(ggplot)
p_load(limma)
p_load(plotly)


library(limma)


set.seed(1234)
options(stringsAsFactors = F)
# ----------------------------------------------------
sample.inf <- read_csv("biosample_metadata.csv")
# ----------------------------------------------------
# USECASE - One vaccine (Fluvirin) Age difference DE genes 
# with 0 as pre-vaccination - ie. case and control cohorts
# ----------------------------------------------------
this.path <- "clean/"
studies <- c("SDY400", "SDY404", "SDY520")
sample.inf.subset <- sample.inf[which(sample.inf$immport_study_accession %in% studies), ]

# read in pdata and fetch age groupings -------------
clean.dat.dirs <- list.dirs(path = this.path)
clean.dat.dirs <- clean.dat.dirs[which(!clean.dat.dirs %in% this.path)]
clean.dat.dirs2 <- clean.dat.dirs[grepl(paste(unique(sample.inf.subset$gse), collapse = "|") , clean.dat.dirs)]

all.age.groups <- lapply(clean.dat.dirs2, function(this.dir){
  dir <- this.dir
  files <- list.files(dir) 
  file.path <- paste(dir, files[grepl("pdata", files)], sep = "/")
  pdat <- read_csv(file.path)
  age.group <- tolower(pdat$`age group:ch1`)
  if("frail" %in% age.group) {
    age.group[which(age.group %in% "frail")] <- "old" 
  }
  gsm <- pdat$X1
  
  return(tibble(gsm=gsm, age_group=age.group))
})
names(all.age.groups) <- gsub(this.path, "", clean.dat.dirs2)
all.age.groups <- do.call(rbind, all.age.groups)
sum(sample.inf.subset$gsm %in% all.age.groups$gsm)
sample.inf.subset <- merge(x = sample.inf.subset, y = all.age.groups, by = "gsm", all.x = TRUE)

# ------------------------------------------------------------------------
all.genes <- lapply(clean.dat.dirs, function(this.dir){
  dir <- this.dir
  files <- list.files(dir) 
  file.path <- paste(dir, files[grepl("final_exprs", files)], sep = "/")
  expr <- read_csv(file.path)
  genes <- expr$X1
  
  return(genes)
})
common.genes <- Reduce(intersect, all.genes)
length(common.genes) 

all.expr <- lapply(clean.dat.dirs, function(this.dir){
  dir <- this.dir
  files <- list.files(dir) 
  file.path <- paste(dir, files[grepl("final_exprs", files)], sep = "/")
  expr <- read_csv(file.path)
  expr <- expr[which(expr$X1 %in% common.genes), ]
  expr <- expr[match(common.genes, expr$X1), ]
  expr <- expr[, -1]
  expr <- as.data.frame(expr)
  return(expr)
})
all.expr.df <- do.call(cbind, all.expr)
rownames(all.expr.df) <- common.genes
dim(all.expr.df)

all.expr.df <- all.expr.df[, which(colnames(all.expr.df) %in% sample.inf$gsm)]
all.expr.df <- all.expr.df[,match(sample.inf$gsm, colnames(all.expr.df))]
# -----------------------------------------------------------------------
covars.formula <- as.formula(paste("~", paste(c("gpl", "gender", "race"), collapse =' + ')))

fit <- limma::eBayes(limma::lmFit(all.expr.df, model.matrix(covars.formula, data = sample.inf)))
all.expr.df.adjusted <- as.matrix(residuals(fit, all.expr.df))

pc.dat <- all.expr.df.adjusted[sample(rownames(all.expr.df.adjusted), 1000), ]

PC <- prcomp(pc.dat, scale.=T, center = T)

plotdata <- data.frame(SampleID=rownames(PC$rotation),
                       PC1=PC$rotation[,1],
                       PC2=PC$rotation[,2])


plotdata$gpl <- sample.inf$gpl
plotdata$study <- sample.inf$immport_study_accession
plotdata$vaccine <- sample.inf$immport_immune_exposure_materia_id


g1 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=gpl)) + geom_point()
g3 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=study)) + geom_point()
g2 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=vaccine)) + geom_point()
subplot(g1, g2, g3)
# -----------------------------------------------------------------------
dat.adjusted <- all.expr.df.adjusted[, which(colnames(all.expr.df.adjusted) %in% sample.inf.subset$gsm)]
dat.adjusted <- dat.adjusted[,match(sample.inf.subset$gsm, colnames(dat.adjusted))]

sample.inf.subset$age_group <- as.factor(sample.inf.subset$age_group)
# ------------------------------------------------------------------------
mdf.mad <- apply(dat, 1, mad)
mdf.mean <- apply(dat, 1, mean)
smoothScatter(x = mdf.mean, y = mdf.mad)

variable.genes <- names(mdf.mad[mdf.mad >  0.6])
length(variable.genes)

dat <- dat.adjusted[which(rownames(dat.adjusted) %in% variable.genes), ]

fit <- limma::eBayes(limma::lmFit(dat, model.matrix(~ age_group + immport_vaccination_time + age_group:immport_vaccination_time, data = sample.inf.subset)))
top.table <- as.data.frame(limma::topTable(fit, number=2000 ,coef=2))

sum(top.table$adj.P.Val <= 0.05)
sum(top.table$P.Value <= 0.05)
summary(abs(top.table$logFC))
sum(top.table$logFC >= 1)
plot(sort(top.table$logFC))
top.table[which(abs(top.table$logFC) >= 1), ]
write_csv(tibble(rownames(top.table)[which(abs(top.table$logFC) >= 1)]), "test.csv")
# ------------------------------------------------------------------------