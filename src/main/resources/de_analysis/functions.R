# -----------------------------------------------------
# Analysis viz-helper functions
# -----------------------------------------------------
get_pca <- function(exp, sliced.pheno){
  # has to be called before variable gene filtration 
  pc.dat <- exp[sample(rownames(exp), 1000), ]
  
  PC <- prcomp(pc.dat, scale.=T, center = T)
  
  plotdata <- data.frame(SampleID=rownames(PC$rotation),
                         PC1=PC$rotation[,1],
                         PC2=PC$rotation[,2])
  
  plotdata$gpl <- sliced.pheno$gpl
  plotdata$study <- sliced.pheno$immport_study_accession
  plotdata$vaccine <- sliced.pheno$vo_vaccine_name
  plotdata$gender <- sliced.pheno$gender
  
  g1 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=gpl)) + geom_point()
  g2 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=vaccine)) + geom_point()
  g3 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=study)) + geom_point()
  g4 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=gender)) + geom_point()
  
  subplot(g1, g2, g3, g4)
}

get_heatmap <- function(pheno.dat, expr.dat, top.table, study, rank.range, ...){
  annotation_col = data.frame(
    gender = factor(pheno.dat[[study]])
  )
  rownames(annotation_col) = pheno.dat$gsm
  
  sig.genes <- top.table[rank.range, ]
  mat <- expr.dat[which(rownames(expr.dat) %in% rownames(sig.genes)), ]
  
  pheatmap(mat, annotation_col = annotation_col)
}