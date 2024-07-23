require(qtl)

# Load the qtl list mapped by QTL.mapping.R
qtl_list <- read.csv('~/Desktop/qtl.csv', header = T, stringsAsFactors = F)

# Load the genotype of segregants and object 'cross' (contains average phenotype for each segregant for each trait).
geno <- read.table('~/Desktop/Yeast/Round/1000seq/data/BYxRM_GenoData.txt', header = T, stringsAsFactors = F)
load('~/Desktop/Yeast/Round/1000seq/data/cross_synII_res.Rdata')

# Estimate the phenotypic variance explained by detected qtl for each trait
for(i in 1:7){
  sug1 <- cross_synII_res
  sug1$pheno <- data.frame(cross_synII_res$pheno[[i]])
  out.em1 <- scanone(sug1, method = 'mr')
  for (j in 1:nrow(out.em1)) {
    out.em1$POS[j] <- strsplit(rownames(out.em1)[j], '[_]')[[1]][3]
    out.em1$marker[j] <- rownames(out.em1)[j]
  }
  chr <- qtl_list[qtl_list$trait == trait[i],]$chr
  POS <- qtl_list[qtl_list$trait == trait[i],]$pos
  qtl <- as.data.frame(cbind(chr, POS))
  qtl <- merge(qtl, out.em1, by = c('chr', 'POS'))
  qtl <- as.data.frame(qtl$marker)
  colnames(qtl) <- 'marker'
  geno_new <- merge(qtl, geno, by = 'marker')
  yy <- sug1$pheno
  yy <- yy[,1]
  xx  = t(geno_new[,-1])
  xx = xx =='R'
  xx = xx*1
  target <- data.frame(r2 = summary(lm(yy~xx))$r.squared, 
                       se = sqrt(((summary(lm(yy~xx))$r.squared)^2)/qchisq(pf(summary(lm(yy~xx))$fstatistic[[1]],summary(lm(yy~xx))$fstatistic[[2]],summary(lm(yy~xx))$fstatistic[[3]],
                       lower.tail=F), 1, lower.tail=F)))
  c <- rbind(c, target)
}
