require(qtl)

# location of gwasFx.R
source('~/Desktop/Yeast/Round/BY4741_synII/gwasFx.R)

# Loads R/QTL object 'cross' (contains average phenotype for each segregant for each trait,
# markers for QTL mapping and genetic map
load('~/Desktop/Yeast/Round/BY4741_synII/synII/data/cross_synII_res.Rdata')
load('~/Desktop/Yeast/Round/BY4741_synII/BY4741/data/cross_BY_res.Rdata')

# Loads genotype data
geno <- read.table('~/Desktop/Yeast/data/BYxRM_GenoData.txt', header = T)

trait <- names(cross_synII_res$pheno)

for(i in 1:7){
  yy_synII <- scale(cross_synII_res$pheno[[i]])
  naidx <- which(is.na(yy_synII))
  if(length(naidx) > 0){
    synII <- gwas(y = yy_synII[-naidx], X = xx[-naidx,])
    }else{synII <- gwas(y = yy_synII, X = xx)}
  yy_BY <- scale(cross_BY_res$pheno[[i]])
  naidx <- which(is.na(yy_BY))
  if(length(naidx) > 0){
    BY4741 <- gwas(y = yy_BY[-naidx], X = xx[-naidx,])
    }else{BY4741 <- gwas(y = yy_BY, X = xx)}
}
