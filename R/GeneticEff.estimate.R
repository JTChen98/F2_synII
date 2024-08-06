# location of gwasFx.R
source('~/Desktop/Yeast/Round/BY4741_synII/gwasFx.R)

# Loads R/QTL object 'cross' (contains average phenotype for each segregant for each trait,
# markers for QTL mapping and genetic map
load('~/Desktop/Yeast/Round/BY4741_synII/synII/data/cross_synII_res.Rdata')
load('~/Desktop/Yeast/Round/BY4741_synII/BY4741/data/cross_BY_res.Rdata')

# Loads genotype data
geno <- read.table('~/Desktop/Yeast/Round/1000seq/data/BYxRM_GenoData.txt', header = T)
