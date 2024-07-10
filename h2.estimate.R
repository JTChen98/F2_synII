require('hglm')

geno <- read.table('~/Desktop/BYxRM_GenoData.txt', header = TRUE)
load('~/Desktop/cross_synII_res.Rdata')

trait <- names(cross_synII_res$pheno)

X <- t(as.matrix(geno[,-1] == 'R'))*1
colnames(X) <- geno$marker
X <- scale(X)
A <- tcrossprod(X)/ncol(X)
svdA <- svd(A)
Z <- svdA$u %*% diag(sqrt(svdA$d))

pheno_synII_res <- cross_synII_res$pheno
c <- target <- c()
for (i in 1:length(trait)) {
  phe <- trait[i]
  naidx <- which(is.na(pheno_synII_res[,phe]))
  if(length(naidx) > 0){
    model <- hglm(y = pheno_synII_res[-naidx,phe], X = matrix(1, nrow(A) - length(naidx), 1), Z = Z[-naidx,])
  } else {
    model <- hglm(y = pheno_synII_res[,phe], X = matrix(1, nrow(A), 1), Z = Z)
  }
  target <- data.frame(h2 = model$varRanef/(model$varRanef + model$varFix), trait = phe)
  c <- rbind(c, target)
}
