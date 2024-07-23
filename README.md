The genotype data for 1,008 F2 intercross segregants can be found in BYxRM_GenoData.txt.


The raw phenotype data for each trait can be found in pheno_raw_BY.Rdata or pheno_raw_synII.Rdata.


We performed multiple regression for the phenotype for each trait with the following covariates: effect of growth on control media YPD, 384-well measurement plate number, row, column, and round. The residuals can be found in cross_BY_res.Rdata or cross_synII_res.Rdata.


Gene expressions from the RNA-sequencing experiment were quantified using XAEM v0.1.2 (https://github.com/WenjiangDeng/XAEM/), the code can be found in XAEM.sh.


Heritability analysis was performed using the R package hglm v2.2.1 (https://cran.r-project.org/web/packages/hglm/), the code can be found in h2.estimate.R.


QTL analysis was performed using the R package qtl v1.66 (https://rqtl.org/), the code can be found in QTL.mapping.R.
