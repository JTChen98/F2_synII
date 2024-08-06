## Data
-   [BYxRM_GenoData.txt](https://github.com/JTChen98/F2_synII/blob/main/Data/BYxRM_GenoData.txt)

    The genotype data for 1,008 F2 intercross segregants.

-   [pheno_raw_BY.Rdata](https://github.com/JTChen98/F2_synII/blob/main/Data/pheno_raw_BY.Rdata)
    
    Contains an R list object called 'pheno_raw_BY'. Each element in the list is a different phenotype. This file contains replicate the measurements for each phenotype and segregant in BY group, where each of the F2 haploid individuals was mated with lab strains derived from BY4741.
    
-   [pheno_raw_synII.Rdata](https://github.com/JTChen98/F2_synII/blob/main/Data/pheno_raw_synII.Rdata)
    
    Contains an R list object called 'pheno_raw_synII'. Each element in the list is a different phenotype. This file contains replicate the measurements for each phenotype and segregant in synII group, where each of the F2 haploid individuals was mated with the same lab strains but whose chromosome II had been replaced by the synthetic version.

-   [cross_BY_res.Rdata](https://github.com/JTChen98/F2_synII/blob/main/Data/cross_BY_res.Rdata)
   
    Contains an R/QTL cross object called 'cross_BY_res', with average phenotype residuals for each segregant in BY group and genotypes calls for each segregant. Also contains the genetic map.

-   [cross_synII_res.Rdata](https://github.com/JTChen98/F2_synII/blob/main/Data/cross_synII_res.Rdata)
   
    Contains an R/QTL cross object called 'cross_synII_res', with average phenotype residuals for each segregant in synII group and genotypes calls for each segregant. Also contains the genetic map.

**Note: We performed multiple regression for the phenotypes with the following covariates: effect of growth on control media YPD, 384-well measurement plate number, row, column, and round. The residuals were used for downstream analysis.**

## Analysis

### QTL mapping

QTL analysis was performed using the R package [qtl v1.66](https://rqtl.org/), the code can be found in [QTL.mapping.R](https://github.com/JTChen98/F2_synII/blob/main/Code/QTL.mapping.R).

1.

2.

3.

### Estimating heritability

Heritability analysis was performed using the R package [hglm v2.2.1](https://cran.r-project.org/web/packages/hglm/), the code can be found in [h2.estimate.R](https://github.com/JTChen98/F2_synII/blob/main/Code/h2.estimate.R).

1.

2.

3.

### Gene expression quantification

Gene expressions from the RNA-sequencing experiment were quantified using [XAEM v0.1.2](https://github.com/WenjiangDeng/XAEM/), the code can be found in [XAEM.sh](https://github.com/JTChen98/F2_synII/blob/main/Code/XAEM.sh).

1.

2.

3.

### Estimating genetic effect

1.

2.

3.

### Estimating Epistasis effect

1.

2.

3.
