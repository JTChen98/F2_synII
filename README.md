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

-   [align.txt](https://github.com/JTChen98/F2_synII/blob/main/Data/align.txt)

    The alignment result between the wild-type chomosome II reference sequence and synII.

**Note: We performed multiple regression for the phenotypes with the following covariates: effect of growth on control media YPD, 384-well measurement plate number, row, column, and round. The residuals were used for downstream analysis.**

## Analysis

### QTL mapping

QTL analysis was performed using the R package [qtl v1.66](https://rqtl.org/), the code can be found in [QTL.mapping.R](https://github.com/JTChen98/F2_synII/blob/main/Code/QTL.mapping.R).

-   Get chromosome offsets (to convert from marker index per chromosome to marker index across genome)

``` r
mindex.split = getMarkerIndexSplit(cross)
chr.mindex.offset = sapply(mindex.split, min)-1
```

-   Extract phenotypes, genotypes, and number of individuals phenotyped per trait

``` r
gdata = extractGenotype(cross)
n.pheno = countStrainsPerTrait(cross$pheno) 
pdata.01 = extractScaledPhenotype(cross, TRUE)
```

-   Calculate empirical FDRs using a permutation approach.

``` r
cross_new <- cross
cross_new$pheno <- as.data.frame(pdata.01)
operm <- scanone(cross_new, method = 'mr', n.perm=1000)
summary(operm, alpha=0.05)
```

-   Calculate LOD scores for each genotypic marker and each trait

``` r
LODS.01 = get.LOD.by.COR(n.pheno, pdata.01, gdata)
LODS.01s = LODmatrix.2.scanone(LODS.01, cross)
peaklist.01 = getChrPeaks(mindex.split, chr.mindex.offset, LODS.01) 
peakArray.01 = getPeakArray(peaklist.01, 3.63)
```

-   Estimate the phenotypic residuals by fitting trait-specific linear models that included the significant QTL genotypes as additive covariates

``` r
pdata.02 = getPhenoResids(pdata.01, gdata, peakArray.01) 
```

**Note: This process of peak detection, calculation of empirical significance thresholds and expansion of the linear model for each trait to include significant QTLs detected at each step was repeated four times.**

### Estimating heritability

Heritability analysis was performed using the R package [hglm v2.2.1](https://cran.r-project.org/web/packages/hglm/), the code can be found in [h2.estimate.R](https://github.com/JTChen98/F2_synII/blob/main/Code/h2.estimate.R).

``` r
geno[1:5,1:5]
```

    ##                 marker A01_01 A01_02 A01_03 A01_04
    ## 27915_chr01_27915_T_C      R      B      R      R
    ## 28323_chr01_28323_G_A      R      B      R      R
    ## 28652_chr01_28652_G_T      R      B      R      R
    ## 29667_chr01_29667_C_A      R      B      R      R
    ## 30756_chr01_30756_C_G      R      B      R      R

-   The design matrix for the fixed effects

``` r
X <- t(as.matrix(geno[,-1] == 'R'))*1
colnames(X) <- geno$marker
X <- scale(X)
A <- tcrossprod(X)/ncol(X)
X <- matrix(1, nrow(A), 1)
```

-   The design matrix for the random effects

``` r
svdA <- svd(A)
Z <- svdA$u %*% diag(sqrt(svdA$d))
```

-   Fitting Hierarchical Generalized Linear Model

``` r
model <- hglm(y = phenotype, X = X, Z = Z)
```

### Gene expression quantification

Gene expressions from the RNA-sequencing experiment were quantified using [XAEM v0.1.2](https://github.com/WenjiangDeng/XAEM/), the code can be found in [XAEM.sh](https://github.com/JTChen98/F2_synII/blob/main/Code/XAEM.sh).

-   Indexing transcripts

``` sh
TxIndexer -t /path/to/transcripts.fa -o /path/to/TxIndexer_idx
```

-   Construction of the X matrix

``` sh
Rscript XAEM_home/R/genPolyesterSimulation.R /path/to/transcripts.fa /path/to/design_matrix
GenTC -i /path/to/TxIndexer_idx -l IU -1 /path/to/design_matrix/sample_01_1.fasta -2 /path/to/design_matrix/sample_01_2.fasta -p 8 -o /path/to/design_matrix
Rscript XAEM_home/R/buildCRP.R in=/path/to/design_matrix/eqClass.txt out=/path/to/design_matrix/X_matrix.RData H=0.025
```

-   Construction of the Y count matrix

``` sh
Rscript Create_count_matrix.R workdir=/path/to/XAEM_project core=8
```

-   Updating the X matrix and isoform expression using AEM algorithm

``` sh
Rscript AEM_update_X_beta.R workdir=/path/to/XAEM_project core=8 design.matrix=X_matrix.RData isoform.out=XAEM_isoform_expression.RData paralog.out=XAEM_paralog_expression.RData merge.paralogs=FALSE isoform.method=average remove.ycount=TRUE
```

## Citation

If you want to cite this paper:

-   Chen, J., Zheng, J., Wang, Y., Wang, Y., Feng, X., Li, T., Chen, X., Liu, L., Fu, X., Shen, Y., Shen, X. The nature of supernature: Unveiling the doppelgänger effect of a synthetic chromosome. (2024)  
