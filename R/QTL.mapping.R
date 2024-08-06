library(lme4)
library(foreach)
library(doMC)
library(qtl)


registerDoMC(cores=8)

# location of QTL_mappingFx.R
source('~/Desktop/Yeast/Round/1000seq/reference/R_code/QTL_mappingFx.R')

# Loads 'pheno_raw' list (contains phenotype measurements for each trait)
#load('~/Desktop/Yeast/Round/1000seq/data/pheno_raw.Rdata')

# Loads R/QTL object 'cross' (contains average phenotype for each segregant for each trait,
# markers for QTL mapping and genetic map
load('~/Desktop/Yeast/Round/BY4741_synII/BY4741/data/cross_BY_res.Rdata')
cross <- cross_BY4741_res

# QTL Mapping ################################################################################################# 

mindex.split = getMarkerIndexSplit(cross)
# get chromosome offsets  (to convert from marker index per chromosome to marker index across genome)
chr.mindex.offset = sapply(mindex.split, min)-1

######extract phenotypes, genotypes, and number of individuals phenotyped per trait 
gdata     = extractGenotype(cross)
n.pheno   = countStrainsPerTrait(cross$pheno) 
pdata.01      = extractScaledPhenotype(cross, TRUE)
cross_new <- cross
cross_new$pheno <- as.data.frame(pdata.01)
operm <- scanone(cross_new, method = 'mr', n.perm=1000)
summary(operm, alpha=0.05)

LODS.01       = get.LOD.by.COR(n.pheno, pdata.01, gdata)
LODS.01s      = LODmatrix.2.scanone(LODS.01, cross)
peaklist.01   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.01) 
#LODS.01.FDR   = getPeakFDR(peaklist.01$chr.peaks.lod, pdata.01, gdata, 1000)
peakArray.01  = getPeakArray(peaklist.01, 3.63)

pdata.02      = getPhenoResids(pdata.01, gdata, peakArray.01) 
cross_new <- cross
cross_new$pheno <- as.data.frame(pdata.02)
operm <- scanone(cross_new, method = 'mr', n.perm=1000)
summary(operm, alpha=0.05)

LODS.02       = get.LOD.by.COR(n.pheno, pdata.02, gdata)
LODS.02s      = LODmatrix.2.scanone(LODS.02, cross, LODS.01s)
peaklist.02   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.02) 
#LODS.02.FDR   = getPeakFDR(peaklist.02$chr.peaks.lod, pdata.02, gdata, 1000)
peakArray.02  = getPeakArray(peaklist.02, 3.64)
#peakArray.02  = rbind(peakArray.01, peakArray.02)

pdata.03      = getPhenoResids(pdata.02, gdata, peakArray.02)
cross_new <- cross
cross_new$pheno <- as.data.frame(pdata.03)
operm <- scanone(cross_new, method = 'mr', n.perm=1000)
summary(operm, alpha=0.05)

LODS.03       = get.LOD.by.COR(n.pheno, pdata.03, gdata)
LODS.03s      = LODmatrix.2.scanone(LODS.03, cross, LODS.01s)
peaklist.03   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.03) 
#LODS.03.FDR   = getPeakFDR(peaklist.03$chr.peaks.lod, pdata.03, gdata, 1000)
peakArray.03  = getPeakArray(peaklist.03, 3.48)
#peakArray.03  = rbind(peakArray.02, peakArray.03)

pdata.04      = getPhenoResids(pdata.03, gdata, peakArray.03)
cross_new <- cross
cross_new$pheno <- as.data.frame(pdata.04)
operm <- scanone(cross_new, method = 'mr', n.perm=1000)
summary(operm, alpha=0.05)

LODS.04       = get.LOD.by.COR(n.pheno, pdata.04, gdata)
LODS.04s      = LODmatrix.2.scanone(LODS.04, cross, LODS.01s)
peaklist.04   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.04) 
#LODS.03.FDR   = getPeakFDR(peaklist.03$chr.peaks.lod, pdata.03, gdata, 1000)
peakArray.04  = getPeakArray(peaklist.04, 3.64)
#peakArray.03  = rbind(peakArray.02, peakArray.03)


pdata.05      = getPhenoResids(pdata.04, gdata, peakArray.04)
cross_new <- cross
cross_new$pheno <- as.data.frame(pdata.05)
operm <- scanone(cross_new, method = 'mr', n.perm=1000)
summary(operm, alpha=0.05)

LODS.05       = get.LOD.by.COR(n.pheno, pdata.05, gdata)
LODS.05s      = LODmatrix.2.scanone(LODS.05, cross, LODS.01s)
peaklist.05   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.05) 
#LODS.03.FDR   = getPeakFDR(peaklist.03$chr.peaks.lod, pdata.03, gdata, 1000)
peakArray.05  = getPeakArray(peaklist.05, 3.61)
#peakArray.03  = rbind(peakArray.02, peakArray.03)


pA1=cbind(peakArray.01, 'J1')
pA2=cbind(peakArray.02, 'J2')
pA3=cbind(peakArray.03, 'J3')
pA4=cbind(peakArray.04, 'J4')

names(pA1)[3]='jump'
names(pA2)[3]='jump'
names(pA3)[3]='jump'
names(pA4)[3]='jump'

peak.index = data.frame(rbind(pA1,pA2,pA3,pA4))
peak.index = peak.index[order(peak.index$trait, peak.index$markerIndex),]
peak.index = split(peak.index, peak.index$trait)

#save(peak.index, file='~/1000BYxRM/QTL/peakindex.bin')
# also ran with refine = TRUE 
#save(fQTLs_FDR05r, file = '~/1000BYxRM/QTL/fQTLs_FDR05r.bin')

fQTLs_FDR05 = foreach( i=1:length(peak.index) ) %dopar% doQTLModel(i,peak.index,cross, LODS.01s, refine=TRUE)
names(fQTLs_FDR05) = names(peak.index)
fQTLs_FDR05 =lapply(fQTLs_FDR05, function(x) {
        x$CIs = data.frame((as.vector(x$fqtl$ests$ests[-1])), x$CIs)
        names(x$CIs)[1]='eff.size'
        return(x)
        })
