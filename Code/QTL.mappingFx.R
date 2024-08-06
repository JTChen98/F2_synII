# functions for QTL Mapping

#### Split marker index by chromosome #################################################################
getMarkerIndexSplit = function(impcross) {
    mr = scanone(impcross, pheno.col=1, method='mr')
    marker.info       = parse.rownames(rownames(mr))
    marker.info$cmPOS = mr$pos
    mindex.split = split(rownames(marker.info), marker.info$NCHR)
    mindex.split = sapply(mindex.split, as.numeric)
    return(mindex.split)
}
#######################################################################################################

##### Extract genotype matrix from cross structure and recode as -1,1 #####################################
extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3 }
###########################################################################################################

##### Count number of strains with data for each phenotype from cross structure ###########################
countStrainsPerTrait = function(pheno) {apply(pheno, 2, function(x){sum(!is.na(x))})}
###########################################################################################################

##### Extract phenotype matrix from cross structure and mean center ... optional standardize Variance######
extractScaledPhenotype=function(impcross,scaleVar=FALSE){apply(impcross$pheno, 2, scale, scale=scaleVar)}
###########################################################################################################

###### calculate LOD score given matrices of phenotype and genotype data ##############################
# n.pheno is count of strains with data
# pheno is matrix (strains X phenotypes)
# gdata is matrix (strains X markers)
#.....  be careful about NAs and only using it for single marker scans 

get.LOD.by.COR = function(n.pheno, pheno, gdata) {
    # Lynch and Walsh p. 454
    return( (-n.pheno*log(1-cor(pheno, gdata, use='pairwise.complete.obs')^2))/(2*log(10)) ) }

#######################################################################################################

#######################################################################################################
getChrPeaks = function(mindex.split, chr.mindex.offset, LODS) {
    chr.peaks.lod    = sapply(mindex.split, function(markers) { apply(LODS[,markers], 1, max) })
    # get marker index of LOD peaks per chromosomes                             
    chr.peaks.index = sapply(mindex.split, function(markers)  { apply(LODS[,markers], 1, which.max) })
    # convert chromosome marker index to genome maker index                             
    chr.peaks.index = t(apply(chr.peaks.index, 1, function(x){x+chr.mindex.offset}))

    return(list(chr.peaks.lod = chr.peaks.lod, chr.peaks.index=chr.peaks.index))
}
#######################################################################################################

########################################################################################################
getPeakArray = function(peaklist, threshold) {
    tryCatch( {
    keepPeaks   = which(peaklist$chr.peaks.lod>threshold, arr.ind=T)
    kP = data.frame(rownames(keepPeaks), peaklist$chr.peaks.index[keepPeaks])
    names(kP)=c('trait', 'markerIndex') 
    kP = kP[order(kP$trait, kP$markerIndex),]
    return(kP)} ,error=function(e) {return(NULL) })
}
########################################################################################################

###### fix QTLs and get residuals phenotypes ###########################################################
getPhenoResids = function(pdata,gdata, peakArray, intercept=FALSE) {
    presids = pdata
    for( i in 1:ncol(pdata) ) {
        spA = peakArray[peakArray$trait==colnames(pdata)[i],]
        if(nrow(spA)>0){
            if(intercept) {
                  rr = residuals(lm(pdata[,i]~gdata[,spA$markerIndex]))
            }else{
                 rr = residuals(lm(pdata[,i]~gdata[,spA$markerIndex]-1))
            }
            presids[as.numeric(names(rr)),i]=rr
        }
    }
    return(presids)
} 
########################################################################################################

########################################################################################################################
doQTLModel = function(i,peak.index, fcross, LODSso, refine=TRUE){
    traitName = names(peak.index)[i]
    print(traitName)
    pheno.colIndex = match(traitName,names(fcross$pheno))

    qtlMarkers = unique(peak.index[[i]][,'markerIndex'])
   
    mqtl=NA
    fqtl=NA
    aqtl=NA
    rqtl=NA
    qtlMarkersCP = NA
    CIs =NA
    nQTL = length(qtlMarkers)
    if(nQTL>0) {
        qtlMarkersNames = rownames(LODSso)[qtlMarkers]
        qtlMarkersCP = LODSso[qtlMarkersNames ,c(1,2)]

        mqtl =  makeqtl(fcross, chr=qtlMarkersCP$chr, pos=qtlMarkersCP$pos, qtl.name=qtlMarkersNames, what='prob')
        if(refine == TRUE) {
            rqtl  = refineqtl(fcross,pheno.col=pheno.colIndex, mqtl, method='hk',
                           model='normal', incl.markers=TRUE,maxit=1)
        } else {
             rqtl = mqtl
        }
        fqtl  = fitqtl(fcross, pheno.col =pheno.colIndex, rqtl, method='hk', get.ests=TRUE, dropone=TRUE)
            
        #scan for interactors amongst QTLs with significant main effects
         if(nQTL==1) { aqtl=NA } 
         if(nQTL>1)  { aqtl = addint(fcross, pheno.col=pheno.colIndex, rqtl, method='hk', qtl.only=TRUE)      }

        # get confidence intervals
        if(refine ==TRUE) {
            CIs = lapply(1:length(qtlMarkers),function(x) {lodint(rqtl,qtl.index=x) })
            CIs = makeCI.df(CIs)
            }
    }
   return(list(qtlMarkers=qtlMarkers, qtlMarkersCP=qtlMarkersCP, mqtl=mqtl, fqtl=fqtl, aqtl=aqtl,rqtl=rqtl,CIs=CIs))
}    
########################################################################################################################

###### convert rownames of cross (marker information) to  a data frame ###########################
parse.rownames = function(rn) {
    rrn = data.frame(do.call('rbind', strsplit(rn, '_')),stringsAsFactors=F)
    rrn[,1] = as.numeric(rrn[,1])
    rrn$nchrom = rrn[,2]
    rrn$nchrom = gsub('chr0', 'chr', rrn$nchrom)
    rrn$nchrom = gsub('chr', '', rrn$nchrom)
    rrn$nchrom = as.numeric(rrn$nchrom)
    rrn[,3] = as.numeric(rrn[,3])
    names(rrn) = c('GPOS', 'CHR', 'POS', 'REF', 'ALT', 'NCHR')
    return(rrn)  }
###################################################################################################

###################################################################################################
LODmatrix.2.scanone= function(LODS, cross, LL=NULL) {
    if(is.null(LL)){
        LODSm = t(as.matrix(LODS))
        LODSs = scanone(cross, pheno.col=3, method='mr')
        LODSso = data.frame(LODSs, LODSm)
        LODSso= LODSso[,-3]
        class(LODSso)=class(LODSs)
        return(LODSso)
    } else { 
        LODSso =  data.frame(LL[,c(1,2)],t(as.matrix(LODS)))
        class(LODSso)  = class(LL)
        return(LODSso)
    }
}
###################################################################################################

##### make CI data frame ##########################################################################
makeCI.df = function(CIs) {
        CIs.df = do.call('rbind',CIs)
        
        CIs.peak  = CIs.df[seq(2,nrow(CIs.df),3),]
        CIs.peak$Marker = rownames(CIs.peak)
        CIs.peak$cpos  = parse.rownames(CIs.peak$Marker)$POS
        CIs.peak$gpos  = parse.rownames(CIs.peak$Marker)$GPOS
        colnames(CIs.peak)= paste('peak',colnames(CIs.peak), sep='')
        
        CIs.left  = CIs.df[seq(1,nrow(CIs.df),3),c(2,3)]
        CIs.left$Marker = rownames(CIs.left)
        CIs.left$cpos  = parse.rownames(CIs.left$Marker)$POS
        CIs.left$gpos  = parse.rownames(CIs.left$Marker)$GPOS
        colnames(CIs.left)= paste('left',colnames(CIs.left), sep='')
        
        CIs.right = CIs.df[seq(3,nrow(CIs.df),3),c(2,3)]
        CIs.right$Marker = rownames(CIs.right)
        CIs.right$cpos  = parse.rownames(CIs.right$Marker)$POS
        CIs.right$gpos  = parse.rownames(CIs.right$Marker)$GPOS
        colnames(CIs.right)= paste('right', colnames(CIs.right), sep='')
        
        CIs = data.frame(CIs.peak, CIs.left, CIs.right)
        return(CIs)
}

#### Convert Peak List to Peak Array ... fix type problems and expand physical marker information #####################
makePeakArray = function(lodPeaks) {
    lodPeakArray         = peak.2.array(lodPeaks)
    lodPeakArray         = na.omit(lodPeakArray)
    lodPeakArray$chr     = as.numeric(lodPeakArray$chr)
    lodPeakArray$lod     = as.numeric(lodPeakArray$lod)
    lodPeakArray$peak.cM = as.numeric(lodPeakArray$peak.cM)
    lodPeakArray$inf.cM  = as.numeric(lodPeakArray$inf.cM)
    lodPeakArray$sup.cM  = as.numeric(lodPeakArray$sup.cM)

    peak.mi =  parse.rownames(lodPeakArray$mname.peak)[,c(1,3)]
    names(peak.mi)=paste('peak', names(peak.mi), sep='')
    inf.mi  =  parse.rownames(lodPeakArray$mname.inf)[,c(1,3)]
    names(inf.mi)=paste('inf', names(inf.mi), sep='')
    sup.mi  =  parse.rownames(lodPeakArray$mname.sup)[,c(1,3)]
    names(sup.mi)=paste('sup', names(sup.mi), sep='')
    lodPeakArray = data.frame(lodPeakArray, peak.mi, inf.mi, sup.mi, stringsAsFactors=FALSE)
}

##### Make QTL Model ###################################################################################################
qtlModel  = function(traitPeakArray) {
    traitName =  unique(traitPeakArray$trait)
    pheno.colIndex = match(traitName,names(fcross$pheno))
    
    mqtls = makeqtl(fcross, chr=traitPeakArray$chr, pos=traitPeakArray$peak.cM, what='prob')
    fqtl  = fitqtl(fcross, pheno.col =pheno.colIndex, mqtls, method='hk', get.ests=TRUE, dropone=TRUE)
    
    pheno = pdata[,pheno.colIndex]
    geno  = gdata[,as.character(traitPeakArray$mname.peak)]

    pearsonR = cor(pheno, geno, use='pairwise.complete.obs')
    pearsonR2 =pearsonR^2

    LOD=(-sum(!is.na(pheno))*log(1-pearsonR2))/(2*log(10)) 
    attr(LOD, 'trait')=traitName
    
    if(nrow(traitPeakArray)==1) {aqtl=NA} 
    else { aqtl = addint(fcross, pheno.col=pheno.colIndex, mqtls, method='hk', qtl.only=TRUE)}
    return(list(mqtls=mqtls, fqtl=fqtl,pearsonR=pearsonR, pearsonR2=pearsonR2, LOD=LOD, aqtl=aqtl))
}
########################################################################################################

chrPeakFinder = function(x,y,z, LOD.threshold, peak.radius ) {
    peak.lods=c()
    peak.ind=c()
    maxLOD=max(y)
    while(maxLOD>LOD.threshold) {
        maxLOD.cind =which.max(y)
        maxLOD.cpos =x[maxLOD.cind]
        maxLOD.ind  =z[maxLOD.cind]
        l.ind = findInterval(maxLOD.cpos-peak.radius, x, all.inside=T, rightmost.closed=T)
        r.ind = findInterval(maxLOD.cpos+peak.radius, x)
        y[l.ind:r.ind]=0
        peak.lods= c(peak.lods, maxLOD)
        peak.ind = c(peak.ind, maxLOD.ind)
        maxLOD = max(y)
    }
   return( cbind(peak.ind, peak.lods) )     
}
