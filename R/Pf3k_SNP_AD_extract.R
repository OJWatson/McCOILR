#------------------------------------------------
#' Pf3k personal REF and MAF AD extraction
#'
#' This function returns the allele depth for a given chromosme and position
#' from the Pf3k database
#' 
#' @param chromosome Numeric giving chromosome number to search within
#' @param positions Numeric giving chromosome positions
#' @param ftp_address String for ftp url address. Default is "ftp://ngs.sanger.ac.uk/production/pf3k/release_5/5.1/"
#' @param qual Boolean detailing if QUAL for each SNP is extracted too. Default = FALSE
#' 
#' @return returns 1 list of n lists, where n is length of positions, and each 
#' of the n lists contains 3 lists (REF & MAF & QUAL) of length m, where m is the number
#' of sequences within Pf3k. Where more than one type of MAF is possible, 
#' the most common MAF is recorded.
#' 
#' @export
#' 

Pf3k_SNP_AD_extraction = function(chromosome, positions,
                                  ftp_address="ftp://ngs.sanger.ac.uk/production/pf3k/release_5/5.1/",
                                  qual = FALSE){
    
    ## Grab the filenames from the ftp address
    filenames <- RCurl::getURL(ftp_address, ftp.use.epsv = FALSE, dirlistonly = TRUE) %>% 
        strsplit("\r\n") %>% 
        unlist()
    
    ## look up the correct file for the chromosome specified
    if(chromosome < 10){
        pattern <- paste0("_0",chromosome,".*gz$")
        chrom <- paste0("Pf3D7_0",chromosome,"_v3")
    } else {
        pattern <- paste0("_",chromosome,".*gz$")   
        chrom <- paste0("Pf3D7_",chromosome,"_v3")
    }

    ## create full file path 
    vcffile <- paste0(ftp_address,filenames[grep(pattern,filenames)])
    
    ## sort positions
    positions <- sort(positions)
    
    ## create vcf params using specified chrom and positions
    param <- VariantAnnotation::ScanVcfParam(
        info=c("AF"),
        geno="AD",
        which = GenomicRanges::GRanges(chrom,IRanges::IRanges(positions,positions))
    )
    
    ## read in the vcf data
    vcf <- VariantAnnotation::readVcf(vcffile,genome=chrom, param=param)
    
    ## turn data into more useful list format
    AD <- as.list.data.frame(as.data.frame.matrix(t(VariantAnnotation::geno(vcf)[["AD"]])))
    
    ## extract REF and MAF
    AD_extraction <- lapply(AD,function(x){
        af <- lapply(x,function(y) y)
        ref <- lapply(af, function(y) y[1]) %>% unlist
        maf <- lapply(af, function(y){
            if(length(y) == 2){
                return(y[2])
            } else {
                return(max(y[2:length(y)]))
            }
        }) %>% unlist
        return(list("REF"=ref,"MAF"=maf))
    }) 
    
    if(qual == TRUE){
        QUAL <- qual(vcf)
        for(i in 1:length(AD)){
            AD[[i]]$QUAL = QUAL[i]
        }
    }
    
    ## return results
    return(AD_extraction)
}


#------------------------------------------------
#' Allele Depth Matrix Creation
#'
#' This function returns list of 2 matrices corrresponding to the allele depth for
#' REF and MAF, calculated from the output of \code{Pf3k_SNP_AD_extraction}
#' 
#' @param lookup object. List of lists, with each list having 2 elements
#' chromosome and positions. 
#' @param AD Output of \code{Pf3k_SNP_AD_extraction}
#' 
#' @return returns 1 list of 2 matrices representing the REF and MAF frequencies
#' required
#' 
#' @export
#' 

AD_mat_creation = function(AD, lookup){
    
    ## first work out how effective the AD extraction was and how many loci were found
    foundloci <- 0
    for(i in 1:14){
        newloci <- sum(!is.na(match(lookup[[i]]$positions,
                                    as.numeric(lapply(strsplit(names(AD[[i]]),":|_"),
                                                      function(x) x[4]) %>% unlist))),
                       na.rm=TRUE)
        
        foundloci <- foundloci + newloci
    }
    
    # set up matrices
    A1 <- A2 <- matrix(data = 0,nrow = length(AD[[1]][[1]]$REF),
                       ncol=foundloci,
                       dimnames = list(names(AD[[1]][[1]]$REF),1:foundloci))
    
    fillstart <- 1
    # Need to do matches as sometimes VCF has a field day and brings in extra SNPs
    
    for(i in 1:14){
        
        for(j in match(lookup[[i]]$positions,
                       as.numeric(lapply(strsplit(names(AD[[i]]),":|_"),function(x) x[4]) %>% unlist))){
            if(!is.na(j)){
                A1[,fillstart] <- AD[[i]][[j]]$REF
                A2[,fillstart] <- AD[[i]][[j]]$MAF
                A3 <- A1[,fillstart] + A2[,fillstart]
                non_reads <- which(A3==0)
                A1[non_reads,fillstart] <- -1
                A2[non_reads,fillstart] <- -1
                colnames(A1)[fillstart] <- colnames(A2)[fillstart] <- names(AD[[i]][j])
                fillstart <- fillstart + 1
            }
        }
    }
    
    
    
    ## return results
    return(list("A1"=A1,"A2"=A2))
}


