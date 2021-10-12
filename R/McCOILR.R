#' Rcpp implementation of THE REAL McCOIL
#' 
#' @name McCOILR
#' @aliases McCOILR package-McCOILR
#' @docType package
#' @title Wrapper for THE REAL McCOIL in Rcpp
#' 
#' @importFrom stats median quantile sd
#' @importFrom utils read.table write.table head
#' @importFrom Rcpp sourceCpp
#' 
#' @useDynLib McCOILR, .registration = TRUE
#' 
#' @description Wrapper for THE REAL McCOIL in Rcpp, so that package can be more easily run on distributed computing services and cluster infrastructure. 
#' 
#' @references 1 Chang H-H, Worby CJ, Yeka A, Nankabirwa J, Kamya MR, Staedke SG, Dorsey G, Murphy M, Neafsey DE, Jeffreys AE, Hubbart C, Rockett KA, Amato R, Kwiatkowski DP, Buckee C, Greenhouse B. 2017. THE REAL McCOIL: A method for the concurrent estimation of the complexity of infection and SNP allele frequency for malaria parasites. PLOS Comput Biol 13: e1005348. doi:10.1371/journal.pcbi.1005348
#'


NULL
