#------------------------------------------------
#' The Real McCOIL proportional method function
#'
#' This function triggers the c code for the proportional method
#' 
#' @param dataA1 The intensity of signals of allele 1 from the SNP assay. Row names are names of samples and column names are names of assays.
#' @param dataA2 The intensity of signals of allele 2 from the SNP assay. Row names are names of samples and column names are names of assays.
#' @param maxCOI Upper bound for COI. The default is 25.
#' @param totalrun The total number of MCMC iterations. The default is 10000.
#' @param burnin The total number of burnin iterations. The default is 1000.
#' @param M0 Initial COI. The default is 15.
#' @param epsilon The level of measurement error (eest). The default is 0.2.
#' @param path The default is the current directory.
#' @param output The name of output file. The default is output.txt.
#' @param err_method The default is 1. 
#' 1: use pre-specified epsilon; 
#' 2: use likelihood-free sampling for epsilon; 
#' 3: update epsilon according to likelihood (for 2 and 3, pre-specified epsilon was used as initial value) 
#' 
#' 
#' @return summary of output as data.frame
#' 
#' @export


McCOIL_proportional = function(dataA1, dataA2, maxCOI = 25, totalrun = 10000, burnin = 1000, 
                               M0 = 15, epsilon = 0.02, err_method = 1, path = getwd(), output = "output.txt") {
    
    
    grid = beta_grid_25
    n = nrow(dataA1)
    k = ncol(dataA1)
    M0 = rep(M0, n)
    P0 = rep(0.5, k)
    A1 = as.vector(t(dataA1))
    A2 = as.vector(t(dataA2))
    
    paramList <- list(max = as.integer(maxCOI), iter = as.integer(totalrun), n0 = as.integer(n), 
                      k0 = as.integer(k), A1 = as.double(A1), A2 = as.double(A2), M0 = as.integer(M0), P0 = as.double(P0), 
                      A = as.double(grid$A), B = as.double(grid$B), c0 = as.double(epsilon), file_index = as.character(output), 
                      path = as.character(path), err_method0 = as.integer(err_method))
    
    
    
    
    if ((n > 10 & k > 10)) {
        McCOIL_proportional_cpp(paramList = paramList)
    } else {
        stop(paste("Sample size is too small (n=", n, ", k=", k, ").", sep = ""))
    }
    
    ## summarize results
    outputMCMC2 = read.table(paste(path, "/", output, sep = ""), head = F)
    meanM = as.numeric(round(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2, mean)))
    meanP = as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, mean))
    medianM = as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2, median))
    medianP = as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, median))
    M975 = as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2, 
                            function(x) quantile(x, probs = 0.975)))
    P975 = as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, 
                            function(x) quantile(x, probs = 0.975)))
    M025 = as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2,
                            function(x) quantile(x, probs = 0.025)))
    P025 = as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, 
                            function(x) quantile(x,  probs = 0.025)))
    sdM = as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, (1:n) + 1], 2, sd))
    sdP = as.numeric(apply(outputMCMC2[(burnin + 1):totalrun, ((1:k) + n + 1)], 2, sd))
    
    
    if (err_method == 3) {
        mean_e3 = as.numeric(mean(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)]))
        median_e3 = as.numeric(median(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)]))
        e3_975 = as.numeric(quantile(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)], probs = 0.975))
        e3_025 = as.numeric(quantile(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)], probs = 0.025))
        sd_e3 = as.numeric(sd(outputMCMC2[(burnin + 1):totalrun, (k + n + 2)]))
    }
    if ((err_method == 1) | (err_method == 2)) {
        output_sum = data.frame(cbind.data.frame(rep(output, (n + k)), c(rep("C", n), rep("P", k)),
                                      c(rownames(dataA1),colnames(dataA1)), 
                                      c(meanM, meanP), c(medianM, medianP),
                                      round(c(sdM, sdP), digits = 5), 
                                      c(M025, P025), c(M975, P975)))
    } else {
        output_sum = data.frame(cbind.data.frame(rep(output, (n + k + 1)), c(rep("C", n), rep("P", k),"epsilon"),
                                      c(rownames(dataA1), colnames(dataA1), "epsilon"), c(meanM, meanP, mean_e3), 
                                      c(medianM, medianP, median_e3), round(c(sdM, sdP, sd_e3), digits = 5), 
                                      c(M025, P025, e3_025), c(M975, P975, e3_975)))
    }
    colnames(output_sum) = c("file", "CorP", "name", "mean", "median", "sd", "quantile0.025","quantile0.975")
    
    write.table(output_sum, paste(path, "/", output, "_summary.txt", sep = ""), sep = "\t", 
                col.names = T, row.names = F, quote = F)
    
    return(output_sum)
    
}

