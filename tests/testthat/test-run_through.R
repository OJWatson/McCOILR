test_that("pipeline", {
   
    # read in demo data and view it
    data0 = read.table(system.file("extdata","cat_input_test.txt",package="McCOILR"), head=T)
    data=data0[,-1]
    rownames(data)=data0[,1]
    
    # view the heterozygosity calls
    head(data)
    
    # create results directory and run analysis
    dir.create(path = "cat_output")
    out_cat <- McCOIL_categorical(data,maxCOI=25, threshold_ind=20, threshold_site=20,
                                  totalrun=1000, burnin=100, M0=15, e1=0.05, e2=0.05,
                                  err_method=3, path="cat_output", output="output_test.txt" )
    
    expect_named(out_cat, expected = c("file","CorP","name","mean","median","sd","quantile0.025","quantile0.975"))
    expect_true(any(grepl("ind", out_cat$name)))
    expect_true(any(grepl("site", out_cat$name)))
    expect_true(any(grepl("e\\d", out_cat$name)))
    
    ## proportional test
    
    # read in demo data and view it
    dataA1i = read.table(system.file("extdata","prop_dataA1_test.txt",package="McCOILR"), head=T)
    dataA2i = read.table(system.file("extdata","prop_dataA2_test.txt",package="McCOILR"), head=T)
    dataA1= dataA1i[,-1]
    dataA2= dataA2i[,-1]
    rownames(dataA1)= dataA1i[,1]
    rownames(dataA2)= dataA2i[,1]
    
    # view the read counts
    head(dataA1)
    
    # create results directory and run analysis
    dir.create(path="prop_output")
    out_prop <- McCOIL_proportional(dataA1, dataA2, maxCOI=25, totalrun=1000, burnin=100,
                                    M0=15, epsilon=0.02, err_method=3, path="prop_output",
                                    output="output_test.txt" )
    
    expect_named(out_prop, expected = c("file","CorP","name","mean","median","sd","quantile0.025","quantile0.975"))
    expect_true(any(grepl("ind", out_prop$name)))
    expect_true(any(grepl("site", out_prop$name)))
    expect_true(any(grepl("e\\d", out_prop$name)))
    
    
})
