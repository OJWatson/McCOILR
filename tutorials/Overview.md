# The REAL McCOIL Rcpp package
OJ Watson  
`r format(Sys.time(), '%d %B, %Y')`  
    
## Overview
    
1. Description of aims of package
2. Testing Rcpp package
3. Reducing SNP sites in Uganda study



## 1. Data from Citygrapher 1 

*McCOILR* is simply an Rcpp implementation of [THEREALMcCOIL] (https://github.com/Greenhouse-Lab/THEREALMcCOIL),
which was written solely to make running the software easier within the cluster
framework I use. All rights refer to the writers of the original c code.

The package can be installed as follows, assuming devtools has been installed. 


```r
## first let's install the package
devtools::install_github("OJWatson/McCOILR")

## Load the package
library(McCOILR)
```

## 2. Testing Rcpp package

The package carries out the same 2 R functions as before, which are demonstrated below: 


```r
## categorical test

# read in demo data and view it
data0 = read.table(system.file("extdata","cat_input_test.txt",package="McCOILR"), head=T)
data=data0[,-1]
rownames(data)=data0[,1]

# view the heterozygosity calls
head(data)
```

```
##      site1 site2 site3 site4 site5 site6 site7 site8 site9 site10 site11
## ind1   1.0   0.5     1     0     1   0.0    -1  -1.0  -1.0      0    1.0
## ind2   0.0   0.5     1     1     1   0.0    -1   0.0  -1.0      0    0.5
## ind3   0.5   0.0     1     1    -1   0.0     1   0.5  -1.0      0    0.5
## ind4   0.5   0.0     1    -1     1   0.0     1   0.0  -1.0      0    1.0
## ind5   0.5   0.0     1     1     1   0.0     1   0.5   0.5      0    0.5
## ind6  -1.0   0.0     1    -1    -1   0.5    -1   0.0  -1.0      0    1.0
##      site12 site13 site14 site15 site16 site17 site18 site19 site20 site21
## ind1      0    1.0    1.0     -1      1      0      1      0     -1    0.5
## ind2     -1    0.5    1.0     -1      1     -1      1     -1     -1    0.0
## ind3      0    0.5    0.0     -1      1      0      1     -1     -1    0.0
## ind4      0    0.0    1.0     -1      1      0      1      0      1    0.0
## ind5      0    0.5    0.5      0      1     -1      1     -1      1    0.5
## ind6     -1   -1.0   -1.0     -1      1     -1     -1     -1     -1   -1.0
##      site22 site23 site24 site25 site26 site27 site28 site29 site30 site31
## ind1   -1.0      1    0.5      1   -1.0   -1.0   -1.0      1    0.5      1
## ind2    1.0      1    0.5      1   -1.0    0.0   -1.0      1    0.5      1
## ind3    1.0      1    0.5      1   -1.0    1.0    1.0      1    0.5      1
## ind4    1.0      1    1.0      1    0.5    1.0    0.5      1    1.0      1
## ind5    0.5      1    1.0      1    0.5    0.5    0.5      1    0.0      1
## ind6    1.0      1   -1.0      1   -1.0   -1.0   -1.0      1    0.5     -1
##      site32 site33 site34 site35 site36 site37 site38 site39 site40 site41
## ind1      1      1    0.5      1   -1.0      1    0.0      0     -1      0
## ind2      1      1    0.5      0   -1.0     -1    0.5      1     -1     -1
## ind3      1      1    0.0      1    0.0     -1    0.5      0      1     -1
## ind4      1      1    0.0      1    0.0     -1    0.5      0     -1      0
## ind5      1      1    1.0      1    0.5      1    0.0      0      0      0
## ind6     -1     -1    1.0     -1   -1.0     -1   -1.0      0     -1     -1
##      site42 site43 site44 site45 site46 site47 site48 site49 site50 site51
## ind1   -1.0      0      1    1.0     -1    0.5      0    0.5     -1      0
## ind2   -1.0      0      1    1.0     -1    1.0      0   -1.0      0      0
## ind3   -1.0      0      1    0.5     -1    1.0      0   -1.0      0      0
## ind4   -1.0      0      1    0.0     -1    0.5      0    1.0      0      0
## ind5    0.5      0      1    0.5      0    0.5      0    1.0      0      0
## ind6   -1.0      0      1   -1.0     -1   -1.0     -1   -1.0      0      0
##      site52 site53 site54 site55 site56 site57 site58 site59 site60 site61
## ind1    0.5      0    0.0    0.5      1      1     -1    1.0     -1    0.0
## ind2    1.0      0    0.5    0.0      1      1     -1    0.5     -1    0.0
## ind3   -1.0      0    0.0    0.0     -1      1     -1    0.5     -1    0.0
## ind4    0.0      0    0.0    0.0      1      1      1    1.0      1    0.0
## ind5    0.5      0    0.5    0.0      1      1      1    0.5      1    0.0
## ind6   -1.0      0   -1.0    0.0     -1      1     -1   -1.0     -1    0.5
##      site62 site63 site64 site65 site66 site67 site68 site69 site70 site71
## ind1      1      0    1.0   -1.0      1      1     -1    0.0    0.0      0
## ind2      1      0    1.0   -1.0      1     -1      0    0.0    0.5      0
## ind3      1      0    1.0    0.5      1     -1      0    0.5    1.0      0
## ind4      1      0    1.0   -1.0      1      1      1    0.0    0.0      0
## ind5      1      0    1.0    0.0      1      1     -1    0.0    0.0      0
## ind6      1      0    0.5   -1.0     -1     -1      0   -1.0    0.5     -1
##      site72 site73 site74 site75 site76 site77 site78 site79 site80 site81
## ind1      1      0     -1      1     -1    0.0    0.5    0.0      0    0.0
## ind2     -1      0     -1      1     -1   -1.0    0.5    1.0     -1    0.5
## ind3     -1      0     -1      1     -1   -1.0    0.5    0.5     -1    0.0
## ind4     -1      0      1      0     -1    1.0    1.0    1.0      0    0.0
## ind5      1      0     -1      1      1    0.5   -1.0    0.5     -1    0.0
## ind6     -1      0     -1      1     -1   -1.0    0.0    0.0     -1    0.0
##      site82 site83 site84 site85 site86 site87 site88 site89 site90 site91
## ind1    0.0    0.5      0      0      0      0      1      1      0      0
## ind2    0.0   -1.0      0      0     -1      0      1     -1     -1     -1
## ind3    0.5   -1.0      0      0     -1      0      1     -1      0     -1
## ind4    0.0    1.0      0      0      0      0      1      1      0     -1
## ind5    0.0    0.5      0      0      0      0      1      1      0      0
## ind6   -1.0   -1.0      0     -1     -1      0     -1     -1     -1     -1
##      site92 site93 site94 site95 site96 site97 site98 site99 site100
## ind1      0   -1.0      1      0      1    0.5   -1.0     -1      -1
## ind2      0   -1.0      1      0      1    0.0    1.0      0       1
## ind3      0   -1.0      1      0      1    1.0    0.5      0       1
## ind4      0    0.5      1      0      1    1.0    1.0      0       1
## ind5      0    0.0      1      0      1    0.5    0.5     -1      -1
## ind6      0   -1.0      1     -1      1    1.0   -1.0      0       1
##      site101 site102 site103 site104 site105
## ind1     0.0       1       1      -1      -1
## ind2     0.0       1       1      -1      -1
## ind3     0.0       1       1      -1      -1
## ind4     0.0       1       1       0      -1
## ind5     0.0       0       1       0       1
## ind6     0.5      -1      -1       0      -1
```

```r
# create results directory and run analysis
dir.create(path = "cat_output")
out_cat <- McCOIL_categorical(data,maxCOI=25, threshold_ind=20, threshold_site=20, totalrun=1000, burnin=100, M0=15, e1=0.05, 
                   e2=0.05, err_method=3, path="cat_output", output="output_test.txt" )
```

```
## Time = 1.00 s
```

```r
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
```

```
##         site1    site2     site3    site4     site5    site6    site7
## ind1 0.773019 0.435251  6.523575 5.529429 10.301879 0.635886 0.000000
## ind2 0.997224 0.336538  0.000000 0.254333  0.642381 0.535586 0.061059
## ind3 4.476989 1.684320  2.130197 4.353612  0.524203 0.249656 0.080993
## ind4 7.815435 3.313494 12.697070 8.285221  3.011241 1.013239 0.037105
## ind5 7.931132 0.442517  0.167996 2.177373  0.587036 0.300895 0.030904
## ind6 2.458607 0.638248  3.887696 6.067082  9.065840 0.905267 0.081333
##         site8    site9    site10    site11   site12   site13   site14
## ind1 4.634959 0.195470 13.927895  2.363394 5.382498 3.024445 5.423714
## ind2 6.542795 0.039037  0.000000 -1.000000 0.754422 3.324416 0.255289
## ind3 2.019209 0.000000 11.818180  0.069298 0.673240 3.098849 1.631299
## ind4 6.032118 0.000000  8.047423  0.000000 7.221847 2.497320 3.533402
## ind5 1.024954 0.081224  6.033347 -1.000000 1.568639 1.279992 0.163047
## ind6 4.586606 0.311690  7.699633  2.643740 2.744867 3.400377 6.825433
##        site15    site16   site17   site18   site19   site20 site21 site22
## ind1 0.000000  7.219236 0.000000 0.321880 1.035684 0.000000      0      0
## ind2 0.001773  0.153461 0.000000 0.191307 0.460990 0.000000      0     -1
## ind3 0.030687  2.452794 0.000000 0.459942 0.837007 0.229977      0      0
## ind4 0.000000 11.100150 0.161285 0.703449 1.271100 0.000000      0      0
## ind5 0.159506  3.586462 0.064790 0.190289 0.379542 0.000000      0      0
## ind6 0.000000  5.841971 0.176968 0.372723 1.110334 0.000000      0      0
##         site23   site24    site25   site26    site27    site28    site29
## ind1  0.000000 6.416304 16.064404 0.197101 10.157699  0.212330  9.499981
## ind2 -1.000000 0.000000  1.726730 0.190073  0.000000 -1.000000 -1.000000
## ind3  0.000000 4.429880 11.425260 0.232608  6.536019  0.312552  0.708249
## ind4  0.103252 7.916089  8.395237 0.035063 10.168160  0.157301  8.480841
## ind5  0.000000 1.037307  4.754813 0.007692  3.858720  0.184479  3.152692
## ind6  0.000000 8.487780 13.421972 0.144583  8.898299  0.268645  3.192409
##        site30   site31   site32   site33   site34   site35   site36
## ind1 9.835422 2.055015 5.926235 7.078459 0.000000 0.789787 0.815353
## ind2 4.855352 0.000000 0.000000 0.384062 0.000000 0.330539 0.350350
## ind3 8.368239 2.397373 9.008590 0.639089 0.019107 0.256107 0.391928
## ind4 6.388393 7.613878 6.477691 2.964806 0.594425 1.565033 0.501707
## ind5 5.864166 0.316654 0.477997 0.202001 0.000000 0.321333 0.275677
## ind6 9.806597 5.342317 7.105007 7.560561 0.000000 4.308457 0.590327
##        site37    site38   site39   site40   site41   site42    site43
## ind1 0.000000  0.000000 0.261659 6.073929 0.000000 0.000000  2.331914
## ind2 0.218452  0.198614 0.646294 0.250292 1.099951 0.041389  2.158835
## ind3 0.000000  0.157776 0.181682 1.581731 0.082319 0.000000  2.443111
## ind4 0.476174  0.016512 2.729292 4.065335 0.599497 0.000000 10.687130
## ind5 7.588091 -1.000000 0.005372 2.848523 0.271334 0.000000  0.201186
## ind6 0.000000  0.000000 0.805857 7.973652 0.000000 0.000000  5.407707
##        site44   site45   site46    site47   site48    site49   site50
## ind1 0.000000 3.208030 1.332707  7.534703 0.332472  8.247594 1.611010
## ind2 0.000000 3.220459 0.940977  7.443709 0.737532  8.433792 2.815520
## ind3 0.000000 3.197772 0.273793  3.723208 0.288094  2.569323 4.338393
## ind4 0.000000 4.860701 0.768297 10.746410 0.663363 11.933540 6.475605
## ind5 2.899116 1.124923 0.000000  4.235002 0.380688  4.981566 3.503463
## ind6 0.000000 3.657817 1.061532  6.957743 0.463582  6.632334 4.134599
##         site51    site52 site53   site54   site55   site56   site57
## ind1  0.211831  2.414855      0 2.047866 5.552769 3.515977 0.000000
## ind2 -1.000000  0.073066      0 2.474965 8.368501 0.000000 0.000000
## ind3  0.107921  0.256071      0 1.815254 2.822285 1.317744 2.633277
## ind4  0.000000  5.014762      0 3.182412 9.823132 5.101557 1.148384
## ind5  0.200261 -1.000000      0 2.672133 5.945422 1.046215 0.000000
## ind6  0.333725  1.689286      0 3.011802 6.636636 2.838735 0.488921
##         site58   site59    site60   site61   site62   site63    site64
## ind1  5.617321 5.567360  5.889082 0.407831 4.735254 0.212227  4.301725
## ind2 -1.000000 4.662254  6.840784 0.643242 0.000000 0.378215  0.000000
## ind3  2.067533 3.208334  4.003303 0.183730 0.435801 0.000000 -1.000000
## ind4  4.601010 9.515062  7.887465 2.409367 2.435037 0.698486 12.837560
## ind5 -1.000000 0.136660  4.014484 0.299023 0.272797 0.719106  0.000000
## ind6  9.737784 8.436536 10.875579 0.642340 5.458116 0.009906  8.947218
##        site65   site66 site67   site68    site69   site70 site71    site72
## ind1 0.086193 0.195682      0 3.923784  3.012980 0.084219      0  4.481197
## ind2 0.000000 1.024617     -1 4.139583 -1.000000 0.683016      0  0.080532
## ind3 0.139899 0.000000      0 3.054424  0.778499 2.766324      0  2.581468
## ind4 0.430866 0.102573      0 5.653807  1.565004 4.702235      0  1.237159
## ind5 0.000000 0.494720      0 0.185001  0.102478 0.000000      0 -1.000000
## ind6 0.184573 0.226434      0 8.391926  5.339827 0.406806      0  7.603121
##        site73    site74    site75   site76   site77    site78   site79
## ind1 0.736700  7.174779  5.692942 0.117824 5.038420 -1.000000 1.222317
## ind2 0.740144  0.000000 -1.000000 0.180125 4.277511  0.000000 0.579069
## ind3 0.544233  4.176406  0.092368 0.000000 0.442165  0.040142 0.361243
## ind4 1.276951 12.122590  7.206281 0.157702 0.510834  0.000000 2.107665
## ind5 1.091841  0.000000 -1.000000 0.000000 3.393604 -1.000000 1.105395
## ind6 1.407879 10.076725  9.867932 0.212253 3.497185  0.012627 1.730980
##        site80   site81   site82    site83   site84   site85   site86
## ind1 0.074298 2.925244 4.445084  1.221600 0.000000 2.218825 2.304839
## ind2 2.361009 0.299759 0.000000 -1.000000 0.035173 0.000000 0.000000
## ind3 0.021137 0.116697 3.949225  0.000000 0.664386 0.000000 0.000000
## ind4 0.097928 0.167922 5.304122  0.142190 0.185392 1.107071 0.219858
## ind5 0.273618 0.204405 1.739071  0.000000 0.125757 1.898828 0.036527
## ind6 6.016359 0.184774 3.817465  0.528938 0.341512 2.321284 2.074875
##        site87   site88   site89   site90   site91   site92   site93
## ind1 0.000000 5.781570 0.731047 2.762129 2.193086 5.669109 4.948302
## ind2 3.341576 1.042625 0.570464 1.512688 0.000000 0.206961 2.146975
## ind3 2.694466 0.672151 2.886937 1.046025 2.864557 5.115930 3.521272
## ind4 3.884768 4.873755 4.610479 3.401538 1.075722 6.348375 5.779720
## ind5 4.599888 4.351057 0.915393 2.209481 0.402474 2.972046 3.850597
## ind6 1.984991 5.338447 5.046478 2.476585 4.482983 9.222386 4.920970
##        site94   site95   site96    site97   site98    site99   site100
## ind1 0.239401 3.890739  0.88387  0.000000 2.115714  1.017555  0.133018
## ind2 0.000000 0.000000 -1.00000 -1.000000 1.289655  0.000000 -1.000000
## ind3 0.113925 0.505523  0.00000  2.245503 1.295489  1.101571  0.198458
## ind4 0.564348 2.553972  0.00000  0.204393 3.732895  0.682663  0.219620
## ind5 0.000000 0.269705 -1.00000  0.000000 2.773473 -1.000000  0.000000
## ind6 3.745815 3.310597  0.00000  1.290676 3.431699  2.073178  0.273429
##       site101  site102   site103  site104   site105  site106  site107
## ind1 5.686296 3.520340  1.448558 5.022188  2.554234 4.375710  0.00000
## ind2 2.571777 0.482079 -1.000000 3.015142 -1.000000 1.412338 -1.00000
## ind3 2.507220 3.057959  4.206432 2.439314  0.000000 2.534611  1.08814
## ind4 5.440795 4.996793  7.513086 2.248669  0.375817 5.752069  0.00000
## ind5 3.429925 0.861018  0.000000 0.024655 -1.000000 3.745864  0.00000
## ind6 4.014097 3.320495  5.106893 2.776721  3.485361 4.850756  0.00000
```

```r
# create results directory and run analysis
dir.create(path="prop_output")
out_prop <- McCOIL_proportional(dataA1, dataA2, maxCOI=25, totalrun=5000, burnin=100, M0=15, epsilon=0.02, err_method=3, 
                                path="prop_output", output="output_test.txt" )
```

```
## Iter 500 out of 5000
## Iter 1000 out of 5000
## Iter 1500 out of 5000
## Iter 2000 out of 5000
## Iter 2500 out of 5000
## Iter 3000 out of 5000
## Iter 3500 out of 5000
## Iter 4000 out of 5000
## Iter 4500 out of 5000
## Iter 5000 out of 5000
## Time = 35.00 s
```

The R functions now return the summary outputs, just for ease of looking at the 
results:


```r
## view summary data.frame for categorical
str(out_cat)
```

```
## 'data.frame':	90 obs. of  8 variables:
##  $ file         : Factor w/ 1 level "output_test.txt": 1 1 1 1 1 1 1 1 1 1 ...
##  $ CorP         : Factor w/ 4 levels "C","e1","e2",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ name         : Factor w/ 90 levels "e1","e2","ind1",..: 3 14 25 30 31 32 33 34 35 4 ...
##  $ mean         : Factor w/ 63 levels "0.0132025622222222",..: 59 59 59 58 60 59 60 60 58 62 ...
##  $ median       : Factor w/ 63 levels "0.0098235","0.010386",..: 59 59 59 58 60 59 60 60 58 62 ...
##  $ sd           : Factor w/ 89 levels "0","0.00163",..: 61 63 73 1 75 72 74 69 64 80 ...
##  $ quantile0.025: Factor w/ 61 levels "0.000145","0.000502",..: 59 59 59 58 59 59 59 59 58 61 ...
##  $ quantile0.975: Factor w/ 70 levels "0.04479","0.05231305",..: 64 65 67 58 67 66 67 66 64 70 ...
```

```r
## Have a look at the categorical output COI distribution
hist(as.numeric(as.character(out_cat$mean[out_cat$CorP=="C"])),main = "Categorical Mean COI", xlab="COI")
```

![](Overview_files/figure-html/View output-1.png)<!-- -->

```r
## view summary data.frame for proportional
str(out_prop)
```

```
## 'data.frame':	144 obs. of  8 variables:
##  $ file         : Factor w/ 1 level "output_test.txt": 1 1 1 1 1 1 1 1 1 1 ...
##  $ CorP         : Factor w/ 3 levels "C","epsilon",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ name         : Factor w/ 144 levels "epsilon","ind1",..: 2 13 24 32 33 34 35 36 37 3 ...
##  $ mean         : Factor w/ 114 levels "0.0100813442857143",..: 109 109 109 113 109 111 114 109 109 109 ...
##  $ median       : Factor w/ 114 levels "0.006793","0.007087",..: 109 109 109 113 109 111 114 109 109 109 ...
##  $ sd           : Factor w/ 139 levels "0","0.00228",..: 119 114 34 138 120 131 139 111 56 123 ...
##  $ quantile0.025: Factor w/ 113 levels "0.000132","0.000233",..: 107 107 107 110 107 109 112 107 107 107 ...
##  $ quantile0.975: Factor w/ 117 levels "0.032866775",..: 111 111 110 117 111 113 109 110 110 111 ...
```

```r
## Have a look at the proportional output COI distribution
hist(as.numeric(as.character(out_prop$mean[out_cat$CorP=="C"])),main = "Proportional Mean COI", xlab="COI")
```

![](Overview_files/figure-html/View output-2.png)<!-- -->

## 3. Reducing SNP sites in Uganda study

The main reason the package was created was so that the software could be easily
run on my cluster framework in my department. This was because I wanted to see
what effect reducing the number of SNP locations would have in terms of 
increasing or decreasing the estimated COI. 


```r
# Read in Uganda dataset
snps <- read.csv(system.file("extdata","SNP.txt",package="McCOILR"),sep="\t")

# let's just look at Kihihi
kihihi_data <-  snps[snps$location=="Kihihi",]

# create 100 sets of 24, 48 and 96 SNP site samples
barcode24 <- replicate(sample(3:105,size = 24,replace = FALSE),n = 100)
barcode48 <- replicate(sample(3:105,size = 48,replace = FALSE),n = 100)
barcode96 <- replicate(sample(3:105,size = 96,replace = FALSE),n = 100)

## The following is code to run the categorical analysis on our cluster but is 
## just looping through the above barcode sets and running the categorical analysis
## on each one

## --------------------------------------------------------------------------##

workdir <- "M:/OJ/McCOILR_Results"
didehpc::didehpc_config_global(workdir=workdir,
                               credentials="C:\\Users\\Oliver\\.smbcredentials",
                               temp=didehpc::path_mapping("tmp","T:","//fi--didef3.dide.ic.ac.uk/tmp","T:"),
                               home=didehpc::path_mapping("OJ","M:","//fi--didef3.dide.ic.ac.uk/Malaria","M:"))
didehpc::web_login()
root <- file.path(workdir, "contexts")
packages.vector <- c("Rcpp","McCOILR")
context::context_log_start()
ctx <- context::context_save(root,
                             packages = packages.vector,
                             package_sources= provisionr::package_sources(github=c("OJWatson/McCOILR")))

config <- didehpc::didehpc_config(use_workers = TRUE)
obj <- didehpc::queue_didehpc(ctx, config = config)

paramList <- list()
length(paramList) <- 300

for(i in 1:300){

    if(i <=100){
        data <- kihihi_data[,barcode24[,i]+2]
        output <- paste0("Kihihi_24BarcodeTest_50maxCOI/24_",i,".txt")
    } else if (i > 100 & i <= 200) {
        data <- kihihi_data[,barcode48[,i-100]+2]
        output <- paste0("Kihihi_48BarcodeTest_50maxCOI/48_",i-100,".txt")
    } else {
        data <- kihihi_data[,barcode96[,i-200]+2]
        output <- paste0("Kihihi_96BarcodeTest_50maxCOI/96_",i-200,".txt")
    }
    

paramList[[i]] <- (list(data = data, maxCOI=50, threshold_ind=20, 
                       threshold_site=20, totalrun=100000, burnin=100,
                       M0=15, e1=0.05, e2=0.05, err_method=3, 
                       path="M:/OJ/McCOILR_Results",
                        output=output))

}

dir.create(paste0(workdir,"/Kihihi_24BarcodeTest_50maxCOI"))
dir.create(paste0(workdir,"/Kihihi_48BarcodeTest_50maxCOI"))
dir.create(paste0(workdir,"/Kihihi_96BarcodeTest_50maxCOI"))

workers <- obj$submit_workers(50)
grp <- queuer::qlapply(X = paramList[1:100],obj = obj, timeout=0, 
                       FUN = function(x){return(McCOIL_categorical(data = x$data,maxCOI = x$maxCOI, threshold_ind = x$threshold_ind,
                                                                   threshold_site = x$threshold_site, totalrun = x$totalrun, 
                                                                   burnin = x$burnin, M0 = x$M0, e1 = x$e1, e2 = x$e2, 
                                                                   err_method = x$err_method, path = x$path, output = x$output)
                                                )},
                       name = "24_barcode")

## --------------------------------------------------------------------------##
```

The code above is executing the categorical method on 100 different randomly chosen
sets of SNP sites, with sets equal to 24, 48, 96, with 100,000 MCMC iterations 
and a maximum possible COI of 50.

The summary results of these MCMC iterations has been copied into the tutorials folder, 
so the following analysis can be verified. Below I'm simply showing the overall 
distribution for mean COI for the sets of 24, 48 and 96 SNP sites. 


```r
# storage for the 24, 48 and 96 site datasets
twenty4 <- list()
fourty8 <- list()
ninety6 <- list()

for(i in 1:100){
    
  twenty4[[i]] <- read.csv(paste0("Kihihi_24_BarcodeTest_50maxCOI/24_",i,".txt_summary.txt"),sep="\t")
  fourty8[[i]] <- read.csv(paste0("Kihihi_48_BarcodeTest_50maxCOI/48_",i,".txt_summary.txt"),sep="\t")
  ninety6[[i]] <- read.csv(paste0("Kihihi_96_BarcodeTest_50maxCOI/96_",i,".txt_summary.txt"),sep="\t")
  
  

}

## Grab the estimated mean COI and lets plot these
mean24 <- lapply(twenty4,function(x){return(as.numeric(as.character(x$mean[x$CorP=="C"])))})
mean48 <- lapply(fourty8,function(x){return(as.numeric(as.character(x$mean[x$CorP=="C"])))})
mean96 <- lapply(ninety6,function(x){return(as.numeric(as.character(x$mean[x$CorP=="C"])))})

par(mfrow=c(1,3))
hist(as.numeric(unlist(mean24)),main="24 Sites",xlab="COI")
hist(as.numeric(unlist(mean48)),main="48 Sites",xlab="COI")
hist(as.numeric(unlist(mean96)),main="96 Sites",xlab="COI")
```

![](Overview_files/figure-html/Plot pseudo-steiner-again-1.png)<!-- -->

Initially it seemed to me that the fewer SNP sites was causing a larger COI to be
estimates. However, upon further inspection of the chains it could be seen that 
using fewer SNP sites resulted in MCMC chains failing to converge more often, which
gave spurious results.


```r
## First let's grab one MCMC output with 24 sites
MCMC_24 <- read.table("M:/OJ/McCOILR_Results/Kihihi_24BarcodeTest_50maxCOI/24_99.txt",head=F)

## First let's grab one MCMC output with 96 sites
MCMC_96 <- read.table("M:/OJ/McCOILR_Results/Kihihi_96BarcodeTest_50maxCOI/96_99.txt",head=F)

## We'll save these in the tutorials folder too for reproducibility
save(list=c("MCMC_24","MCMC_96"),file = "MCMC24&96.RData")

## Now let's plot some random COI estimate chains

# for MCMC24 we will plot 9 random COI and allele frequency chains with some thinning
par(mfrow=c(3,3))
for(i in floor(seq(2,length(MCMC_24)-2,length.out = 9))){
    plot(MCMC_24[[i]][seq(1,100000,50)],type="l",
         ylab = "COI", main = twenty4[[99]]$name[i-1])
}
```

![](Overview_files/figure-html/Plot example chains-1.png)<!-- -->

```r
# and the same for MCMC96
par(mfrow=c(3,3))
for(i in floor(seq(2,length(MCMC_96)-2,length.out = 9))){
    plot(MCMC_96[[i]][seq(1,100000,50)],type="l",
         ylab = ifelse(length(grep(pattern = "assay",as.character(ninety6[[99]]$name[i-1]))==1),"Freq","COI"),
         main = ninety6[[99]]$name[i-1])
}
```

![](Overview_files/figure-html/Plot example chains-2.png)<!-- -->

We can also see that when we use the 96 SNP sites, the overall mean COI is equal to
2.60, with a median equal to 2 which agrees very well with the analysis presented of 
the categorical method for Kihihi within the supplementary material for THE REAL
McCOIL, whereas for 24 sites these do not agree:


```r
## Mean for 24
print( paste0("Mean for 24 = ", mean(as.numeric(unlist(mean24)))))
```

```
## [1] "Mean for 24 = 8.94219554030875"
```

```r
## Median for 24
print( paste0("Median for 24 = ", median(as.numeric(unlist(mean24)))))
```

```
## [1] "Median for 24 = 5"
```

```r
## Mean for 96
print( paste0("Mean for 96 = ", mean(as.numeric(unlist(mean96)))))
```

```
## [1] "Mean for 96 = 2.59792"
```

```r
## Median for 96
print( paste0("Median for 24 = ", median(as.numeric(unlist(mean96)))))
```

```
## [1] "Median for 24 = 2"
```

As such it makes sense in the future when using THE REAL McCOIL on new datasets
to use more sites within computational reason that display both polymorphic and 
intermediate/high allele frequencies in multiple populations, such as the 128 SNPs 
used within the publication attached to THE REAL McCOIL. 

In addition I will be adding functionality to McCOILR to pipeline chain convergence
assessment, such that when x number of repetitions on the same dataset do not yield 
convergent results using Gelman and Rubin convergence diagnostics additional repetitions
are automated until a convergent result is achieved. 
