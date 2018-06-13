
#include "llfunction_prop.h"
#include <Rcpp.h>
#include <string>

//' @title
//' McCOIL_proportional_Rcpp
//' @description
//' Proportional cpp code that thins outputs, and stores outputs rather than 
//' writing to file. 
//' 
//' @param paramList A list of parameters created with equivalent R function
//' 
//' @details
//' \code{McCOIL_proportional_Rcpp} implements THE REAL McCOIL proportional
//'   method
//' 
//' @return Rcpp List
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List McCOIL_proportional_Rcpp(Rcpp::List paramList) {
    
    ////////can be changed/////////
    double varP=0.1;  
    double varS=0.1;
    double varE=0.0001;
    double upper_bound_c=0.1; //epsilon
    
    ///////don't change////////
    // both unused so commented out
    // double temp_l_est=0;
    // double temp_l_true=0;
    
    int i=0, j=0, prime=0, x=0, y=0;
    double sumori=0, sumcan=0; 
    
    std::cout << "1" << std::endl;
    
    ////read in parameter values////
    int max_moi= Rcpp::as<int>(paramList["max"]);	
    int iter = Rcpp::as<int>(paramList["iter"]);
    int n = Rcpp::as<int>(paramList["n0"]); 
    int k = Rcpp::as<int>(paramList["k0"]); 
    std::vector<double> A1 = Rcpp::as<std::vector<double> >(paramList["A1"]);
    std::vector<double> A2 = Rcpp::as<std::vector<double> >(paramList["A2"]);
    std::vector<int> M0 = Rcpp::as<std::vector<int> >(paramList["M0"]);
    std::vector<int> P0 = Rcpp::as<std::vector<int> >(paramList["P0"]);
    std::vector<double> A = Rcpp::as<std::vector<double> >(paramList["A"]);
    std::vector<double> B = Rcpp::as<std::vector<double> >(paramList["B"]);
    double c = Rcpp::as<double>(paramList["c0"]);
    int err_method = Rcpp::as<int>(paramList["err_method0"]); //1: use pre-specified e1 and e2; 2: use likelihood-free sampling for e1 and e2; 3: update e1 and e2 according to likelihood (for 2 and 3, pre-specified e1 and e2 were used as initial values) 
    int saves = Rcpp::as<int>(paramList["saves"]);
    
    
    std::cout << "2" << std::endl;
    
    ////read in the values////
    std::vector<int> Maccept(n+1), M(n+1), Mcan(n+1);
    std::vector<double> P(k+1), Pcan(k+1); 
    std::vector<int> Paccept(k+1);
    std::vector<std::vector<double> > ll(n+1, std::vector<double>(k+1));
    //double ll[(n+1)][(k+1)];
    std::vector<std::vector<double> > llcan(n+1, std::vector<double>(k+1));
    //double llcan[(n+1)][(k+1)];
    std::vector<std::vector<double> > gridA(26, std::vector<double>(51));
    //double gridA[26][51];
    std::vector<std::vector<double> > gridB(26, std::vector<double>(51));
    //double gridB[26][51];
    std::vector<std::vector<double> > dataA1(n+1, std::vector<double>(k+1));
    //double dataA1[(n+1)][(k+1)];
    std::vector<std::vector<double> > dataA2(n+1, std::vector<double>(k+1));
    //double dataA2[(n+1)][(k+1)];
    std::vector<std::vector<double> > Strue(n+1, std::vector<double>(k+1));
    //double Strue[(n+1)][(k+1)];
    std::vector<std::vector<double> > Strue_can(n+1, std::vector<double>(k+1));
    //double Strue_can[(n+1)][(k+1)];
    std::vector<std::vector<int> > Strue_accept(n+1, std::vector<int>(k+1));
    //int Strue_accept[(n+1)][(k+1)];
    
    int c_accept=0;
    double c_can;
    double q1=0.0, q2=0.0; // unused variables to be commented out TODO: Comment all unuseds
    
    for (i=1;i<=n;i++){
        M[i]= M0[i-1];
        Mcan[i]= M[i];
        Maccept[i]= 0;
        for (j=1;j<=k;j++){
            dataA1[i][j]= A1[(i-1)*k+j-1]; 
            dataA2[i][j]= A2[(i-1)*k+j-1]; 
            Strue[i][j]= dataA1[i][j]/(dataA1[i][j]+dataA2[i][j]);
            Strue_can[i][j]= Strue[i][j];
            Strue_accept[i][j]=0;
            if (i==1) {
                P[j]= P0[j-1];
                Pcan[j]= P[j];
                Paccept[j]= 0;
            }
        }
    }
    
    std::cout << "4" << std::endl;
    
    
    for (i=2;i<=25; i++){
        for (j=1; j<=50;j++){
            
            gridA[i][j]=A[(i-2)*50+j-1];
            gridB[i][j]=B[(i-2)*50+j-1];
            if (j==1) {
                gridA[i][0]= gridA[i][j];
                gridB[i][0]= gridB[i][j];
            }
        }
    }
    
    std::time_t t1, t2;
    t1 = time(NULL); // time 1
    
    std::cout << "5" << std::endl;
    
    
    // ////output////
    // char var_file[1000];
    // sprintf(var_file, "%s/%s", &path[0], &file_index[0]); 
    // FILE *V0 = fopen(var_file, "w");
    // 
    
    // set up results
    int storage_counter = 0;
    int n_saves = (iter/saves) + 1;
    std::vector<std::vector<int> > M_res(n); for(auto a : M_res) a.reserve(n_saves);
    std::vector<std::vector<double> > P_res(k); for(auto a : P_res) a.reserve(n_saves);
    std::vector<std::vector<double> > c_res(1); c_res[0].reserve(n_saves);
    
    // set up a progress bar
    int progress = iter/100;
    std::cout << "|";
    for(int prog = 0; prog < progress ; prog++) std::cout << "-";
    std::cout << "|\n|";
    
    ////MCMC////
    
    std::cout << "6" << std::endl;
    
    //calculate likelihood of initial values
    for (i=1;i<=n;i++){
        for (j=1;j<=k;j++){
            ll[i][j]= logLike(M[i], P[j], dataA1[i][j], dataA2[i][j], Strue[i][j], gridA, gridB, c);
            llcan[i][j]= ll[i][j];
        }
    }
    
    std::cout << "7" << std::endl;
    
    for (i=1;i<=iter;i++){
        if (i%(iter/10)==0) Rprintf("Iter %d out of %d\n", i, iter);
        //update M
        for (j=1; j<=n; j++){
            prime = R::rbinom(1,0.5);
            if (prime==0) prime=-1;
            Mcan[j]= M[j] + prime;
            if ((Mcan[j]<=max_moi) && (Mcan[j]>0)){
                sumcan=0;
                sumori=0;
                for (x=1;x<=k;x++){
                    llcan[j][x]=logLike(Mcan[j], P[x], dataA1[j][x], dataA2[j][x], Strue[j][x], gridA, gridB, c);
                    sumcan+=llcan[j][x];
                    sumori+=ll[j][x];
                }
                //accept
                if (std::log(R::runif(0.0,1.0)) < (sumcan-sumori)) {
                    M[j]=Mcan[j];
                    Maccept[j]++;
                    for (x=1;x<=k;x++){
                        ll[j][x]=llcan[j][x];
                    }
                }
                //reject
                else {
                    Mcan[j]=M[j];
                    for (x=1;x<=k;x++){
                        llcan[j][x]=ll[j][x];
                    }
                }
            }
            else {
                Mcan[j]=M[j];
            }
        }
        //update PP
        for (j=1; j<=k; j++){
            Pcan[j]= R::rnorm(P[j],varP);
            if ((Pcan[j]<1) && (Pcan[j]>0)){
                sumcan=0;
                sumori=0;
                for (y=1;y<=n;y++){
                    llcan[y][j]=logLike(M[y], Pcan[j], dataA1[y][j], dataA2[y][j], Strue[y][j], gridA, gridB, c);
                    sumcan+=llcan[y][j];
                    sumori+=ll[y][j];
                }
                //accept
                if (std::log(R::runif(0.0,1.0)) < (sumcan-sumori)) {
                    P[j]=Pcan[j];
                    Paccept[j]++;
                    for (y=1;y<=n;y++){
                        ll[y][j]=llcan[y][j];
                    }
                }
                //reject
                else {
                    Pcan[j]=P[j];
                    for (y=1;y<=n;y++){
                        llcan[y][j]=ll[y][j];
                    }
                }
            }
            else {
                Pcan[j]=P[j];
            }
        }
        
        if (c>0){
            //update St
            for (j=1;j<=n;j++){
                for (x=1;x<=k;x++){
                    Strue_can[j][x]= R::rnorm(Strue[j][x],varS);
                    q1= R::dnorm(Strue[j][x],Strue_can[j][x],varS,0);
                    q2= R::dnorm(Strue_can[j][x],Strue[j][x],varS,0);
                    if (Strue_can[j][x]<0) {
                        q2= R::pnorm(0, Strue[j][x], varS,1,0);
                        Strue_can[j][x]=0;
                        q1= R::dnorm(Strue[j][x],Strue_can[j][x],varS,0);
                    }
                    else if (Strue_can[j][x]>1) {
                        q2= 1- R::pnorm(1, Strue[j][x], varS,1,0);
                        Strue_can[j][x]=1;
                        q1= R::dnorm(Strue[j][x],Strue_can[j][x],varS,0);			
                    }
                    if (Strue[j][x]==0){
                        q1= R::pnorm(0, Strue_can[j][x], varS,1,0);
                    }
                    else if (Strue[j][x]==1){
                        q1= 1- R::pnorm(1, Strue_can[j][x], varS,1,0);
                    }
                    
                    llcan[j][x]=logLike(M[j], P[x], dataA1[j][x], dataA2[j][x], Strue_can[j][x], gridA, gridB, c);
                    //accept
                    if ((llcan[j][x]-ll[j][x]+std::log(q1)-std::log(q2))>std::log(R::runif(0.0,1.0))){
                        Strue[j][x]= Strue_can[j][x];
                        ll[j][x]= llcan[j][x];
                        Strue_accept[j][x]++;
                    }
                    //reject
                    else{
                        Strue_can[j][x]= Strue[j][x];
                        llcan[j][x]= ll[j][x];
                    }
                } 
            }
        }
        
        //update c
        if (err_method==2){
            c = R::runif(0.0,upper_bound_c);
        }
        if (err_method==3){
            c_can= R::rnorm(c,varE);
            if (c_can>=0){
                //calculate likelihood
                sumcan=0;
                sumori=0;
                for (y=1;y<=n;y++){
                    for (x=1;x<=k;x++){
                        llcan[y][x]=logLike(M[y], P[x], dataA1[y][x], dataA2[y][x], Strue[y][x], gridA, gridB, c_can);
                        sumcan+=llcan[y][x];
                        sumori+=ll[y][x];
                    }
                }
                //accept
                if (std::log(R::runif(0.0,1.0)) < (sumcan-sumori)) {
                    c= c_can;
                    c_accept++;
                    for (y=1;y<=n;y++){
                        for (x=1;x<=k;x++){
                            ll[y][x]=llcan[y][x];
                        }
                    }
                }
                //reject
                else {
                    c_can= c;
                    for (y=1;y<=n;y++){
                        for (x=1;x<=k;x++){
                            llcan[y][x]=ll[y][x];
                        }
                    }
                }
            }
            else{
                //reject
                c_can= c;
            } //end update c
        }
        
        if(i % saves == 0 | i == 1)
        {

        for (x=1;x<=n;x++) M_res[x][storage_counter] = M[x];
        for (x=1;x<=k;x++) P_res[x][storage_counter] = P[x];
        if (err_method==3) c_res[0][storage_counter] = c;
        if(i % progress == 0) std::cout << "-";
        storage_counter++;
        }
    }
    std::cout << "|\n";
    
    t2 = time(NULL); // time 2
    
    // create list results
    Rcpp::List res = Rcpp::List::create(Rcpp::Named("M") = M_res,
                                        Rcpp::Named("P") = P_res,
                                        Rcpp::Named("c") = c_res,
                                        Rcpp::Named("M_accept") = Maccept,
                                        Rcpp::Named("P_accept") = Paccept,
                                        Rcpp::Named("c_accept") = c_accept,
                                        Rcpp::Named("Time") = difftime(t2, t1));
    
    return(res);
}
