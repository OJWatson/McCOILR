#include "llfunction_prop.h"

double logLike(int M, double P, double dataA1, double dataA2, double Strue, double gridA[][51], double gridB[][51], double c){

	double llobs = 0;
	double lltrue = 0.0;
	double alpha = 0.0, beta = 0.0;
	int order_P, k, l = 0;
		
	if ((dataA1<0)||(dataA2<0)){
		return 0;
	}
	else{
	//llobs
		if (c==0) llobs = 0;
		else{
			double So= dataA1/(dataA1+dataA2);
			double intensity= std::sqrt(pow(dataA1,2) + std::pow(dataA2,2));
			double sd= std::sqrt(c/intensity);
			if (So==0){
				llobs = R::pnorm(0, Strue, sd, 1, 0);
			}
			else if (So==1){
				llobs = 1- R::pnorm(1, Strue, sd, 1, 0);
			}
			else {
				llobs = R::dnorm(So, Strue, sd, 0);
			}
			llobs = std::log(llobs);
		}

		//lltrue
		if (M==1){
			if (Strue ==1) lltrue= std::log(P);
			else if (Strue ==0) lltrue= std::log(1-P);
			else lltrue=-9999999;
		} 
		else if (Strue ==1){
			lltrue= M* std::log(P);
		}
		else if (Strue ==0){
			lltrue= M* std::log(1-P);
		}
		else{
			if (P<=0.5){
				order_P = round(P/0.01);
				if (order_P==0) order_P=1;  
				alpha= gridA[M][order_P];
				beta= gridB[M][order_P];
				lltrue= R::dbeta(Strue, alpha, beta, 1) + std::log(1-(double)std::pow(P,M)-(double)std::pow((1-P),M));
			}
			else{
				P=1-P;
				Strue=1-Strue;
				order_P = round(P/0.01);	
				if (order_P==0) order_P=1;
				alpha= gridA[M][order_P];
				beta= gridB[M][order_P];
				lltrue= R::dbeta(Strue, alpha, beta, 1) + std::log(1-(double)std::pow(P,M)-(double)std::pow((1-P),M));
			}				
		}
		return (lltrue+llobs);
	}

}

