#include "llfunction_cat.h"

double logLike_het(int M, double P, double S2, double e1, double e2){
	double lltrue = 0.0;
	if (S2 ==1){
		lltrue= std::log( (1-e1)*((double)std::pow(P,M)) + e2/2*(1- (double)std::pow(P,M)- (double)std::pow((1-P),M)) );
	}
	else if (S2 ==0){
		lltrue= std::log( (1-e1)*((double)std::pow((1-P),M)) + e2/2*(1- (double)std::pow(P,M)- (double)std::pow((1-P),M)) );
	}
	else if (S2 ==0.5){
		if ((M==1) && (e1==0)) lltrue=-9999999;
		else lltrue= std::log( e1*((double)std::pow(P,M)) + e1*((double)std::pow((1-P),M)) + (1-e2)*(1- (double)std::pow(P,M)- (double)std::pow((1-P),M)) );
	}
	else lltrue=0.0;
	return (lltrue);

}

