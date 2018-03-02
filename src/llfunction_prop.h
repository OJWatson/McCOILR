#include <stdio.h>
#include <Rcpp.h>
#include <string>
using namespace std;

double logLike(int M, double P, double dataA1, double dataA2, double Strue,
               std::vector<std::vector<double> > gridA, std::vector<std::vector<double> > gridB, 
               double c);


