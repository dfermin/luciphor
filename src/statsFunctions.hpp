
#ifndef STATISTICSFUNCTIONS_H_
#define STATISTICSFUNCTIONS_H_


#include <list>
#include <map>
using namespace std;

double getMean(list<double> *x);

double getVar(list<double> *x);
double getVar_muZero(list<double> *x);
double getVar_givenMu(list<double> *x, double mu);

double getMode(list<double> *x);

double log_gaussianProb(double mu, double sigma2, double x);
double log_mix_gaussianProb(double mu, double pi, double var1, double var2, double x);

double log_laplaceProb(double mu, double sigma2, double x);
double log_mix_laplaceProb(double mu, double pi, double var1, double var2, double x);
double log_uniformProb(double a, double b);

double runif();
double randNorm(double mu, double sigma);
double getListMinMax(list<double> *Lptr, int getWhat);
double getInterQuartileVariance(list<double> *Lptr);

void pruneList(list<double> *ptr, double trimPercent);

double cum_binomial_prob(int N, int K, double pr);
double binomial_prob(int N, int K, double pr);

#endif /* STATISTICSFUNCTIONS_H_ */

