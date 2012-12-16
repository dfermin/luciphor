/*
 * statsFunctions.cpp
 *
 *  Created on: Apr 16, 2011
 *      Author: dfermin
 */


#include <list>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <boost/filesystem/operations.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "globals.hpp"
#include "statsFunctions.hpp"

using namespace std;

const double PI = 4.0 * atan(1);


// Function returns a random number between 0 - 1 that is from a uniform distribution
double runif() {
	return ( rand() / double(RAND_MAX) );
}



// Function to return the mean of a given list of values
double getMean(list<double> *x) {
	double ret = 0.0;
	double N = (signed) x->size();
	list<double>::iterator l_iter;

	if(N == 0) return 0;

	for(l_iter = x->begin(); l_iter != x->end(); l_iter++)  ret += *l_iter;

	ret /= N;
	return ret;
}


// Function to return the variance of a given list of values
double getVar(list<double> *x) {
	double ret = 0.0;
	double mu = 0.0;
	double N = (signed) x->size();
	double v = 0.0;
	list<double>::iterator L;

	if(g_IS_HCD) mu = getMode(x);
	else mu = getMean(x);

	for(L = x->begin(); L != x->end(); L++) {
		v += pow( (*L - mu), 2.0 );
	}
	N -= 1.0;
	ret = v / N;
	return ret;
}



// Function to return the variance of a given list of values. In this function
// the mean is provided, not calculated
double getVar_givenMu(list<double> *x, double mu) {
	double ret = 0.0;
	double N = (signed) x->size();
	double v = 0.0;
	list<double>::iterator L;

	for(L = x->begin(); L != x->end(); L++) {
		v += pow( (*L - mu), 2.0 );
	}
	N -= 1.0;
	ret = v / N;
	return ret;
}



// Function returns the variance of the given list but assumes the mean of the
// data is zero
double getVar_muZero(list<double> *x) {
	double ret = 0.0;
	double mu = 0.0;
	double N = (signed) x->size();
	double v = 0.0;
	double i = 0.0;
	list<double>::iterator L;

	for(L = x->begin(); L != x->end(); L++) {
		i = dbl_isnan(*L) ? 0 : *L;
		v += pow( (i - mu), 2.0 );
	}
	N -= 1.0;
	ret = v / N;
	return ret;
}



// Function computes the mode for the given list of values
double getMode(list<double> *x) {
	double mode = 0;
	double j, max_j, max_id, a, b;
	double Nbins = 2000;
	double binWidth = 0.0001;
	double LIMIT = 0.1;
	vector<double> v;
	list<double>::iterator L;

	v.resize(Nbins, 0); // initialize vector to have 200 bins, all at zero

	for(L = x->begin(); L != x->end(); L++) {

		for(j = 0; j < (Nbins - 1); j++) {
			a = -LIMIT + (j * binWidth);
			b = a + binWidth;

			if( (*L >= a) & (*L < b) ) {
				v[j]++;
				break; // done, get out of this iteration
			}
		}
	}

	// Now find the largest value in 'v'
	max_j = 0;
	for(int i = 0; i < Nbins; i++) if( v[i] > max_j ) max_j = v[i];

	// now find the element of 'v' that holds the max value
	max_id = 0;
	for(int i = 0; i < Nbins; i++) if( v[i] == max_j ) max_id = i;


	mode = -LIMIT + (max_id*binWidth) + (binWidth/2);

	return mode;
}




// Function returns the probability of a given peak intensity (in log scale)
// under the gaussian pdf with given mu and sigma^2
double log_gaussianProb(double mu, double sigma2, double x) {

	double ret = 0.0;
	ret = - 0.5 * pow(x-mu,2.0) / sigma2 - 0.5 * log(2.00 * PI * sigma2);
	return ret;
}




// Function computes the probability of a given HCD peak distance (in log scale)
// under the Laplace pdf with given mu and sigma^2
double log_laplaceProb(double mu, double sigma2, double x) {

	double ret = 0.0;

	double b = sqrt( (sigma2/2.0) );

	double x1 = -log( (2.0*b) );

	double x2 = -1 * ( (fabs(x-mu))/ b );

	ret = x1 + x2;

	return ret;
}



// Function computes a probability of a given HCD peak distance (in log scale)
// from a *mixture* of 2 laplace distributions.
double log_mix_laplaceProb(double mu, double pi, double var1, double var2, double x) {

	// mu is the same for both distributions
	// pi is the proportion for one of the two distributions
	// var1 is the variance for distro #1 (inter-quartile range)
	// var2 is the variance for distro #2

	double ret = 0.0;
	double x1, x2, x_max;

	x1 = log( pi ) + log_laplaceProb(mu, var1, x);
	x2 = log( 1 - pi ) + log_laplaceProb(mu, var2, x);

	x_max = fmax(x1, x2);

	x1 -= x_max;
	x2 -= x_max;

	ret = x_max + log( exp(x1) + exp(x2) );

	return ret;
}


// Function computes a probability of a given HCD peak distance in log scale
// from a *mixture* of 2 gaussian distributions
double log_mix_gaussianProb(double mu, double pi, double var1, double var2, double x) {

	// mu is the same for both distributions
	// pi is the proportion for one of the two distributions
	// var1 is the variance for distro #1
	// var2 is the variance for distro #2

	double ret = 0.0;
	double x1, x2, x_max;

	x1 = log( pi ) + log_gaussianProb(mu, var1, x);
	x2 = log( 1 - pi ) + log_gaussianProb(mu, var2, x);

	x_max = fmax(x1, x2);

	x1 -= x_max;
	x2 -= x_max;

	ret = x_max + log( exp(x1) + exp(x2) );

	return ret;
}




// Function computes the probability of a given HCD peak distance (in log scale)
// under the uniform pdf. The uniform distribution is assumed for unmatched peaks
double log_uniformProb(double a, double b) {
	// b > a is a requirement

	return -log( b - a );
}




// Function to trim the top/bottom 10% of the passed list
void pruneList(list<double> *ptr, double trimPercent) {
	int i, x, N, fromFront, fromBack;
	list<double>::iterator L;
	list<double> new_L;
	vector<double> v;

	double alpha = trimPercent * 0.5; // we want to take 1/2 of this percentage from top & bottom

	N = (signed) ptr->size();
	x = (int) ((double)N * alpha);

	fromFront = x;
	fromBack = N - x - 1;

	ptr->sort(); // sorted from low to high

	for(L = ptr->begin(); L != ptr->end(); L++) v.push_back(*L);

	for(i = fromFront; i < fromBack; i++) { new_L.push_back( v[i] ); }
	*ptr = new_L;
}




// Function uses the 'box_muller' polar method to generate a random number
// that follows a guassian (Normal) distribution. The code for this was taken
// from: http://www.taygeta.com/random/gaussian.html
double randNorm(double mu, double sigma) {

    double x1, x2, w, y1, t, ret;
    static double y2;
    static int use_last = 0;

    if(use_last) { // use the value of y2 created in last call
		y1 = y2;
		use_last = 0;
    }
    else {
		do {
			x1 = 2.0 * runif() - 1.0;
			x2 = 2.0 * runif() - 1.0;
			w = (x1 * x1) + (x2 * x2);
		} while( w >= 1.0 );

		t = sqrt( (-2.0 * log(w)) / w );
		y1 = x1 * t;
		y2 = x2 * t;
    }

    ret = (mu + y1 * sigma );

    return ret;
}




// Function returns the cumulative binomial probabilty for N trials and K successes
double cum_binomial_prob(int N, int K, double pr) {
  int k;
  double result;
  result = 0.0;

  if(K > 0) {
	  for(k=K;k<=N;k++) {
		result += binomial_prob(N, k, pr);
	  }
  }

  return result;
}



// Function returns the exact binomial probability for N trials and K successes
double binomial_prob(int N, int K, double pr) {
  int i;
  double num, denom, result;

  if(pr <= 0.0 || pr >= 1.0) {
    result = 0.0;
    return result;
  }
  else {
    num = 0.0;
    for(i=1;i<=N;i++) num += log( ((double) i) );

    denom = 0.0;
    for(i=1;i<=K;i++) denom += log( ((double) i) );
    for(i=1;i<=(N-K);i++) denom += log( ((double) i) );

    result = num - denom + ((double) K) * log(pr) + ((double) (N-K)) * log(1.0-pr);
    result = exp(result);
  }

  return result;
}



// Function returns the minimum or maximum value in a list of doubles
double getListMinMax(list<double> *Lptr, int getWhat) {
	double ret = 0.0;
	list<double>::iterator L;
	const int MAX = 1;

	Lptr->sort(); // by default, order low-to-high

	if(getWhat == MAX) Lptr->reverse(); // order from high-to-low

	ret = Lptr->front();
	return ret;
}



// Function to get the interquartile range variance from the given list reference
double getInterQuartileVariance(list<double> *Lptr) {
	double ret = 0.0;
	list<double>::iterator L;
	list<double> g;
	vector<double> v;
	int qLow, qHi, p;
	double quantile1, quantile2;
	int N = (signed) Lptr->size();
	double i;

	quantile1 = 0.1;
	quantile2 = 0.9;

	v.reserve( Lptr->size() );

	Lptr->sort(); // sorted low-to-high

	for(L = Lptr->begin(); L != Lptr->end(); L++) v.push_back(*L);

	i = quantile1 * (double)N;
	qLow = (int) round_dbl(i, 0);

	i = quantile2 * (double)N;
	qHi = (int) round_dbl(i, 0);


	for(int p = qLow; p < qHi; p++) g.push_back( v.at(p) );

	ret = getVar(&g); // get variance

	return ret;
}

