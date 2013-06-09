/*
 * MSProductClass.cpp
 *
 *  Created on: Apr 16, 2011
 *      Author: dfermin
 */


#include <algorithm>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <deque>
#include "structs.hpp"
#include "nonParamDist.hpp"
#include "statsFunctions.hpp"

using namespace std;
using std::ofstream;

const double PI = 4*atan(1);


// HW
double getMean2(deque<double> *x) {
        double ret = 0.0;
        double N = (signed) x->size();
        deque<double>::iterator iter;

        for(iter = x->begin(); iter != x->end(); iter++)
        	ret += *iter;

        ret /= N;
        return ret;
}


double getVariance2(deque<double> *ptr) {
    double ret = 0.0;
    double mu = 0.0;
    double N = (signed) ptr->size();
    double v = 0.0;
    deque<double>::iterator iter;
    double D = 0.0;

    mu = getMean2(ptr);

    for(iter = ptr->begin(); iter != ptr->end(); iter++) {
        D = *iter;
    	v += pow( (D - mu), 2.0 );
    }

    N -= 1.0;
    ret = v / N;
    return ret;
}





double getLogNPdensityInt_b(double x, modelParamStruct *paramPtr) {

	double a = 0, b = 0;
	double sum = 0;
	double tmp1 = 0, tmp2 = 0;
	double fx = 0;
	int ntick = paramPtr->ntick_int;

	double xx = x; 
 	double start_tick = paramPtr->tickMarks_int_b.at(0);
 	double end_tick = paramPtr->tickMarks_int_b.at(ntick-1);
 	double start_val = paramPtr->f_int_b.at(0);
 	double end_val = paramPtr->f_int_b.at(ntick-1);

	// we are trying to figure out which two bins encompass the value of 'x'
	for(int j = (ntick-1); j >= 1; j--) { // working backwards here
		int i = j - 1;
		a = paramPtr->tickMarks_int_b.at(i);
		b = paramPtr->tickMarks_int_b.at(j);

		if(xx >= a && xx < b) {
			// We have reached the bin that contains 'x'
			// Now compute the area under the curve upto the point that includes 'x'
			tmp1 = (b-xx) / (b-a);
			tmp2 = (xx-a) / (b-a);

			fx = (tmp1 * paramPtr->f_int_b.at(i)) + (tmp2 * paramPtr->f_int_b.at(j));
			sum = fx;
			break; // leave loop since you are done
		}
	}

	if(xx <= start_tick) sum = start_val;
	else if(xx >= end_tick) sum = end_val;
	else {}	
	sum = log(sum);
	return sum;
}




double getLogNPdensityInt_y(double x, modelParamStruct *paramPtr) {

	double a = 0, b = 0;
	double sum = 0;
	double tmp1 = 0, tmp2 = 0;
	double fx = 0;
	int ntick = paramPtr->ntick_int;

	double xx = x; 
 	double start_tick = paramPtr->tickMarks_int_y.at(0);
 	double end_tick = paramPtr->tickMarks_int_y.at(ntick-1);
 	double start_val = paramPtr->f_int_y.at(0);
 	double end_val = paramPtr->f_int_y.at(ntick-1);

	// we are trying to figure out which two bins encompass the value of 'x'
	for(int j = (ntick-1); j >= 1; j--) { // working backwards here
		int i = j - 1;
		a = paramPtr->tickMarks_int_y.at(i);
		b = paramPtr->tickMarks_int_y.at(j);

		if(xx >= a && xx < b) {
			// We have reached the bin that contains 'x'
			// Now compute the area under the curve upto the point that includes 'x'
			tmp1 = (b-xx) / (b-a);
			tmp2 = (xx-a) / (b-a);

			fx = (tmp1 * paramPtr->f_int_y.at(i)) + (tmp2 * paramPtr->f_int_y.at(j));
			sum = fx;
			break; // leave loop since you are done
		}
	}

	if(xx <= start_tick) sum = start_val;
	else if(xx >= end_tick) sum = end_val;
	else {}	
	sum = log(sum);
	return sum;
}


double getLogNPdensityInt_U(double x, modelParamStruct *paramPtr) {

	double a = 0, b = 0;
	double sum = 0;
	double tmp1 = 0, tmp2 = 0;
	double fx = 0;
	int ntick = paramPtr->ntick_int;

	double xx = x; 
 	double start_tick = paramPtr->tickMarks_int_U.at(0);
 	double end_tick = paramPtr->tickMarks_int_U.at(ntick-1);
 	double start_val = paramPtr->f_int_U.at(0);
 	double end_val = paramPtr->f_int_U.at(ntick-1);

	// we are trying to figure out which two bins encompass the value of 'x'

	for(int j = (ntick-1); j >= 1; j--) { // working backwards here
		int i = j - 1;
		a = paramPtr->tickMarks_int_U.at(i);
		b = paramPtr->tickMarks_int_U.at(j);

		if(xx >= a && xx < b) {
			// We have reached the bin that contains 'x'
			// Now compute the area under the curve upto the point that includes 'x'
			tmp1 = (b-xx) / (b-a);
			tmp2 = (xx-a) / (b-a);

			fx = (tmp1 * paramPtr->f_int_U.at(i)) + (tmp2 * paramPtr->f_int_U.at(j));
			sum = fx;
			break; // leave loop since you are done
		}
	}
	
	if(xx <= start_tick) sum = start_val;
	else if(xx >= end_tick) sum = end_val;
	else {}	

	sum = log(sum);
	return sum;
}




double getLogNPdensityDist(double x, modelParamStruct *paramPtr) {

	double a = 0, b = 0;
	double sum = 0;
	double tmp1 = 0, tmp2 = 0;
	double fx = 0;
	int ntick = paramPtr->ntick_dist;
	double xx = x; // mz_adjust_dist(x);
 	double start_tick = paramPtr->tickMarks_dist.at(0);
 	double end_tick = paramPtr->tickMarks_dist.at(ntick-1);
 	double start_val = paramPtr->f_dist.at(0);
 	double end_val = paramPtr->f_dist.at(ntick-1);

	// we are trying to figure out which two bins encompass the value of 'x'
	for(int j = (ntick-1); j >= 1; j--) { // working backwards here
		int i = j - 1;
		a = paramPtr->tickMarks_dist.at(i);
		b = paramPtr->tickMarks_dist.at(j);

		if(xx >= a && xx < b) {
			// We have reached the bin that contains 'x'
			// Now compute the area under the curve upto the point that includes 'x'
			tmp1 = (b-xx) / (b-a);
			tmp2 = (xx-a) / (b-a);

			fx = (tmp1 * paramPtr->f_dist.at(i)) + (tmp2 * paramPtr->f_dist.at(j));
			sum = fx;
			break; // leave loop since you are done
		}
	}

	if(xx <= start_tick) sum = start_val;
	else if(xx >= end_tick) sum = end_val;
	else {}	
	sum = log(sum);
	return sum;
}



double getLogNPdensityDist_U(double x, modelParamStruct *paramPtr) {

	double a = 0, b = 0;
	double sum = 0;
	double tmp1 = 0, tmp2 = 0;
	double fx = 0;
	int ntick = paramPtr->ntick_dist;
	double xx = x; // mz_adjust_dist(x);
 	double start_tick = paramPtr->tickMarks_dist_U.at(0);
 	double end_tick = paramPtr->tickMarks_dist_U.at(ntick-1);
 	double start_val = paramPtr->f_dist_U.at(0);
 	double end_val = paramPtr->f_dist_U.at(ntick-1);

	// we are trying to figure out which two bins encompass the value of 'x'

	for(int j = (ntick-1); j >= 1; j--) { // working backwards here
		int i = j - 1;
		a = paramPtr->tickMarks_dist_U.at(i);
		b = paramPtr->tickMarks_dist_U.at(j);

		if(xx >= a && xx < b) {
			// We have reached the bin that contains 'x'
			// Now compute the area under the curve upto the point that includes 'x'
			tmp1 = (b-xx) / (b-a);
			tmp2 = (xx-a) / (b-a);

			fx = (tmp1 * paramPtr->f_dist_U.at(i)) + (tmp2 * paramPtr->f_dist_U.at(j));
			sum = fx;
			break; // leave loop since you are done
		}
	}
	
	if(xx <= start_tick) sum = start_val;
	else if(xx >= end_tick) sum = end_val;
	else {}	

	sum = log(sum);
	return sum;
}



double normalDensity(double curTickMark, double curScore, double h) {
   double res = 0;
   double x = (curTickMark - curScore) / h;

   res =  1.0 / sqrt(2.0 * PI) * exp( -0.5 * x * x );
   return res;
}



void percentileTrim(deque<double> *ptr, double percent) {
	deque<double> ret;
	int I = 0;
	int J = 0;
	double N = (signed) ptr->size();
	double half = percent * 0.5;
	double x = round_dbl((half * N), 0);

	I = (int) x;
	J = (int) (N - x);

	for(int i = I; i < J; i++) ret.push_back( ptr->at(i) );

	ptr->clear();
	*ptr = ret;
}





// Function "trims" the m/z distances of the passed list of m/z distances
// so that extreme outlier values are truncated to a given value.
// This was Hyungwon's idea for controlling for extreme outliers.
double mz_adjust_dist(double x) {
	double res = x;
	if(g_IS_HCD) {
		if(res >= 0) {
			while(fabs(res) > 0.5) res -= 1;
		}
		else {
			while(fabs(res) > 0.5) res += 1;
		}
	}	

	return res;
}




void estimateNonparamDist(list<double> *ptr, modelParamStruct *paramPtr, double mz_err) {

	deque<double> *X = new deque<double>;
	deque<double> *tmpf = new deque<double>;
	deque<double>::iterator tic;
	deque<double>::iterator curPeak;

	// This is a more efficient way of converting the data from a list to a deque
	int N = (signed)ptr->size();
	X->resize( N );
	std::copy(ptr->begin(), ptr->end(), X->begin());

	int ntick = 2000; // interval is from -mz_err to mz_err
	double nmatched = 0;
	paramPtr->ntick_dist = ntick;
	double min_dist = -0.5;
	double max_dist = 0.5;
	paramPtr->tickMarks_dist.resize(paramPtr->ntick_dist);
	int ii = 0;
	for(tic = paramPtr->tickMarks_dist.begin(); tic != paramPtr->tickMarks_dist.end(); tic++) {
		paramPtr->tickMarks_dist.at(ii) = min_dist + ((double) ii) * (max_dist - min_dist) / ((double) (paramPtr->ntick_dist-1));
		ii++;
	}

	nmatched = (signed) X->size();

	double sigma = sqrt( getVar(ptr) );
	paramPtr->bw_dist = 1.06 * (sigma / pow(nmatched, 0.2) );
	paramPtr->bw_dist *= 0.1;
	cerr << "Distance Matched NP Bandwidth: " << paramPtr->bw_dist << "\t(N = " << N << ")\n";

	deque<double>::iterator curScore;
	paramPtr->f_dist.resize(paramPtr->ntick_dist);

	double tmp_sum = 0;

	ii = 0;
	for(tic = paramPtr->tickMarks_dist.begin(); tic != paramPtr->tickMarks_dist.end(); tic++) {
		tmp_sum = 0;
		// iterate over real observations
		for(curScore = X->begin(); curScore != X->end(); curScore++) {
			tmp_sum += normalDensity( *tic, *curScore, paramPtr->bw_dist );
		}
		tmp_sum /= ( ((double) nmatched) * paramPtr->bw_dist ) ;	

		paramPtr->f_dist.at(ii) = tmp_sum > TINY_NUM ? tmp_sum : TINY_NUM;
		ii++;
	}

	delete(X); X = NULL;
	ptr->clear(); // free up memory
}



void estimateNonparamDist_U(list<double> *ptr, modelParamStruct *paramPtr, double mz_err) {

	deque<double> *X = new deque<double>;
	deque<double>::iterator tic;
	deque<double>::iterator curPeak;
	int jj;

	//This is a more efficient way of converting the data from a list to a deque
	int N = (signed)ptr->size();
	X->resize( N );
	std::copy(ptr->begin(), ptr->end(), X->begin());

	int ntick = 2000; // interval is from -mz_err to mz_err
	double nmatched = 0;
	paramPtr->ntick_dist = ntick;
	double min_dist = -0.5;
	double max_dist = 0.5;
	paramPtr->tickMarks_dist_U.resize(paramPtr->ntick_dist);
	int ii = 0;
	for(tic = paramPtr->tickMarks_dist_U.begin(); tic != paramPtr->tickMarks_dist_U.end(); tic++) {
		paramPtr->tickMarks_dist_U.at(ii) = min_dist + ((double) ii) * (max_dist - min_dist) / ((double) (paramPtr->ntick_dist-1));
		ii++;
	}

	nmatched = (signed) X->size();

	double sigma = sqrt( getVar(ptr) );   // need writing another function for variance!!
	paramPtr->bw_dist_U = 1.06 * (sigma / pow(nmatched, 0.2) );
	cerr << "Distance Unmatched NP Bandwidth: " << paramPtr->bw_dist_U << "\t(N = " << N << ")\n";;


	deque<double>::iterator curScore;
	paramPtr->f_dist_U.resize(paramPtr->ntick_dist);

	double tmp_sum = 0;


	jj = 0;
	/******************************
	// iterate over real observations
	for(curScore = X->begin(); curScore != X->end(); curScore++) {
		jj++;
	}
	nmatched = jj;
	*******************************/
	nmatched = (signed) X->size();

	ii = 0;
	for(tic = paramPtr->tickMarks_dist_U.begin(); tic != paramPtr->tickMarks_dist_U.end(); tic++) {
		tmp_sum = 0;
		// iterate over real observations
		for(curScore = X->begin(); curScore != X->end(); curScore++) {
			tmp_sum += normalDensity( *tic, *curScore, paramPtr->bw_dist_U );
		}
		tmp_sum /= ( ((double) nmatched) * paramPtr->bw_dist_U ) ;	

		paramPtr->f_dist_U.at(ii) = tmp_sum > TINY_NUM ? tmp_sum : TINY_NUM;
		paramPtr->f_dist_U.at(ii) = TINY_DEN;
		ii++;
	}

	delete(X); X = NULL;
	ptr->clear(); // free up memory
}




/************************** Intensity Density Estimation ************************/



void estimateNonparamInt_b(list<double> *ptr, modelParamStruct *paramPtr, double mz_err) {

	deque<double> *X = new deque<double>;
	deque<double>::iterator tic;
	deque<double>::iterator curPeak;
	list<double>::iterator L;


	// This is a more efficient way of converting the data from a list to a deque
	int N = (signed)ptr->size();
	X->resize( N );
	std::copy(ptr->begin(), ptr->end(), X->begin());


	//for(L = ptr->begin(); L != ptr->end(); L++) X->push_back(*L);

	int ntick = 2000; // interval is from -mz_err to mz_err
	double nmatched = 0;
	paramPtr->ntick_int = ntick;

	L = ptr->begin();
	double min_int = *L;
	double max_int = *L;
	for(L = ptr->begin(); L != ptr->end(); L++) {
		if(*L > max_int) max_int = *L; 
		else if(*L < min_int) min_int = *L; 
		else { }
	}	

	min_int = -3.0;
	max_int = 8.0;

	paramPtr->tickMarks_int_b.resize(paramPtr->ntick_int);
	int ii = 0;
	for(tic = paramPtr->tickMarks_int_b.begin(); tic != paramPtr->tickMarks_int_b.end(); tic++) {
		paramPtr->tickMarks_int_b.at(ii) = min_int + ((double) ii) * (max_int - min_int) / ((double) (paramPtr->ntick_int-1));
		ii++;
	}

	nmatched = X->size();

	double sigma = sqrt( getVar(ptr) );   // need writing another function for variance!!
	paramPtr->bw_int_b = 1.06 * (sigma / pow(nmatched, 0.2) );
	paramPtr->bw_int_b *= 0.5;

	cerr << "b-ion Intensity Matched NP Bandwidth: " << paramPtr->bw_int_b << "\t(N = " << N << ")\n";

	deque<double>::iterator curScore;
	paramPtr->f_int_b.resize(paramPtr->ntick_int);

	double tmp_sum = 0;

	ii = 0;
	for(tic = paramPtr->tickMarks_int_b.begin(); tic != paramPtr->tickMarks_int_b.end(); tic++) {
		tmp_sum = 0;
		// iterate over real observations
		for(curScore = X->begin(); curScore != X->end(); curScore++) {
			tmp_sum += normalDensity( *tic, *curScore, paramPtr->bw_int_b );
		}
		tmp_sum /= ( ((double) nmatched) * paramPtr->bw_int_b ) ;	

		paramPtr->f_int_b.at(ii) = tmp_sum > TINY_NUM ? tmp_sum : TINY_NUM;
		ii++;
	}
}



void estimateNonparamInt_y(list<double> *ptr, modelParamStruct *paramPtr, double mz_err) {

	deque<double> *X = new deque<double>;
	deque<double>::iterator tic;
	deque<double>::iterator curPeak;
	list<double>::iterator L;

	// This is a more efficient way of converting the data from a list to a deque
	int N = (signed)ptr->size();
	X->resize( N );
	std::copy(ptr->begin(), ptr->end(), X->begin());

	//for(L = ptr->begin(); L != ptr->end(); L++) X->push_back(*L);

	int ntick = 2000; // interval is from -mz_err to mz_err
	double nmatched = 0;
	paramPtr->ntick_int = ntick;

	L = ptr->begin();
	double min_int = *L;
	double max_int = *L;
	for(L = ptr->begin(); L != ptr->end(); L++) {
		if(*L > max_int) max_int = *L; 
		else if(*L < min_int) min_int = *L; 
		else { }
	}	

	min_int = -3.0;
	max_int = 8.0;

	paramPtr->tickMarks_int_y.resize(paramPtr->ntick_int);
	int ii = 0;
	for(tic = paramPtr->tickMarks_int_y.begin(); tic != paramPtr->tickMarks_int_y.end(); tic++) {
		paramPtr->tickMarks_int_y.at(ii) = min_int + ((double) ii) * (max_int - min_int) / ((double) (paramPtr->ntick_int-1));
		ii++;
	}

	nmatched = X->size();

	double sigma = sqrt( getVar(ptr) );   // need writing another function for variance!!
	paramPtr->bw_int_y = 1.06 * (sigma / pow(nmatched, 0.2) );
	paramPtr->bw_int_y *= 0.5;

	cerr << "y-ion Intensity Matched NP Bandwidth: " << paramPtr->bw_int_y << "\t(N = " << N << ")\n";;

	deque<double>::iterator curScore;
	paramPtr->f_int_y.resize(paramPtr->ntick_int);

	double tmp_sum = 0;

	ii = 0;
	for(tic = paramPtr->tickMarks_int_y.begin(); tic != paramPtr->tickMarks_int_y.end(); tic++) {
		tmp_sum = 0;
		// iterate over real observations
		for(curScore = X->begin(); curScore != X->end(); curScore++) {
			tmp_sum += normalDensity( *tic, *curScore, paramPtr->bw_int_y );
		}
		tmp_sum /= ( ((double) nmatched) * paramPtr->bw_int_y ) ;	

		paramPtr->f_int_y.at(ii) = tmp_sum > TINY_NUM ? tmp_sum : TINY_NUM;
		ii++;
	}
}


void estimateNonparamInt_U(list<double> *ptr, modelParamStruct *paramPtr, double mz_err) {

	deque<double> *X = new deque<double>;
	deque<double>::iterator tic;
	deque<double>::iterator curPeak;
	list<double>::iterator L;
	
	int jj;

	// This is a more efficient way of converting the data from a list to a deque
	int N = (signed)ptr->size();
	X->resize( N );
	std::copy(ptr->begin(), ptr->end(), X->begin());

//	for(L = ptr->begin(); L != ptr->end(); L++) X->push_back(*L);



	int ntick = 2000; // interval is from -mz_err to mz_err
	double nmatched = 0;
	paramPtr->ntick_int = ntick;
	// double min_dist = -mz_err - 0.001;
	// double max_dist = mz_err + 0.001;

	L = ptr->begin();
	double min_int = *L;
	double max_int = *L;
	for(L = ptr->begin(); L != ptr->end(); L++) {
		if(*L > max_int) max_int = *L; 
		else if(*L < min_int) min_int = *L; 
		else { }
	}	

	min_int = -3.0;
	max_int = 8.0;

	paramPtr->tickMarks_int_U.resize(paramPtr->ntick_int);
	int ii = 0;
	for(tic = paramPtr->tickMarks_int_U.begin(); tic != paramPtr->tickMarks_int_U.end(); tic++) {
		paramPtr->tickMarks_int_U.at(ii) = min_int + ((double) ii) * (max_int - min_int) / ((double) (paramPtr->ntick_int-1));
		ii++;
	}

	nmatched = X->size();

	double sigma = sqrt( getVar(ptr) );   // need writing another function for variance!!
	paramPtr->bw_int_U = 1.06 * (sigma / pow(nmatched, 0.2) );
	paramPtr->bw_int_U *= 0.5;

	cerr << "Intensity Unmatched NP Bandwidth: " << paramPtr->bw_int_U << "\t(N = " << N << ")\n";


	deque<double>::iterator curScore;
	paramPtr->f_int_U.resize(paramPtr->ntick_int);

	double tmp_sum = 0;


	jj = 0;
	// iterate over real observations
	for(curScore = X->begin(); curScore != X->end(); curScore++) {
			jj++;
	}
	nmatched = jj;

	ii = 0;
	for(tic = paramPtr->tickMarks_int_U.begin(); tic != paramPtr->tickMarks_int_U.end(); tic++) {
		tmp_sum = 0;
		// iterate over real observations
		for(curScore = X->begin(); curScore != X->end(); curScore++) {
				tmp_sum += normalDensity( *tic, *curScore, paramPtr->bw_int_U );
		}
		tmp_sum /= ( ((double) nmatched) * paramPtr->bw_int_U ) ;	

		paramPtr->f_int_U.at(ii) = tmp_sum > TINY_NUM ? tmp_sum : TINY_NUM;
		ii++;
	}
}




void printDensity(modelParamStruct *paramPtr) {

	deque<double>::iterator tic;
	int ii;

	// distance
	ofstream fd("density_distance.txt", ios::out);
	fd << "tic\tm\tu\n";
	ii = 0;
	for(tic = paramPtr->tickMarks_dist.begin(); tic != paramPtr->tickMarks_dist.end(); tic++) {
		fd << *tic <<  "\t" << paramPtr->f_dist.at(ii) << "\t" <<  TINY_DEN << endl;
		ii++;
	}
	fd.close();			



	/*********
	 *
	 * intensity
	 *
	**********/
	ofstream fi("density_intensity.txt", ios::out);
	fi << "tic\tclass\tintensity\n";
	ii = 0;
	for(tic = paramPtr->tickMarks_int_b.begin(); tic != paramPtr->tickMarks_int_b.end(); tic++) {
		fi << *tic <<  "\tb\t" << paramPtr->f_int_b.at(ii) << endl;
		ii++;
	}

	ii = 0;
	for(tic = paramPtr->tickMarks_int_y.begin(); tic != paramPtr->tickMarks_int_y.end(); tic++) {
		fi << *tic << "\ty\t" << paramPtr->f_int_y.at(ii) << endl;
		ii++;
	}


	ii = 0;
	for(tic = paramPtr->tickMarks_int_U.begin(); tic != paramPtr->tickMarks_int_U.end(); tic++) {
		fi << *tic << "\tU\t" << paramPtr->f_int_U.at(ii) << endl;
		ii++;
	}
	fi.close();


}





