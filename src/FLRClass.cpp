/*
 * fdrEstimations.cpp
 *
 *  Created on: Jun 22, 2012
 *      Author: dfermin
 */



#include <iostream>
#include <deque>
#include <string>
#include <cmath>
#include <list>

#include "globals.hpp"
#include "structs.hpp"
#include "FLRClass.hpp"

using namespace std;

const int NMARKS = 10001; // we want N+1 bins; here N = 500

const double PI = 4*atan(1);

// Default contructor
FLRClass::FLRClass(deque<flrStruct> *ptr, double d) {
	allPSMs.clear();
	tickMarks.clear();
	maxDeltaScore = 0;
	Nreal = 0;
	Ndecoy = 0;
	Nreal2 = 0;
	Ndecoy2 = 0;
	bw_real = 0;
	bw_decoy = 0;



	allPSMs = *ptr;
	maxDeltaScore = d * 1.001; // need maxDeltaScore to be slightly larger for binning;

	deque<flrStruct>::iterator iter;
	for(iter = allPSMs.begin(); iter != allPSMs.end(); iter++) {
		if(iter->deltaScore > 0.1) { 
			if(iter->isDecoy) decoyDeq.push_back( *iter );
			else realDeq.push_back( *iter );
		}

		if(iter->deltaScore > 0.1) {
			if(iter->isDecoy) Ndecoy2++;
			else Nreal2++;
		}
	}

	Nreal = (signed) realDeq.size();
	Ndecoy = (signed) decoyDeq.size();

    // initialize tickMarks deque
	tickMarks.resize(NMARKS);
	tickMarks[0] = 0.0;
	for(int i = 0; i < NMARKS; i++) {
		tickMarks[i] = ((double) i) * maxDeltaScore / (double) NMARKS;
	}

	f0.resize(NMARKS);
	f1.resize(NMARKS);

	localFDR.resize(Nreal);
	globalFDR.resize(Nreal);

	getBandWidth("real");
	getBandWidth("decoy");

}


//Function computes the bandwidth for the given deque of observations
void FLRClass::getBandWidth(string x) {
	deque<flrStruct> *ptr = NULL;
	deque<flrStruct>::iterator iter;
	double sigma = 0;
	double N = 0;
	double ret = 0;


	if(x.compare("real") == 0) ptr = &realDeq;
	else ptr = &decoyDeq;

	N = (signed) ptr->size();

	sigma = sqrt( getVariance(ptr) );

	ret = 1.06 * (sigma / pow(N, 0.2) );

	if(x.compare("real")== 0) ret *= 0.1;
	else ret *= 0.02;

	if(x.compare("real")== 0) bw_real = ret;
	else bw_decoy = ret;
}



// Function to return the mean of a given deque of values
double FLRClass::getMean(deque<flrStruct> *x) {
        double ret = 0.0;
        double N = (signed) x->size();
        deque<flrStruct>::iterator iter;

        for(iter = x->begin(); iter != x->end(); iter++)
        	ret += iter->deltaScore;

        ret /= N;
        return ret;
}



// Function returns the variance of deltaScore in the given deque
double FLRClass::getVariance(deque<flrStruct> *ptr) {
    double ret = 0.0;
    double mu = 0.0;
    double N = (signed) ptr->size();
    double v = 0.0;
    deque<flrStruct>::iterator iter;
    double D = 0.0;

    mu = getMean(ptr);

    for(iter = ptr->begin(); iter != ptr->end(); iter++) {
        D = iter->deltaScore;
    	v += pow( (D - mu), 2.0 );
    }

    N -= 1.0;
    ret = v / N;
    return ret;
}


// Function evaluates normal density function at a specific point
double FLRClass::normalDensity(double curTickMark, double curScore, double h) {
   double res = 0;
   double x = (curTickMark - curScore) / h;

   res =  1.0 / sqrt(2.0 * PI) * exp( -0.5 * x * x );
   return res;
}



// Function evaluates each deltaScore at each tick mark storing the result in
// either f0 (decoys) or f1 (real)
void FLRClass::evalTickMarks() {
	deque<double>::iterator tic;
	deque<flrStruct>::iterator curScore;
	double sum = 0;
	int i = 0;

	/*
	** For real data
	*/
	i = 0;
	for(tic = tickMarks.begin(); tic != tickMarks.end(); tic++) {
		sum = 0;
		// iterate over real observations
		for(curScore = realDeq.begin(); curScore != realDeq.end(); curScore++) {
			sum += normalDensity( *tic, curScore->deltaScore, bw_real );
		}
		sum *= ( 1.0 / ( ((double)Nreal) * bw_real ) );

		f1.at(i++) = sum > TINY_NUM ? sum : TINY_NUM;
	}


	/*
	** For decoy data
	*/
	i = 0;
	for(tic = tickMarks.begin(); tic != tickMarks.end(); tic++) {
		sum = 0;
		// iterate over decoy observations
		for(curScore = decoyDeq.begin(); curScore != decoyDeq.end(); curScore++) {
			sum += normalDensity( *tic, curScore->deltaScore, bw_decoy );
		}
		sum *= ( 1.0 / ( ((double)Ndecoy) * bw_decoy ) );

		f0.at(i++) = sum > TINY_NUM ? sum : TINY_NUM;
	}
}




// Function computes the local and global FDRs
void FLRClass::calcBothFDRs() {
	double tmp_score;
	double AUC_rev_0 = 0; // Area-Under the-Curve from end of tick marks working backwards (f0 data)
	double AUC_rev_1 = 0; // Area-Under the-Curve from end of tick marks working backwards (f1 data)
	double ratio = 0;
	int i = 0;

	deque<flrStruct>::iterator curScore;

	i = 0;
	for(curScore = realDeq.begin(); curScore != realDeq.end(); curScore++) {
		tmp_score = curScore->deltaScore;
		if(tmp_score < 0.1) tmp_score = 0.1;		
		// local FDR
		ratio = 0;
		AUC_rev_0 = getLocalAUC(tmp_score, 0);
		AUC_rev_1 = getLocalAUC(tmp_score, 1);

		//ratio = ( ((double)Ndecoy) / ((double)Nreal) ) * ( AUC_rev_0 / AUC_rev_1 );
		ratio = ( Ndecoy2 / Nreal2 ) * ( AUC_rev_0 / AUC_rev_1 );
		localFDR.at( i ) = ratio;


		// GLOBAL FDR
		ratio = 0;
		AUC_rev_0 = getGlobalAUC(tmp_score, 0);
		AUC_rev_1 = getGlobalAUC(tmp_score, 1);

		//ratio = ( ((double)Ndecoy) / ((double)Nreal) ) * ( AUC_rev_0 / AUC_rev_1 );
		ratio = ( Ndecoy2 / Nreal2 ) * ( AUC_rev_0 / AUC_rev_1 );
		globalFDR.at( i ) = ratio;
		i++;
	}

	raw_globalFDR = globalFDR;
	raw_localFDR = localFDR;
}



// Function computes the area under the curve (AUC) for the bin that contains the passed score
double FLRClass::getLocalAUC(double x, int whichF) {

	deque<double> *Fptr = NULL;

	if(whichF == 0) Fptr = &f0;
	else Fptr = &f1;

	double a = 0, b = 0;
	double sum = 0;
	double tmp1 = 0, tmp2 = 0;
	double fx = 0;

 	double start_tick = tickMarks.at(0);
 	double end_tick = tickMarks.at(NMARKS-1);
 	double start_val = Fptr->at(0);
 	double end_val = Fptr->at(NMARKS-1);

	// we are trying to figure out which two bins encompass the value of 'x'
	for(int j = (NMARKS - 1); j >= 1; j--) { // working backwards here
		int i = j - 1;
		a = tickMarks.at(i);
		b = tickMarks.at(j);

		if(x >= a) {
			// We have reached the bin that contains 'x'
			// Now compute the area under the curve upto the point that includes 'x'
			tmp1 = (b-x) / (b-a);
			tmp2 = (x-a) / (b-a);

			fx = (tmp1 * Fptr->at(i)) + (tmp2 * Fptr->at(j));
			sum += fx;
			break; // leave loop since you are done
		}
	}

	if(x <= start_tick) sum = start_val;
	else if(x >= end_tick) sum = end_val;
	else {}	

	return sum;
}




// Function computes the area under the curve (AUC) for the bin that contains the passed score.
// This function is identical to getLocalAUC except that it computes the area before and after
// the bin contianing x
double FLRClass::getGlobalAUC(double x, int whichF) {

	deque<double> *Fptr = NULL;

	if(whichF == 0) Fptr = &f0;
	else Fptr = &f1;

	double a = 0, b = 0;
	double sum = 0;
	double tmp1 = 0, tmp2 = 0;
	double fx = 0;

	// we are trying to figure out which two bins encompass the value of 'x'
	sum = 0;
	for(int j = (NMARKS - 1); j >= 1; j--) { // working backwards here
		int i = j - 1;
		a = tickMarks.at(i);
		b = tickMarks.at(j);

		if(x < a ) sum += (b - a) * ( 0.5 * (Fptr->at(j) + Fptr->at(i)) );
		else {
			// We have reached the bin that contains 'x'
			// Now compute the area under the curve upto the point that includes 'x'
			tmp1 = (b-x) / (b-a);
			tmp2 = (x-a) / (b-a);

			fx = (tmp1 * Fptr->at(i)) + (tmp2 * Fptr->at(j));

			sum += (b - x) * ( 0.5 * (fx + Fptr->at(j)) );
			break; // leave loop since you are done
		}
	}

	return sum;
}




// Function pepares the minorMaps that are needed to perform minorization on
// the FDR data (both global and local)
void FLRClass::setMinorMaps() {
	deque<double> *deqPtr = NULL;
	map<double, vector<double> > localMap;
	map<int, vector<double> > *mapPtr = NULL;
	map<double, vector<double> >::iterator M;
	list<double> scoreList;
	list<double>::iterator L;
	flrStruct curScore;
	int i = 0;

	vector<double> pair; // 0 = score, 1 = global fdr for score

	double s = 0, g = 0; // score, global fdr for score

	// 0 = globalFDR, 1 = localFDR
	for(int iter = 0; iter < 2; iter++) {

		// change the pointers as necessary.
		if(iter == 0) {
			deqPtr = &globalFDR;
			mapPtr = &minorMapGlobal;
		}
		else {
			deqPtr = &localFDR;
			mapPtr = &minorMapLocal;
		}

		scoreList.clear();
		localMap.clear();

		for(i = 0; i < Nreal; i++) {
			curScore = realDeq.at(i);

			s = curScore.deltaScore;

			g = deqPtr->at(i);

			scoreList.push_back(s);

			pair.clear();
			pair.resize(2);
			pair.at(0) = s;
			pair.at(1) = g;
			localMap[ s ] = pair;
		}

		scoreList.sort(); // sort scores from low to high

		i = 0;
		mapPtr->clear();
		for(L = scoreList.begin(); L != scoreList.end(); L++) {
			M = localMap.find( *L ); // find this raw score

			// store the pair(score, globalFDR) into the ret_ptr map.
			// The map is indexed from 0 to N where 0 is the lowest raw score
			// and it's associated globalFDR
			if(M != localMap.end()) mapPtr->insert(std::make_pair(i, M->second));
			i++;
		}
	}
}


// Function performs minorization on the FDR data (both global and local)
void FLRClass::performMinorization() {
	deque<double> *deqPtr = NULL;
	deque<double>::iterator v;
	map<int, vector<double> > *mapPtr = NULL;
	map<int, vector<double> >::iterator M;

	int i,j,N;
	int curStart, curEnd;
	double f_expect;
	double slope;
	bool cont;
	flrStruct *local_flr = NULL;;


	deque<double> x;
	deque<double> f;
	deque<double> fcopy;
	deque<double> forig;
	deque<bool> isMinorPoint;


	// 0 = global FDR, 1 = local FDR
	for(int iter = 0; iter < 2; iter++) {

		// change the pointers as necessary.
		if(iter == 0) {
			deqPtr = &globalFDR;
			mapPtr = &minorMapGlobal;
		}
		else {
			deqPtr = &localFDR;
			mapPtr = &minorMapLocal;
		}

		N = (signed) mapPtr->size();


		// initialize deques for next iteration
		x.clear();
		f.clear();
		fcopy.clear();
		forig.clear();
		isMinorPoint.clear();


		for(i = 0; i < N; i++) {
			M = mapPtr->find(i);

			x.push_back( M->second.at(0) );
			f.push_back( M->second.at(1) );
			forig.push_back( M->second.at(1) );
			fcopy.push_back(0);
			isMinorPoint.push_back(false);
		}


		// find minimum and set all the following points to the tail regression line
		int minId = 0;
		double minVal = f[0];
		for(i = 1; i < N; i++) {
			if( f[i] < minVal ) {
				minVal = f[i];
				minId = i;
			}
		}
		slope = (0.0 - minVal) / ( maxDeltaScore * 1.1 - x[minId] );
		i = minId;
		while(i < N) {
			f[i] = minVal + slope *  (x[i] - x[ minId ]);
			i++;
		}


		// find maximum and set all the preceding points to the minimum
		int maxId = 0;
		double maxVal = f[0];
		i = 1;
		while(x[i] < (x[N-1] / 2.0)) {
		   if(f[i] >= maxVal) {
			  maxVal = f[i];
			  maxId = i;
		   }
		   i++;
		}

	    slope = (maxVal) / ( x[maxId] - x[N-1] );
	    i = maxId - 1;
	    while(i >= 0) {
	      f[i] = maxVal - slope * (x[maxId] - x[i]);
	      i--;
	    }


	    for(i = 0; i < maxId; i++) isMinorPoint[i] = true;
	    curStart = maxId;
	    curEnd = maxId + 1;

	    while( x[curStart] >= x[curEnd] ) curEnd++;


	    while( curStart < (N-1) ) {

	    	slope = ( f[curEnd] - f[curStart] ) / ( x[curEnd] - x[curStart] );

	    	if(slope > 0.0) curEnd++;
	    	else {
				slope = ( f[curEnd] - f[curStart] ) / ( x[curEnd] - x[curStart] );

				cont = true;
				i = curStart+1;
				while( cont && (i < N) && (i < curEnd) ) {

					f_expect = f[curStart] + slope * (x[i] - x[curStart]);

					if( f_expect > f[i] ) cont = false;
					else cont = true;
					i++;
				}

				if(cont) {
					isMinorPoint[curEnd] = true;
					curStart = curEnd;
					curEnd = curStart;
					while( (x[curStart] >= x[curEnd]) && (curEnd < N) ) curEnd++;
				}
				else curEnd++;
	    	}
	    }

	    isMinorPoint[ N-1 ] = true; // check if this value is bigger than the second last minor point
	    for(i = 0; i < N; i++) fcopy[i] = f[i];

	    curStart = 0;
	    curEnd = curStart + 1;
	    while(isMinorPoint[curEnd] == false) curEnd++;

	    while( (curStart < (N-1)) && (curEnd < N) ) {
	    	i = curStart + 1;
	    	slope = (fcopy[curEnd] - fcopy[curStart]) / (x[curEnd] - x[curStart]);

	    	while(i < curEnd) {
	    		f_expect = fcopy[curStart] + slope * (x[i] - x[curStart]);
	    		if(fcopy[i] > f_expect) f[i] = f_expect;
	    		i++;
	    	}

	    	curStart = curEnd;
	    	curEnd = curStart + 1;
	    	while( (isMinorPoint[curEnd] == false) && (curEnd < N) ) curEnd++;
	    }


	    // map back to globalFDR: f[] has minorized FDR values in increasing order of x[], which is deltaScore
	    v = deqPtr->begin();

	    for(i = 0; i < N; i++) {
	    	for(j = 0; j < N; j++) {
	    		local_flr = new flrStruct;
	    		*local_flr = realDeq.at(i);

	    		if(local_flr->deltaScore == x[j]) {
	    			*v = f[j];
	    			break;
	    		}
	    		delete(local_flr); local_flr = NULL;
	    	}
	    	v++;
	    }

	} // end iter for loop
}


// Function to assign FDR values to the specIDs they belong to. Final values
// are assigned to the 'retMap' for the FLRClass
void FLRClass::assignFDRvalues() {
	deque<flrStruct>::iterator real_iter;
	deque<double>::iterator globalFDR_iter, localFDR_iter;
	flrStruct *F = NULL;
	string curID;
	double local_flr = 0;

	real_iter = realDeq.begin();
	globalFDR_iter = globalFDR.begin();
	localFDR_iter = localFDR.begin();

	retMap.clear();

	while(real_iter != realDeq.end()) {

		F = new flrStruct;
		local_flr = *localFDR_iter;
		if(local_flr > 1) local_flr = 1;


		F->globalFLR = *globalFDR_iter;  // global FDR
		F->localFLR = local_flr;        // local FDR
		F->prob = 1.0 - local_flr;       // probability version of luciphor score
		F->specId = real_iter->specId;
		curID = F->specId;
		retMap.insert( pair<string, flrStruct >(curID, *F) );

		delete(F);

		globalFDR_iter++;
		localFDR_iter++;
		real_iter++;
	}
}



// Function returns the FLR results for the given specID value
flrStruct FLRClass::getFLR(string specId) {
	map<string, flrStruct>::iterator m;
	flrStruct ret;

	m = retMap.find(specId);
	if(m != retMap.end()) ret = m->second;
	else {
		// This PSM's best assignment is a decoy so we don't have FLR data for it
		// Fill in the FLR variables with a place holder value
		ret.deltaScore = -1000;
		ret.globalFLR = 1;
		ret.localFLR = 1;
		ret.isDecoy = 1;
		ret.prob = -1;
	}
	return ret;
}




