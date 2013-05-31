/*
 * FLRClass.hpp
 *
 *  Created on: Jun 22, 2012
 *      Author: dfermin
 */

#ifndef FLRCLASS_HPP_
#define FLRCLASS_HPP_

#include <string>
#include <deque>
#include <iostream>
#include <vector>
#include <map>
#include "globals.hpp"
#include "structs.hpp"

using namespace std;

class FLRClass;

class FLRClass {

private:
	deque<flrStruct> allPSMs;
	deque<flrStruct> realDeq;
	deque<flrStruct> decoyDeq;
	deque<double> tickMarks;
	deque<double> globalFDR;
	deque<double> localFDR;
	deque<double> raw_globalFDR;
	deque<double> raw_localFDR;
	deque<double> f0, f1;
	map<int, vector<double> > minorMapGlobal;
	map<int, vector<double> > minorMapLocal;
	map<string, flrStruct > retMap;


	double maxDeltaScore;
	int Nreal;
	int Ndecoy;
	double Nreal2;
	double Ndecoy2;
	double bw_real;
	double bw_decoy;

public:

	FLRClass(deque<flrStruct> *ptr, double d); // default ctor
	double getVariance(deque<flrStruct> *ptr);
	double getMean(deque<flrStruct> *x);
	double normalDensity(double curTickMark, double curScore, double h);
	double getLocalAUC(double x, int whichF);
	double getGlobalAUC(double x, int whichF);
	flrStruct getFLR(string specId);
	void getBandWidth(string x);
	void evalTickMarks();
	void calcBothFDRs();
	void setMinorMaps();
	void performMinorization();
	void assignFDRvalues();

};



#endif /* FLRCLASS_HPP_ */
