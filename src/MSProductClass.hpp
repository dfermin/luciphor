/*
 * MSProductClass.hpp
 *
 *  Created on: Apr 16, 2011
 *      Author: dfermin
 */

#ifndef MSPRODUCTCLASS_HPP_
#define MSPRODUCTCLASS_HPP_

#include <cmath>
#include <list>
#include <set>
#include <string>
#include <map>
#include <algorithm> // need this for random_shuffle()
#include <fstream>
#include <iostream>
#include "globals.hpp"
#include "structs.hpp"

using namespace std;

class MSProductClass;

class MSProductClass {
private:
	string seq;
	string specId; // don't really need this variable but useful for debugging
	int charge;
	int totNumIons; // holds total number of fragment ions considered for this peptide sequence
	double nterm_mass;
	double seq_mass;
	double min_I;
	double match_threshold; // all peaks with intensities > this will be considered at final scoring phase
	double mz_err; // used for peak matching
	double frac_matched; // # of matched ions / totNumIons
	fstream debug_ionScores; // used to spit out the scores for ions for debugging


	// Multiple charge states are stored in the vector, element 0 = +1
	// element 4 = +5, etc..
	set<string> b_ion_set, y_ion_set; // holds ion strings in order
	map<string, double> b_ions, y_ions;

	map<double, peakStruct> matchedPeaks;
	map<double, peakStruct> unmatchedPeaks;

	map<string, double> ppmErrMap; // k = fragment ion, v = ppm-level err
	map<double, double> local_spectrum;

	SpecStruct s;

public:
	MSProductClass(string the_specId, string txt, int Z, double ntm);

	void makeIons();
	void generateIonsMZ(string ion, char ion_type);
	void generate_NL_ionsMZ(string ion, char ion_type);
	void reduceNeutralLoss(map<double, double> *obs_spectrum_ptr);
	void printIons();
	void recordMatchPeaks(bool forModeling);
	void recordUnmatchedPeaks();
	void keepOnlySiteDetermIons(set<string> &ions);

	void addPeakData(map<double, peakStruct> *targetPtr, char whichMap);
	void getMatchedPeaks(map<double, peakStruct> *Mptr);

	void getUnmatchedPeaks(map<double, double> *srcMapPtr, map<double, peakStruct> *Mptr, map<double, peakStruct> *Uptr);
	double calcSpectrumScore(map<double, peakStruct> *Mpeaks);
	double calcSpectrumScore_HCD(map<double, peakStruct> *Mpeaks);

	void percentilePeaks(double alpha);
	bool assignSpectrumMap(map<double, double> src);
	void updateMZerr(double new_mz_err) { mz_err = new_mz_err; }
	void findMinIntensityThreshold();
	void randomizePeaks();
	void swapSpectra();
	double getMass();
	double getIonMass(string srcStr);
	void assignFragmentIonsMZ(list<double> &ret);
	void calc_ppm_err();


	double get_min_I() { return min_I; }
	double get_match_threshold() { return match_threshold; }
	double getFractionMatched() { return frac_matched; }


	map<double, peakStruct> getPeakMap(char whichMap) {
		if(whichMap == 'm') return matchedPeaks;
		else return unmatchedPeaks;
	}

	map<string, vector<double> > get_NL_peaks();
	scoreStruct scorePermutation();
};


#endif /* MSPRODUCTCLASS_HPP_ */
