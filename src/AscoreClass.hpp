/*
 * AscoreClass.hpp
 *
 *  Created on: Nov 4, 2011
 *      Author: dfermin
 */

#ifndef ASCORECLASS_HPP_
#define ASCORECLASS_HPP_


#include "globals.hpp"
#include "structs.hpp"
#include <string>
#include <set>
#include <map>
#include <vector>
using namespace std;


class AscoreClass;

// This class is used to perform the Ascore calculations for Luciphor.

class AscoreClass {

private:
	string seq;
	int charge;
	int numPhosphorylations;
	int numPotentialSites;
	double nterm_mass;
	double seq_mass;
	double score;
	double minMZ, maxMZ;
	double minIntensity, maxIntensity;
	double mz_err; // based upon charge state of peptide, used for peak matching



	// Multiple charge states are stored in the vector, element 0 = +1
	// element 4 = +5, etc..
	set<string> b_ion_set, y_ion_set; // holds ion strings in order
	map<string, double> b_ions, y_ions; // holds ion strings and their m/z values

	map<double, double> local_spectrum;
	map<string, double> neutralLossMap; // holds neutral loss masses
	map<double, peakStruct> matched_spectrum; // holds the peaks that were matched (for output)

public:
	AscoreClass(string txt, int Z, double ntm);

	void makeIons();
	void generateIonsMZ(string ion, char ion_type);
	void generate_NL_ionsMZ(string ion, char ion_type);
	void recordMZrange();
	void assignSpectrumMap(map<double, vector<double> > *ptr); // { local_spectrum = src; }
	void getSpectrumAt_X_depth(int numPeaks, map<double, double> *Mptr);
	void recordBestSpectrum(int bestPeakDepth);
	void recordSiteDetermIons(set<string> &ions);
	void removePrecursorNL();
	void generateNL_ionsMZ(string ion, char ion_type);
	map<double, peakStruct> getMatchedSpectrumMap() { return matched_spectrum; }
	map<double, double> getRawSpectrum() { return local_spectrum; }
	double getMass();
	int getNumModSTYs(string txt);
	deque<double> getPeptideScore();
	double getFinalAscore(int optimalPeakDepth);
	double getPhosphoRS_pr();
	double getPeakProb(double i);
	vector<double> getPeptideScoreVec();
};


#endif /* ASCORECLASS_HPP_ */
