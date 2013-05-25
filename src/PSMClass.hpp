/*
 * PSMClass.hpp
 *
 *  Created on: Mar 18, 2011
 *      Author: dfermin
 */

#ifndef PSMCLASS_HPP_
#define PSMCLASS_HPP_

#include <string>
#include <map>
#include <vector>
#include <deque>
#include <list>
#include <set>
#include <fstream>
#include "structs.hpp"
using namespace std;

class PSMClass;


class PSMClass {
private:
	string specId;
	string peptide;
	string modPeptide;
	string origModPeptide; // the TPP's predicted sequence for this pepptide
	int charge;
	int scanNum;
	int peakType; // 0 = raw, 1 = scaled (0-100), 2 = median normalized in log scale
	int numPhosphoSites;
	int numSupportingSpectra; // nss = how many spectra match this phospho-peptide
	double mass; // calculated neutral mass of peptide
	double precursor_mass; // neutral precursor mass of peptide
	double numPotentialSites;
	double numPermutations;
	double numDecoyPermutations;
	double iniProb;
	double delta_top2peaks; // holds the difference between the top 2 most intense peak of spectrum
	double min_intensity;
	double max_intensity; // maximum intensity observed from RAW peaks
	double max_intensity_orig; // maximum peak intensity observed in RAW *BEFORE* removing neutral loss peaks
	double median_intensity; // the median peak intensity for the spectrum
	double meanI;
	double varI;
	double nterm_mass;
	double luciphor_deltaScore;
	double scoreTime;
	double maxSiteLevelScore;
	double mz_err; // based upon charge state of peptide, used for peak matching
	double NLprob; // probability score for identified neutral loss peak

	// These are derived from the FLRclass
	double localFLR;
	double globalFLR;
	double luciphorProb;


	bool is_valid_phosphoPSM;
	bool use_for_model;
	bool is_unambiguous;
	bool is_randomized_seq;

	ascoreFinalStruct afs; // holds result for Ascore run

	vector<int> styPos; // array to hold the positions of STY in the peptide sequence

	map<int, double> siteLevelScoreMap; //holds site-level scores, k = site position, v = score

	set<string> phosphoVersionSet;
	set<string> decoySet;

	// these are the different spectrum data formats we use.
	// we use 'spectrum' to point to the right map in the right function
	map<double, vector<double> > raw_spectrum; // spectrum v[0] = raw, v[1] = scaled (0-100), v[2] = log(median normalized)

	list<double> mz;
	list<double> intensities;
	list<double> M_ints_b, M_ints_y, U_ints_local;
	list<double> M_dist_b, M_dist_y, U_dist_local;
	list<double> all_theoPeaks;

	map<double, peakStruct> matchedPeakMap; // holds peaks that were matched
	map<double, peakStruct> unmatchedPeakMap; // holds unmatched peaks
	map<double, peakStruct> ascoreMatchedSpectrum1, ascoreMatchedSpectrum2;
	matchedSpectrumStruct bestSpectrum, nextBestSpectrum;
	scoreStruct bestScore_final, nextBestScore_final;
	modelParamStruct local_params;

public:

	PSMClass(){}; // default ctor
	PSMClass(matchDataStruct *ptr);

	void identify_STY_sites();
	void printVariables(bool tabular);

	void recordSpectrum(SpecStruct spec);
	void normalizeSpectrum();
	void identifyNoisyPeaks();
	void writeSpectrumToDisk();
	void generatePermutations();
	void scorePermutations();
	void classifyPeaks();
	void randomizeSeq();
	void getPrecursorMass();

	void calcUnmatchedPeakParams();
	void calcScore();
	void get_rawScores();
	void reduceNeutralLossPeak();
	void reduceMaxPeak();
	void calcSpectralStats();
	void recordMaxIntensity();
	void runAscore();
	void recordBestSpectrumScores(int J);
	void medianNormalizeIntensities();
	void genDecoys();
	void calcNumDecoyPermutations();

	void pickScores(deque<scoreStruct> &v);
	void processTopHits();
	void calcSiteLevelScores(deque<scoreStruct> &scoreDeq);
	void updateSpectrumMap(map<double, double> newSpec, string whichMap);
	void setSpectrumPtr(string whichMap);
	void write_results(ofstream &outf);
	void write_ascore_results(ofstream &outf);
	void setFLR( flrStruct *ptr );
	void writeAscoreSpectrum();
	void getMinDistance_unmatched(double mz, double intensity);
	void threaded_recordModelingParameters_matched();
	void threaded_recordModelingParameters_UNmatched();
	void threaded_scorePSM();
	void clear();

	bool isKeeper() { return is_valid_phosphoPSM; }
	bool useForModeling() { return use_for_model; }
	bool isUnambiguous() { return is_unambiguous; }
	string getSpecId() { return specId; }
	string getPeptideSeq() { return peptide; }
	double getIniProb() { return iniProb; }
	double getMass() { return mass; }
	double getDeltaScore() { return luciphor_deltaScore; }
	double getMaxIntensity() { return max_intensity; }
	void setParamStruct( modelParamStruct *paramPtr ) { local_params = *paramPtr; }
	void setNSS(int x) { numSupportingSpectra = x; }
	void updateMZerr( double new_mz ) { mz_err = new_mz; }
	int getCharge() { return charge; }
	int getNSS() { return numSupportingSpectra; }
	int getNumMatchedPeaks() { return ((signed) matchedPeakMap.size()); }
	int size() { return ((signed) raw_spectrum.size()); }
	int getScanNumber() { return scanNum; }
	scoreStruct getScoreStruct() { return bestScore_final; }

	map<double, double> getSpectrumMap(string whichMap);

	flrStruct getFLRdata();
	string getSpectrumFileName();
	double getNumPerms();
	list<double> getIntensities(char x, char ionType);
	list<double> getDistances(char x, char ionType);
	list<double> getPeaks(int whichType);
	list<double>* getParamList(char matchType, char ionType, char dataType);
	//list<double> getParamList(char matchType, char ionType, char dataType);

	void debug_printStats();
};


#endif /* PSMCLASS_HPP_ */
