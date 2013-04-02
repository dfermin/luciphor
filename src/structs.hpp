/*
 * structs.hpp
 *
 *  Created on: Mar 18, 2011
 *      Author: dfermin
 */

#ifndef STRUCTS_HPP_
#define STRUCTS_HPP_

// All structures go into this file

#include <string>
#include <map>
#include <vector>
#include <deque>
#include <list>
#include "boost/bimap.hpp"

using namespace std;


// Struct to hold all of the data relevant to a peptide-to-spectra-match (PSM)
typedef struct {
	string peptide;
	string modPeptide;
	string specId;
	int scanNum;
	int charge;
	double mass;
	double iniProb; // holds PeptideProphet Probability
	map<int, double> mods; // holds the modificaitons of this peptide
} matchDataStruct;


// Struct mapping file and scan numbers to TPP spectrum_query ids
typedef struct {
	string spectrumFile;
	int scanNumber;
	string specId;
} scanIdStruct;


// Struct holds spectrum taken from pwiz library
typedef struct {
	vector<double> mz;
	vector<double> intensity;
} SpecStruct;


// Struct to hold a fragment of a spectrum. This is used in MSProductClass
typedef struct {
	typedef boost::bimap<double, double> bm_type;
	bm_type bm;

	double mz_start;
	double mz_end;
	list<double> intensity;

} PeakBinStruct;


// Struct to hold the scores obtained for a given spectrum interpretation
// This is used in MSProductClass
typedef struct {
	int topNpeaksConsidered;
	int matchedPeaks;
	int unmatchedPeaks;
	int charge;
	double peakFraction; // matchedPeaks / totalPeakCount
	double score;
	double scoreByPeaks; // score / topNpeaksConsidered
	double fracMatched; // # fraction of theoretical peaks that were matched
	string seq;

} scoreStruct;



// Struct holds model paramters for a charge state.
typedef struct {
	double matched_mean_b;
	double matched_var_b;

	double matched_mean_y;
	double matched_var_y;

	double unMatched_mean;
	double unMatched_var;

	double matched_dist_mean_b;
	double matched_dist_var_b;

	double matched_dist_mean_y;
	double matched_dist_var_y;

	double matched_dist_var_IQR_b;  // for HCD, inter-quartile range
	double matched_dist_var_IQR_y;  // for HCD, inter-quartile range

	double unMatched_dist_mean;
	double unMatched_dist_var;

	// used for non-parametric model in HCD, m/z distance
	// These were required by Hyungwon's nonParamDist.cpp code
	int ntick_dist;
	deque<double> f_dist;
	deque<double> f_dist_U;
	double bw_dist;
	double bw_dist_U;
	deque<double> tickMarks_dist;
	deque<double> tickMarks_dist_U;

	int ntick_int;
	deque<double> f_int_b;
	deque<double> f_int_y;
	deque<double> f_int_U;
	double bw_int_b;
	double bw_int_y;
	double bw_int_U;
	deque<double> tickMarks_int_b;
	deque<double> tickMarks_int_y;
	deque<double> tickMarks_int_U;


} modelParamStruct;



// Struct holds data for consolidating PSMClass objects that all share the same
// phospho-peptide sequence
typedef struct {
	string repSpecId;
	double bestScore;
	vector<string> allSpecIds; // specIds for supporting PSM
} repPSMStruct;



// Struct holds the intensity of a matched peak and it's distance from the theoretical peak
typedef struct {
	double intensity;
	double norm_intensity;
	double MZdistance;
	double mz;
	char ionType; // b or y
	bool hasSTY; // used for determining if the peak should be modeled or not
	string ionStr;
} peakStruct;


// Struct holds the spectrum for a PSMClass object with the
// matched peaks highlighted
typedef struct {
	map<double, peakStruct> matched_ions; // holds all info about matched peaks
	map<string, double> IscoreMap; // holds intensity score for matched peaks
	map<string, double> DscoreMap; // holds mz distance score for matched peak
	map<string, double> FinalScoreMap; // holds the final score for this matched peak
} matchedSpectrumStruct;



// Struct to hold information about the PSMs that will be used for modeling
// This is usedin PepXMLClass::prunePSMdeq
typedef struct {
	int totalZ; // total number of PSMs for the given charge state
	int modelZ; // total number of PSMs for the given charge state that qualify for modeling
	double maxScore; // the best TPP score observed from all of the PSMs of this charge state
}prunePSMstruct;



// Struct is used to hold intermediate data for the Ascore calculation
typedef struct {
	int numPeaksPerBin;
	double negLogProb;
	int numPeaksTotal;
	int numMatchedPeaks;
	double peptideScore;
} ascoreStruct;



// Struct used to store the final ascore results for a given phospho permutation
typedef struct {
	int numPeaksPerBin;
	int numMatchedPeaks1;
	int numMatchedPeaks2;
	double peptideScore;
	double maxScoreDiff;
	string seqBest;
	string nextSeq;
} ascoreFinalStruct;



// Struct to hold the necessary variables to use Hyungwon's FDR code
typedef struct {
	string specId;
	bool isDecoy;
	double deltaScore;
	double localFLR;
	double globalFLR;
	double prob;
} flrStruct;






#endif /* STRUCTS_HPP_ */
