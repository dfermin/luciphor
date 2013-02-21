/*
 * AscoreClass.cpp
 *
 *  Created on: Nov 4, 2011
 *      Author: dfermin
 */


#include <map>
#include <string>
#include <cstdlib>
#include <cmath>
#include <list>
#include <iostream>
#include <boost/regex.hpp>
#include "PSMClass.hpp"
#include "AscoreClass.hpp"
#include "globals.hpp"
#include "boost/bimap.hpp"
#include "statsFunctions.hpp"
using namespace std;
using namespace boost;

const double H = 1.00727;
const double OH = 17.002745;
const double e = 0.00054858026;
const double H2O = 18.010015;

const double MZ_WIN = 100.00;
const double MZ_TOLERANCE = 0.6;   // this is constant as defined in the Ascore paper

// default Constructor
AscoreClass::AscoreClass(string txt, int Z, double ntm) {

	seq = txt;
	charge = Z;
	nterm_mass = ntm;
	maxMZ = 0.0;
	minMZ = 0.0;
	maxIntensity = 0;
	minIntensity = 0;

	mz_err = MZ_TOLERANCE;

	neutralLossMap[ "-H2O" ] =   -1.0 * H2O;
	neutralLossMap[ "+H2O" ] = H2O;
	neutralLossMap[ "-H3PO4" ] = -97.976895;
	neutralLossMap[ "-HPO3" ] =  -79.966330;
	neutralLossMap[ "-NH3" ] =   -17.026549;


	// figure out the number of reported and potential phosphorylation sites
	// in the sequence
	numPotentialSites = 0;
	numPhosphorylations = 0;
	for(int i = 0; i < (signed)seq.length(); i++) {
		char c = seq.at(i);
		char C = toupper(c);
		if( (C == 'S') || (C == 'T') || (C == 'Y') ) numPotentialSites++;
		if( (c == 's') || (c == 't') || (c == 'y') ) numPhosphorylations++;
		if( !isalpha(c) ) numPhosphorylations++; // decoy modification
	}


	seq_mass = getMass() + nterm_mass;
	makeIons(); // fragment the given sequence into B and Y ions
}

// Function returns the mass of the current peptide string recorded in 'seq'
double AscoreClass::getMass() {
	double ret = H + H2O;
	int N = (signed)seq.length();

	for(int i = 0; i < N; i++) {
		ret += AAmass[ seq.at(i) ];
	}

	return(ret);
}



// Function fragments the sequence and stores its ions
void AscoreClass::makeIons() {
	string b, y;

	int N = (signed)seq.length();

	for(int i = 0; i < N; i++) {
		b = seq.substr(0, (i+1) );
		y = seq.substr( (i+1) );

		if( (signed)b.length() == N ) continue;
		if( (signed)y.length() == N ) continue;

		generateIonsMZ(b, 'b');

		generateIonsMZ(y, 'y');
	}
}


// Function analyzes the given ion string to determine if it can undergo a
// neutral loss of any kind. If it can, that potential neutral loss is recorded.
//void AscoreClass::generate_NL_ionsMZ(string ion, char ion_type) {
//	double mass, mz_value;
//	int N = (signed) ion.length();
//	string Nstr = int2string(N);
//	string new_ion;
//	bool hasPhospho = false;
//
//	int S, T, Y;
//
//	// compute mass of ion
//	mass = 0.0;
//	for(int i = 0; i < N; i++) { mass += AAmass[ ion.at(i) ]; }
//
//	// determine if the ion contains an STY letter, if it does, then the loss
//	// of a phospho-group is possible. Otherwise, the ion can only lose
//	// NH3 or H2O/H
//	S = 0; T = 0; Y = 0; // count the number of STY in the ion
//	for(int i = 0; i < N; i++) {
//		if( ion.at(i) == 's') S++;
//		if( ion.at(i) == 't') T++;
//		if( ion.at(i) == 'y') Y++;
//	}
//	if( (S+T+Y) > 0 ) hasPhospho = true; // has at least 1 phosphorylated AA
//
//
//	double extraProton = 0; // Computes: ( H * (z-1) ) <- you need this for dealing with multiple charge states;
//	for(int z = 1; z < charge + 1; z++) {
//
//		extraProton = ( H * (z-1) );
//
//		if(ion_type == 'b') {
//
//			// this ion contains a phospho group that can undergo neutral loss
//			if(hasPhospho) {
//				new_ion.clear();
//				new_ion = "b^" + Nstr + ":" + ion + "-H3PO4";
//				if(z > 1) new_ion += "/+" + int2string(z);
//
//				mz_value = mass + H - e + extraProton;
//				mz_value += neutralLossMap[ "-H3PO4" ];
//				mz_value /= z;
//
//				// record the -H3PO4 neutral loss
//				if( (mz_value > MIN_MZ) && (N > 1) ) {
//					b_ions[ new_ion ] = mz_value;
//					b_ion_set.insert(new_ion);
//				}
//			} // end if(hasPhospho)
//		} // end b-ion
//		else if(ion_type == 'y') {
//
//			// this ion contains a phospho group that can undergo neutral loss
//			if(hasPhospho) {
//				new_ion.clear();
//				new_ion = "y^" + Nstr + ":" + ion + "-H3PO4";
//				if(z > 1) new_ion += "/+" + int2string(z);
//
//				mz_value = mass + H2O + H + extraProton;
//				mz_value += neutralLossMap[ "-H3PO4" ];
//				mz_value /= z;
//
//				// record the -H3PO4 neutral loss
//				if( (mz_value > MIN_MZ) && (N > 1) ) {
//					y_ions[ new_ion ] = mz_value;
//					y_ion_set.insert(new_ion);
//				}
//			} // end if(hasPhospho)
//		} // end y-ion
//	}
//}


// function generates the theoretical m/z value for the given ion string.
// All ion variations (neutral losses, etc..) are stored in the b_ions and y_ions
// maps of the MSProductClass object
void AscoreClass::generateIonsMZ(string ion, char ion_type) {
	double mass;
	double mz_value; // the m/z value of the current ion string
	int N = (signed)ion.length();
	string Nstr = int2string(N);
	string new_ion;

	// The code '( H * (z-1) )' is used to add a proton to the ion
	// for multiple charge states. For every additional increase in an ion's
	// charge state, you need to add 1 proton. That additional proton is where
	// the charge state comes from.


	for(int z = 1; z < (charge + 0); z++) { // was (charge + 1) 2011.12.05
		mass = 0.0;
		mz_value = 0.0;

		for(int i = 0; i < N; i++) {
			mass += AAmass[ ion.at(i) ];
		}

		if(ion_type == 'b') {
			mz_value = mass + H - e + ( H * (z-1) ) + nterm_mass;
			mz_value /= z;
			if( (mz_value > MIN_MZ) && (N > 1)) {

				new_ion.clear();
				new_ion = "b^" + Nstr + ":" + ion;
				if(z > 1) { new_ion += "/+" + int2string(z); }

				b_ions[ new_ion ] = mz_value;
				b_ion_set.insert(new_ion);
			}
		}
		else if(ion_type == 'y') {
			mz_value = mass + H2O + H + ( H * (z-1) );
			mz_value /= z;

			if( (mz_value > MIN_MZ) && (N > 1)) {

				new_ion.clear();
				new_ion = "y^" + Nstr + ":" + ion;
				if(z > 1) { new_ion += "/+" + int2string(z); }

				y_ions[ new_ion ] = mz_value;
				y_ion_set.insert(new_ion);
			}
		}
	}
}





// Function records the maximum M/Z value in the stored spectrum
void AscoreClass::recordMZrange() {
	map<double, double>::iterator curPeak;
	list<double> *Lptr1, *Lptr2;
	int a1, remainder;
	double x1;

	// Get the max m/z for spectrum, also record spectrum in bi-directional map
	Lptr1 = new list<double>;
	Lptr2 = new list<double>;

	for(curPeak = local_spectrum.begin(); curPeak != local_spectrum.end(); curPeak++) {
		Lptr1->push_back( curPeak->first );
		Lptr2->push_back( curPeak->second );
	}
	Lptr1->sort();
	Lptr2->sort();

	// get minMZ
	a1 = (int) Lptr1->front();
	remainder = a1 % 100;
	x1 = a1 - remainder;
	minMZ = fmax( MIN_MZ, x1 );

	// get maxMZ
	a1 = (int) Lptr1->back();
	remainder = a1 % 100;
	x1 = a1 - remainder + 100;
	maxMZ = x1;


	minIntensity = Lptr2->front();
	maxIntensity = Lptr2->back();

	delete(Lptr1);
	delete(Lptr2);
}



// Function creates an artificial spectrum from local_spectrum consisting of only the
// top 'numPeaks' found within each MZ_WIN of local_spectrum.
// This returned map is what will be used for calculating Ascore
void AscoreClass::getSpectrumAt_X_depth(int numPeaks, map<double, double> *Mptr) {
	double curMZ, nextMZ;
	double mz, intensity;
	map<double, double>::iterator curPeak;
	map<double, double> bestPeaksMap;
	map<double, list<double> > deepMap;
	map<double, list<double> >::iterator dmIter;
	list<double> *Lptr = NULL;
	list<double>::iterator iterL;

	double start_mz = round_dbl(minMZ, 0);
	double end_mz   = round_dbl(maxMZ, 0);

	for(curMZ = start_mz; curMZ < (end_mz + MZ_WIN); curMZ += MZ_WIN) {
		nextMZ = curMZ + MZ_WIN;

		mz = 0;
		intensity = 0;
		deepMap.clear();
		bestPeaksMap.clear();

		for(curPeak = local_spectrum.begin(); curPeak != local_spectrum.end(); curPeak++) {
			// according to Ascore paper "only one peak per one m/z unit was allowed"
			mz = round_dbl(curPeak->first, 0);
			intensity = curPeak->second;

			if( (mz >= curMZ) && (mz < nextMZ) ) { // candiate peak
				dmIter = deepMap.find(mz);
				if(dmIter == deepMap.end()) deepMap[ mz ].push_back(intensity);
				else  dmIter->second.push_back(intensity);
			}
		}

		Lptr = new list<double>;

		//select the most intense peak in each list of deepMap
		for(dmIter = deepMap.begin(); dmIter != deepMap.end(); dmIter++) {
			mz = dmIter->first;
			dmIter->second.sort();
			bestPeaksMap[ mz ] = dmIter->second.back();
			Lptr->push_back( bestPeaksMap[ mz ] );
		}

		// bestPeakMap now contains the most intense peaks for this mz window range
		// pick the top 'i' most intense peaks in this window
		Lptr->sort();
		Lptr->reverse();
		int i = numPeaks;
		iterL = Lptr->begin();
		while(i > 0) {
			intensity = *iterL;
			for(curPeak = bestPeaksMap.begin(); curPeak != bestPeaksMap.end(); curPeak++) {
				if(curPeak->second == intensity) {
					mz = curPeak->first;
					Mptr->insert(pair<double, double>(mz, intensity));
					break;
				}
			}
			iterL++;
			i--;
		}

		delete(Lptr); Lptr = NULL;
	}
}





// Function records the spectrum corresponding to the given peak depth
void AscoreClass::recordBestSpectrum(int bestPeakDepth) {
	vector<double> ret;
	map<double, double>::iterator curPeak;
	map<string, double>::iterator curTheoPeak;
	map<double, double> *curPeakMap = NULL;
	ascoreStruct *ASS = NULL;
	peakStruct *psPtr = NULL;

	double expected_mz, a, b, mz, intensity;

	// holds all the peaks that match a particular m/z value. The
	map<double, vector<peakStruct> > matched_ions;

	// create the map for this iteration
	curPeakMap = new map<double, double>;
	getSpectrumAt_X_depth(bestPeakDepth, curPeakMap);

	// prepare struct for storing results
	ASS = new ascoreStruct();
	ASS->negLogProb = 0.0;
	ASS->numMatchedPeaks = 0.0;
	ASS->numPeaksPerBin = (signed) curPeakMap->size();

	matched_ions.clear(); // prep for next iteration

	// iterate over b-ions
	for(curTheoPeak = b_ions.begin(); curTheoPeak != b_ions.end(); curTheoPeak++) {
		expected_mz = curTheoPeak->second;
		a = expected_mz - mz_err;
		b = expected_mz + mz_err;

		for(curPeak = curPeakMap->begin(); curPeak != curPeakMap->end(); curPeak++) {
			mz = curPeak->first;
			intensity = curPeak->second;

			if( (mz >= a) && (mz <= b) ) {
				psPtr = new peakStruct;
				psPtr->ionStr = curTheoPeak->first;
				psPtr->ionType = 'b';
				psPtr->MZdistance = mz; // we store the m/z value for this peak in MZdistance in this function
				psPtr->intensity = intensity;
				matched_ions[ expected_mz ].push_back( *psPtr );
				delete(psPtr);
			}
		}
	}


	// iterate over y-ions
	for(curTheoPeak = y_ions.begin(); curTheoPeak != y_ions.end(); curTheoPeak++) {
		expected_mz = curTheoPeak->second;
		a = expected_mz - mz_err;
		b = expected_mz + mz_err;

		for(curPeak = curPeakMap->begin(); curPeak != curPeakMap->end(); curPeak++) {
			mz = curPeak->first;
			intensity = curPeak->second;

			if( (mz >= a) && (mz <= b) ) {
				psPtr = new peakStruct;
				psPtr->ionStr = curTheoPeak->first;
				psPtr->ionType = 'y';
				psPtr->MZdistance = mz; // we store the m/z value for this peak in MZdistance in this function
				psPtr->intensity = intensity;
				matched_ions[ expected_mz ].push_back( *psPtr );
				delete(psPtr);
			}
		}
	}


	/*
	 * The matched_ions map _MAY_ contain multiple peaks assigned to the same
	 * theoretical peak. Here we filter those out and choose the most intense
	 * peak as the representative one for this PSM
	 */
	map<double, vector<peakStruct> >::iterator curPeakVec;
	vector<peakStruct>::iterator v_iter;
	list<double> *Lptr = NULL;
	int N = 0;
	double maxI = 0.0;
	for(curPeakVec = matched_ions.begin(); curPeakVec != matched_ions.end(); curPeakVec++) {

		N = (signed) curPeakVec->second.size();
		if(N == 1) {
			psPtr = new peakStruct;
			*psPtr = curPeakVec->second.at(0);
		}
		else { // more than one peak in the running
			Lptr = new list<double>;
			for(v_iter = curPeakVec->second.begin(); v_iter != curPeakVec->second.end(); v_iter++)
				Lptr->push_back( v_iter->intensity );

			Lptr->sort();
			maxI = Lptr->back();
			delete(Lptr);

			for(v_iter = curPeakVec->second.begin(); v_iter != curPeakVec->second.end(); v_iter++) {
				if(v_iter->intensity == maxI) {
					psPtr = new peakStruct;
					*psPtr = *v_iter;
				}
			}
		}

		matched_spectrum[ psPtr->MZdistance ] = *psPtr;
		delete(psPtr);
	}
}



// Function returns the peptide score for the currently assigned phospho permutation
double AscoreClass::getPeptideScore() {

	map<double, double>::iterator curPeak;
	map<string, double>::iterator curTheoPeak;
	map<double, double> *curPeakMap = NULL;

	double expected_mz, a, b, mz;
	set<double> matched_ions;
	double ret = 0.0;
	double negLogProb = 0.0;
	double prob = 0.0;
	int numMatchedPeaks = 0;

	double depthWt[] = { 0, 0.5, 0.75, 1, 1, 1, 1, 0.75, 0.5, 0.25, 0.25 };

	// skip this PSM, it is an unambiguous case
	if(numPotentialSites == numPhosphorylations) return 1000;
	
	
	for(int numPeaks = 1; numPeaks <= 10; numPeaks++) { // number of peaks per m/z window

		// create the map for this iteration
		curPeakMap = new map<double, double>;
		getSpectrumAt_X_depth(numPeaks, curPeakMap);

		matched_ions.clear(); // prep for next iteration

		// iterate over b-ions
		for(curTheoPeak = b_ions.begin(); curTheoPeak != b_ions.end(); curTheoPeak++) {
			expected_mz = curTheoPeak->second;
			a = expected_mz - mz_err;
			b = expected_mz + mz_err;

			for(curPeak = curPeakMap->begin(); curPeak != curPeakMap->end(); curPeak++) {
				mz = curPeak->first;
				if( (mz >= a) && (mz <= b) ) matched_ions.insert(expected_mz);
			}
		}


		// iterate over y-ions
		for(curTheoPeak = y_ions.begin(); curTheoPeak != y_ions.end(); curTheoPeak++) {
			expected_mz = curTheoPeak->second;
			a = expected_mz - mz_err;
			b = expected_mz + mz_err;

			for(curPeak = curPeakMap->begin(); curPeak != curPeakMap->end(); curPeak++) {
				mz = curPeak->first;
				if( (mz >= a) && (mz <= b) ) matched_ions.insert(expected_mz);
			}
		}

		numMatchedPeaks = (signed)matched_ions.size();

		// compute the cumulative binomial probability at the current peak depth
		int k = numMatchedPeaks;
		if(k == 0) { negLogProb = 0; }
		else {
			int N = ((signed) b_ions.size()) + ((signed) y_ions.size());

			double pr = (((double) numPeaks) / MZ_WIN) * g_MZ_ERR; // this equation is from Alexey's old python code
			//double pr = getPhosphoRS_pr();

			prob = cum_binomial_prob(N, k, pr);
			negLogProb = ( -10.0 * log(prob) ) * depthWt[ numPeaks ];
		}

		ret += negLogProb;
		delete(curPeakMap); curPeakMap = NULL;
	}

	return(ret);
}


// Function returns a vector of the peptide score (without weights) for each peak depth
vector<double> AscoreClass::getPeptideScoreVec() {
	vector<double> ret;
	map<double, double>::iterator curPeak;
	map<string, double>::iterator curTheoPeak;
	map<double, double> *curPeakMap = NULL;
	ascoreStruct *ASS = NULL;

	double expected_mz, a, b, mz;
	set<double> matched_ions;

	for(int numPeaks = 1; numPeaks <= 10; numPeaks++) { // number of peaks per m/z window

		// create the map for this iteration
		curPeakMap = new map<double, double>;
		getSpectrumAt_X_depth(numPeaks, curPeakMap);

		// prepare struct for storing results
		ASS = new ascoreStruct();
		ASS->negLogProb = 0.0;
		ASS->numMatchedPeaks = 0.0;
		ASS->numPeaksPerBin = numPeaks;

		matched_ions.clear(); // prep for next iteration

		// iterate over b-ions
		for(curTheoPeak = b_ions.begin(); curTheoPeak != b_ions.end(); curTheoPeak++) {
			expected_mz = curTheoPeak->second;
			a = expected_mz - mz_err;
			b = expected_mz + mz_err;

			for(curPeak = curPeakMap->begin(); curPeak != curPeakMap->end(); curPeak++) {
				mz = curPeak->first;
				if( (mz >= a) && (mz <= b) ) matched_ions.insert(expected_mz);
			}
		}


		// iterate over y-ions
		for(curTheoPeak = y_ions.begin(); curTheoPeak != y_ions.end(); curTheoPeak++) {
			expected_mz = curTheoPeak->second;
			a = expected_mz - mz_err;
			b = expected_mz + mz_err;

			for(curPeak = curPeakMap->begin(); curPeak != curPeakMap->end(); curPeak++) {
				mz = curPeak->first;
				if( (mz >= a) && (mz <= b) ) matched_ions.insert(expected_mz);
			}
		}

		ASS->numMatchedPeaks = (signed)matched_ions.size();

		// compute the cumulative binomial probability
		int k = ASS->numMatchedPeaks;
		if(k == 0) { ASS->negLogProb = 0; }
		else {
			int N = ((signed) b_ions.size()) + ((signed) y_ions.size());
			double pr = (((double) numPeaks) / MZ_WIN) * g_MZ_ERR; // this equation is from Alexey's old python code
			//double pr = getPhosphoRS_pr();
			//double pr = ((double)numPeaks) / MZ_WIN;  // used in original Ascore manuscript

			double prob = cum_binomial_prob(N, k, pr);
			ASS->negLogProb = -10.0 * log(prob);
		}
		ret.push_back( ASS->negLogProb );
		delete(ASS); ASS = NULL;
		delete(curPeakMap); curPeakMap = NULL;
	} // done iterating over peak depths

	return ret;
}




// Function computes the Ascore for the permutation. This is the final ascore value
double AscoreClass::getFinalAscore(int optimalPeakDepth) {
	double ret;
	map<double, double>::iterator curPeak;
	map<string, double>::iterator curTheoPeak;
	map<string, double> *ionMapPtr = NULL;
	map<double, double> *curPeakMap = NULL;
	map<double, peakStruct>::iterator iter_m;
	peakStruct *matchedPeak = NULL;
	int matched_ions = 0;
	double mz, a, b, intensity, expected_mz;
	string curIon;

	// a special condition for our Ascore to handle
	if(numPotentialSites == numPhosphorylations) return 1000;


	curPeakMap = new map<double, double>;
	getSpectrumAt_X_depth(optimalPeakDepth, curPeakMap);

	getSiteDetermIons();

	/********************* Begin peak matching part ***************************
	 *
	 * Because Ascore rounds peak m/z values, some theoretical peaks will match
	 * multiple observed peaks because of their close proximity.
	 * Therefore we have to pick the best observed peak to be assigned to a
	 * theoretical peak.
	 */
	map<double, deque<peakStruct> > candMatchedPeaks; // k = theo_peak v = all cand peaks
	map<double, deque<peakStruct> >::iterator cmpIter;
	deque<peakStruct> *tmpDeq = NULL;
	deque<peakStruct>::iterator iterD;

	// b-ions
	candMatchedPeaks.clear();
	for(curTheoPeak = b_ions.begin(); curTheoPeak != b_ions.end(); curTheoPeak++) {
		expected_mz = curTheoPeak->second;
		curIon = curTheoPeak->first;
		a = expected_mz - mz_err;
		b = expected_mz + mz_err;

		for(curPeak = curPeakMap->begin(); curPeak != curPeakMap->end(); curPeak++) {
			mz = curPeak->first;
			intensity = curPeak->second;

			if( (mz >= a) && (mz <= b) ) {
				matchedPeak = new peakStruct;
				matchedPeak->intensity = intensity;
				matchedPeak->norm_intensity = (intensity / maxIntensity) * 100.00;
				matchedPeak->ionStr = curIon;
				matchedPeak->MZdistance = mz - expected_mz;
				matchedPeak->hasSTY = containsSTY(curIon);
				matchedPeak->mz = mz;

				matchedPeak->ionType = 'b';

				// see if this entry already exists in candMatchedPeaks
				cmpIter = candMatchedPeaks.find(expected_mz);
				if(cmpIter == candMatchedPeaks.end()) {
					tmpDeq = new deque<peakStruct>;
					tmpDeq->push_back( *matchedPeak );
					candMatchedPeaks[ expected_mz ] = *tmpDeq;
					delete(tmpDeq);
				}
				else candMatchedPeaks[ expected_mz ].push_back(*matchedPeak);

				delete(matchedPeak);
			}
		}
	}

	// pick best observed peak for each theoretical peak in candMatchedPeaks
	for(cmpIter = candMatchedPeaks.begin(); cmpIter != candMatchedPeaks.end(); cmpIter++) {
		int N = (signed) cmpIter->second.size();

		if(N == 1) { // only one peak for this theoretical m/z value
			matchedPeak  = new peakStruct;
			*matchedPeak = cmpIter->second.at(0);

			// see if we have already recorded this as a matched peak
			iter_m = matched_spectrum.find(matchedPeak->mz);
			if(iter_m == matched_spectrum.end())
				matched_spectrum[ matchedPeak->mz ] = *matchedPeak;
			else {
				if(matchedPeak->intensity > iter_m->second.intensity)
					matched_spectrum[ matchedPeak->mz ] = *matchedPeak;
			}
			delete(matchedPeak);
		}
		else { // at least 2 of the observed peaks are assigned to this theoretical peak
			   // pick the one with the highest intensity value
			list<double> *Lptr = new list<double>; // holds intensities
			map<double, peakStruct> *Mptr = new map<double, peakStruct>;
			for(iterD = cmpIter->second.begin(); iterD != cmpIter->second.end(); iterD++) {
				Lptr->push_back( iterD->intensity );
				Mptr->insert(pair<double, peakStruct>( iterD->intensity, *iterD) );
			}

			Lptr->sort(); // sorted low to high
			intensity = Lptr->back(); // get highest intensity

			matchedPeak = new peakStruct;
			iter_m = Mptr->find(intensity);
			*matchedPeak = iter_m->second; // get the peak associated with this intensity

			// see if we have already recorded this as a matched peak
			iter_m = matched_spectrum.find(matchedPeak->mz);
			if(iter_m == matched_spectrum.end())
				matched_spectrum[ matchedPeak->mz ] = *matchedPeak;
			else {
				if(matchedPeak->intensity > iter_m->second.intensity)
					matched_spectrum[ matchedPeak->mz ] = *matchedPeak;
			}
			delete(matchedPeak);
			delete(Lptr);
			delete(Mptr);
		}
	}


	// y-ions
	candMatchedPeaks.clear();
	for(curTheoPeak = y_ions.begin(); curTheoPeak != y_ions.end(); curTheoPeak++) {
		expected_mz = curTheoPeak->second;
		curIon = curTheoPeak->first;
		a = expected_mz - mz_err;
		b = expected_mz + mz_err;

		for(curPeak = curPeakMap->begin(); curPeak != curPeakMap->end(); curPeak++) {
			mz = curPeak->first;
			intensity = curPeak->second;

			if( (mz >= a) && (mz <= b) ) {
				matchedPeak = new peakStruct;
				matchedPeak->intensity = intensity;
				matchedPeak->norm_intensity = (intensity / maxIntensity) * 100.00;
				matchedPeak->ionStr = curIon;
				matchedPeak->MZdistance = mz - expected_mz;
				matchedPeak->hasSTY = containsSTY(curIon);
				matchedPeak->mz = mz;

				matchedPeak->ionType = 'y';

				// see if this entry already exists in candMatchedPeaks
				cmpIter = candMatchedPeaks.find(expected_mz);
				if(cmpIter == candMatchedPeaks.end()) {
					tmpDeq = new deque<peakStruct>;
					tmpDeq->push_back( *matchedPeak );
					candMatchedPeaks[ expected_mz ] = *tmpDeq;
					delete(tmpDeq);
				}
				else candMatchedPeaks[ expected_mz ].push_back(*matchedPeak);

				delete(matchedPeak);
			}
		}
	}

	// pick best observed peak for each theoretical peak in candMatchedPeaks
	for(cmpIter = candMatchedPeaks.begin(); cmpIter != candMatchedPeaks.end(); cmpIter++) {
		int N = (signed) cmpIter->second.size();

		if(N == 1) { // only one peak for this theoretical m/z value
			matchedPeak  = new peakStruct;
			*matchedPeak = cmpIter->second.at(0);

			// see if we have already recorded this as a matched peak
			iter_m = matched_spectrum.find(matchedPeak->mz);
			if(iter_m == matched_spectrum.end())
				matched_spectrum[ matchedPeak->mz ] = *matchedPeak;
			else {
				if(matchedPeak->intensity > iter_m->second.intensity)
					matched_spectrum[ matchedPeak->mz ] = *matchedPeak;
			}
			delete(matchedPeak);
		}
		else { // at least 2 of the observed peaks are assigned to this theoretical peak
			   // pick the one with the highest intensity value
			list<double> *Lptr = new list<double>; // holds intensities
			map<double, peakStruct> *Mptr = new map<double, peakStruct>;
			for(iterD = cmpIter->second.begin(); iterD != cmpIter->second.end(); iterD++) {
				Lptr->push_back( iterD->intensity );
				Mptr->insert(pair<double, peakStruct>( iterD->intensity, *iterD) );
			}

			Lptr->sort(); // sorted low to high
			intensity = Lptr->back(); // get highest intensity

			matchedPeak = new peakStruct;
			iter_m = Mptr->find(intensity);
			*matchedPeak = iter_m->second; // get the peak associated with this intensity

			// see if we have already recorded this as a matched peak
			iter_m = matched_spectrum.find(matchedPeak->mz);
			if(iter_m == matched_spectrum.end())
				matched_spectrum[ matchedPeak->mz ] = *matchedPeak;
				//matched_spectrum.insert(pair<double, peakStruct>(matchedPeak->mz, *matchedPeak));
			else {
				if(matchedPeak->intensity > iter_m->second.intensity)
					matched_spectrum[ matchedPeak->mz ] = *matchedPeak;
			}
			delete(matchedPeak);
			delete(Lptr);
			delete(Mptr);
		}
	}
	/*********************** End peak matching part ***************************/


	if(g_DEBUG_MODE) {
		cerr << endl << seq << endl;
		cerr << "m/z\tion string\tintensity\n";
		for(map<double, peakStruct>::iterator I = matched_spectrum.begin(); I != matched_spectrum.end(); I++) {
			cerr << I->first << "\t" << I->second.ionStr << "\t" << I->second.intensity << endl;
		}
	}

	matched_ions = (signed) matched_spectrum.size();

	if(matched_ions == 0) ret = 0.0;
	else {
		double k = matched_ions;
		int N = (signed)b_ion_set.size() + (signed)y_ion_set.size();

		// this is what is published in the Ascore paper
		double pr = ( ((double)optimalPeakDepth) / MZ_WIN ) * g_MZ_ERR; // Oct 24, 2012 our pr value is based upon Alexey's old python code
		//double pr = ( ((double)optimalPeakDepth) / MZ_WIN ); // used in original Ascore paper
		//double pr = getPhosphoRS_pr();




		double prob = cum_binomial_prob(N, k, pr);
		ret = -10.0 * log(prob);
	}

	delete(curPeakMap);
	return ret;
}





// Function figures out the site-determining ions for the given phospho-peptide sequence
int AscoreClass::getSiteDetermIons() {
	vector<string> B;
	vector<string> Y;
	vector<int> M_b, M_y;
	list<int> phosphoPos;
	int N = 0;
	int ret = 0;
	int ctr = 0;
	string b_ion_str, y_ion_str, tmp, rev_str;
	size_t found;
	string modChars = "sty234567890@#$%&;?~";

	boost::regex ion_regex("^([by].\\d+:\\w+)(\/\+\\d+)?");
	boost::smatch matches;

	string srchStr = seq; // assign the modPeptide sequence of this AscoreClass instance

	N = (signed)srchStr.length(); // length of sequence

	// initialize vectors to all zeros
	for(int i = 0; i < N; i++) {
		M_b.push_back(0);
		M_y.push_back(0);
	}

	// Label M_b with the cumulative sum of phospho sites in the peptide along
	// the b-ion series
	for(int i = 0; i < N; i++) {
		char c = srchStr.at(i);
		found = modChars.find(c);
		if(found != string::npos) ctr++;
		M_b.at(i) = ctr;
	}



	// Label M_y with the cumulative sum of phospho sites in the peptide along
	// the y-ion series. We reverse the sequence for the y-ions. It's easier
	// than figuring out the math in reverse.
	ctr = 0;
	for(int i = N-1; i > -1; i--) {
		char c = srchStr.at(i);
		found = modChars.find(c);
		if(found != string::npos) ctr++;
		M_y.at(i) = ctr;
	}


	/*
	 * Wherever the M_b and M_y vectors have a number equal to 'numPhosphorylations'
	 * for this peptide, that is the ion ladder number for a site-determining ion
	 */
	// B-ions
	for(int i = 1; i < N-1; i++) { // we do not include/keep the last ion (which is the whole peptide)
		if( M_b.at(i) == numPhosphorylations ) {
			b_ion_str = srchStr.substr(0, (i+1) );
			tmp = "b^" + int2string((i+1)) + ":" + b_ion_str;
			B.push_back(tmp);
		}
	}


	// Y-ions (we have to work from the C-term (right-hand-side) of the string)
	for(int i = N-1; i > 0; i--) {
		if(M_y.at(i) == numPhosphorylations) {
			int j = N - i;
			y_ion_str = srchStr.substr(i, j);
			tmp = "y^" + int2string(j) + ":" + y_ion_str;
			Y.push_back(tmp);
		}
	}


	/*
	 * Now that we have the site specific ion strings, we will clean out
	 * the 'b/y_ion_set' sets and store only these cases
	 */
	map<string, double> *tmpIonMap = NULL;
	set<string> *tmpStrSet = NULL;
	set<string>::iterator setIter;
	vector<string>::iterator v;

	// b-ions
	tmpIonMap = new map<string, double>;
	tmpStrSet = new set<string>;
	for(v = B.begin(); v != B.end(); v++) {
		// find this site-determining ion inside of the original b_ion_set
		// Then record it's various forms (ie: different charge states)
		for(setIter = b_ion_set.begin(); setIter != b_ion_set.end(); setIter++) {

			boost::regex_match(*setIter, matches, ion_regex);
			tmp.clear();
			tmp.assign(matches[1].first, matches[1].second);

			if(tmp.compare(*v) == 0) {
				tmpIonMap->insert(pair<string, double>(*setIter, b_ions[ *setIter ]));
				tmpStrSet->insert(*setIter);
			}
		}

	}
	b_ion_set.clear();
	b_ion_set = *tmpStrSet;
	delete(tmpStrSet);
	b_ions.clear();
	b_ions = *tmpIonMap;
	delete(tmpIonMap);

	// y-ions
	tmpIonMap = new map<string, double>;
	tmpStrSet = new set<string>;
	for(v = Y.begin(); v != Y.end(); v++) {
		// find this site-determining ion inside of the original y_ion_set
		// Then record it's various forms (ie: different charge states)
		for(setIter = y_ion_set.begin(); setIter != y_ion_set.end(); setIter++) {

			boost::regex_match(*setIter, matches, ion_regex);
			tmp.clear();
			tmp.assign(matches[1].first, matches[1].second);

			if(tmp.compare(*v) == 0) {
				tmpIonMap->insert(pair<string, double>(*setIter, y_ions[ *setIter ]));
				tmpStrSet->insert(*setIter);
			}
		}
	}
	y_ion_set.clear();
	y_ion_set = *tmpStrSet;
	delete(tmpStrSet);
	y_ions.clear();
	y_ions = *tmpIonMap;
	delete(tmpIonMap);


	if(g_DEBUG_MODE) {
		cerr << "\nSite determining B-ions (if any):\n";
		for(setIter = b_ion_set.begin(); setIter != b_ion_set.end(); setIter++) {
			cerr << *setIter << endl;
		}

		cerr << "\nSite determining Y-ions (if any):\n";
		for(setIter = y_ion_set.begin(); setIter != y_ion_set.end(); setIter++) {
			cerr << *setIter << endl;
		}
		cerr << endl;
	}


	ret = (signed) y_ion_set.size() + (signed) b_ion_set.size();
	return ret;
}



// Function returns the number of modified STY residues in the given string
int AscoreClass::getNumModSTYs(string txt) {
	int ret = 0;
	size_t found;
	string modChars = "sty234567890@#$%&;?~";

	int N = (signed)txt.length();

	for(int i = 0; i < N; i++) {
		char c = txt.at(i);
		found = modChars.find(c);
		if(found != string::npos) ret++;
	}

	return ret;
}



// Function returns the probability 'pr' for the binomial distribution
// as defined in the PhosphoRS paper
double AscoreClass::getPhosphoRS_pr() {
	double ret = 0;
	double delta_mz;

	double N = (signed) local_spectrum.size();

	delta_mz = maxMZ - minMZ;

	ret = (N * g_MZ_ERR) / delta_mz;

	return ret;
}



