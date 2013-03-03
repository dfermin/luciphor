/*
 * MSProductClass.cpp
 *
 *  Created on: Apr 16, 2011
 *      Author: dfermin
 */



#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "MSProductClass.hpp"
#include "nonParamDist.hpp"
#include "statsFunctions.hpp"
#include "globals.hpp"

#include <boost/regex.hpp>
#include "boost/bimap.hpp"

using namespace std;

const double H = 1.00727;
const double OH = 17.002745;
const double e = 0.00054858026;
const double H2O = 18.010015;



// default constructor
MSProductClass::MSProductClass(string the_specId, string txt, int Z, double ntm) {

	seq = txt;
	charge = Z;
	min_I = 0;
	nterm_mass = ntm;
	totNumIons = 0;
	frac_matched = 0;
	specId = the_specId;

	mz_err = g_MZ_ERR * 0.5;
	if(g_usePPM) mz_err = g_MZ_ERR * PPM;


	neutralLossMap[ "-H2O" ] =   -1.0 * H2O;
	neutralLossMap[ "+H2O" ] = H2O;
	neutralLossMap[ "-H3PO4" ] = -97.976895;
	neutralLossMap[ "-HPO3" ] =  -79.966330;
	neutralLossMap[ "-NH3" ] =   -17.026549;


	seq_mass = getMass() + nterm_mass;
	makeIons(); // fragment the given sequence into B and Y ions
	totNumIons = ((signed)b_ions.size()) + ((signed)y_ions.size());
}


// Function returns the mass of the current peptide string recorded in 'seq'
double MSProductClass::getMass() {
	double ret = H + H2O;
	int N = (signed)seq.length();

	for(int i = 0; i < N; i++) {
		ret += AAmass[ seq.at(i) ];
	}

	return(ret);
}



// Function returns mass of the given fragment ion
double MSProductClass::getIonMass(string srcStr) {
	double ret = 0;
	int N = 0;
	double z = 1;
	string tmp, ionStr, zStr;
	size_t f1, f2, f3;

	f1 = srcStr.find(":") + 1;
	ionStr = srcStr.substr(f1);

	f2 = ionStr.find("/+");
	if(f2 != string::npos) {
		tmp = ionStr.substr( 0, f2 );
		zStr = ionStr.substr( (f2+2) );
		z = str2dbl(zStr);
		ionStr = tmp;
	}

	f3 = ionStr.find("-");
	if(f3 != string::npos) {
		tmp = ionStr.substr(0, f3);
		ionStr = tmp;
		ret += neutralLossMap[ "-H3PO4" ];
	}




	N = (signed)ionStr.size();

	if(srcStr.at(0) == 'b') { // b-ion
		ret += nterm_mass + (H*(z-1)) + H;
	}
	else { // y-ion
		ret += nterm_mass + (H * (z-1)) + H + H2O;
	}

	for(int i = 0; i < N; i++) ret += AAmass[ ionStr.at(i) ];

	return ret;
}




// Function fragments the sequence and stores its ions
void MSProductClass::makeIons() {
	string b, y;

	int N = (signed)seq.length();

	for(int i = 0; i < N; i++) {
		b = seq.substr(0, (i+1) );
		y = seq.substr( (i+1) );

		if( (signed)b.length() == N ) continue;
		if( (signed)y.length() == N ) continue;

		generateIonsMZ(b, 'b');
		if( !g_NO_NL_PEAKS ) generate_NL_ionsMZ(b, 'b');

		generateIonsMZ(y, 'y');
		if( !g_NO_NL_PEAKS ) generate_NL_ionsMZ(y, 'y');
	}
}


// Function analyzes the given ion string to determine if it can undergo a
// neutral loss of any kind. If it can, that potential neutral loss is recorded.
void MSProductClass::generate_NL_ionsMZ(string ion, char ion_type) {
	double mass, mz_value;
	int N = (signed) ion.length();
	string Nstr = int2string(N);
	string new_ion;
	bool hasPhospho = false;

	int S, T, Y;

	// compute mass of ion
	mass = 0.0;
	for(int i = 0; i < N; i++) { mass += AAmass[ ion.at(i) ]; }

	// determine if the ion contains an STY letter, if it does, then the loss
	// of a phospho-group is possible. Otherwise, the ion can only lose
	// NH3 or H2O/H
	S = 0; T = 0; Y = 0; // count the number of STY in the ion
	for(int i = 0; i < N; i++) {
		if( ion.at(i) == 's') S++;
		if( ion.at(i) == 't') T++;
		if( ion.at(i) == 'y') Y++;

	}
	if( (S+T+Y) > 0 ) hasPhospho = true; // has at least 1 phosphorylated AA

	double extraProton = 0; // Computes: ( H * (z-1) ) <- you need this for dealing with multiple charge states;
	for(int z = 1; z < charge; z++) {

		extraProton = ( H * (z-1) );

		if(ion_type == 'b') {

			// this ion contains a phospho group that can undergo neutral loss
			if(hasPhospho) {
				new_ion.clear();
				new_ion = "b^" + Nstr + ":" + ion + "-H3PO4";
				if(z > 1) new_ion += "/+" + int2string(z);

				mz_value = mass + H - e + extraProton;
				mz_value += neutralLossMap[ "-H3PO4" ];
				mz_value /= z;

				// record the -H3PO4 neutral loss
				if( (mz_value > MIN_MZ) && (N > 1) ) {
					b_ions[ new_ion ] = mz_value;
					b_ion_set.insert(new_ion);
				}
			} // end if(hasPhospho)
		} // end b-ion
		else if(ion_type == 'y') {

			// this ion contains a phospho group that can undergo neutral loss
			if(hasPhospho) {
				new_ion.clear();
				new_ion = "y^" + Nstr + ":" + ion + "-H3PO4";
				if(z > 1) new_ion += "/+" + int2string(z);

				mz_value = mass + H2O + H + extraProton;
				mz_value += neutralLossMap[ "-H3PO4" ];
				mz_value /= z;

				// record the -H3PO4 neutral loss
				if( (mz_value > MIN_MZ) && (N > 1) ) {
					y_ions[ new_ion ] = mz_value;
					y_ion_set.insert(new_ion);
				}
			} // end if(hasPhospho)
		} // end y-ion
	}
}


// function generates the theoretical m/z value for the given ion string.
// All ion variations (neutral losses, etc..) are stored in the b_ions and y_ions
// maps of the MSProductClass object
void MSProductClass::generateIonsMZ(string ion, char ion_type) {
	double mass;
	double mz_value; // the m/z value of the current ion string
	int N = (signed)ion.length();
	string Nstr = int2string(N);
	string new_ion;

	// The code '( H * (z-1) )' is used to add a proton to the ion
	// for multiple charge states. For every additional increase in an ion's
	// charge state, you need to add 1 proton. That additional proton is where
	// the charge state comes from.


	for(int z = 1; z < (charge); z++) {
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



// Function just copies the contents of the passed spectrum map to the
// local_spectrum variable
bool MSProductClass::assignSpectrumMap(map<double, double> src) {
	bool ret = false;
	local_spectrum = src;

	if(local_spectrum.empty()) ret = true;

	return ret;
}





// Function returns an intensity threshold based upon the bottom 10% of the peaks
// for the spectrum in 'local_spectrum'
void MSProductClass::findMinIntensityThreshold() {
	map<double, double>::iterator curPeak;
	list<double> I_list;
	list<double>::iterator curI;
	double intensity;
	int N = 0, tenPercent;

	for(curPeak = local_spectrum.begin(); curPeak != local_spectrum.end(); curPeak++) {
		intensity = curPeak->second;
		I_list.push_back(intensity);
	}

	I_list.sort(); // sorted low-to-high
	N = (signed) I_list.size();

	tenPercent = (int) ((double)N * 0.05);

	int i = 0;
	curI = I_list.begin();
	while(i < tenPercent) {
		curI++;
		i++;
	}
	min_I = *curI;

}



// Function tries to match the theoretical peaks in the b_ion and the y_ion maps
// to the observed spectrum stored in local_spectrum. The matched peaks are stored in
// the matchPtr map pointer that is passed to the function. This function is used
// only for building the model parameters.
void MSProductClass::recordMatchPeaks(bool forModeling ) {

	map<double, double>::iterator curObsPeak, matchIter;
	map<string, double>::iterator theoPeak;
	map<string, double> *ionPtr = NULL;
	map<double, peakStruct>::iterator mIter;
	double a, b, mz, intensity;
	double theo_mz, err_tol;
	double tmpD;
	string ionSeq;
	size_t found;
	list<double> I;
	peakStruct *peakPtr = NULL;
	bool isNLpeak = false;

	// define the bimap type we will use
	// For our BIMAP left bimap: k = mz , v = intensity
	//              right bimap: k = intensity, v = mz
	typedef boost::bimap<double,double> bm_type;
	bm_type bm;
	bm_type::right_const_iterator r_iter;

	matchedPeaks.clear();

	for(int iter = 0; iter < 2; iter++) {

		if(iter == 0) ionPtr = &b_ions;
		else ionPtr = &y_ions;

		for(theoPeak = ionPtr->begin(); theoPeak != ionPtr->end(); theoPeak++) {
			ionSeq = theoPeak->first;
			theo_mz = theoPeak->second;

			a = theo_mz - mz_err;
			b = theo_mz + mz_err;

			if(g_usePPM) {
				a = theo_mz - ppmErrMap[ ionSeq ];
				b = theo_mz + ppmErrMap[ ionSeq ];
			}

			bm.clear();
			I.clear();

			// check to see if this theoretical peak corresponds to a neutral loss
			// peak or not
			isNLpeak = false;
			found = ionSeq.find("-"); // the presence of a dash implies a neutral
									  // loss is annotated on the string
			if(found != string::npos) isNLpeak = true;

			if( forModeling && g_NO_NL_PEAKS_MODEL && isNLpeak ) continue;
			else if( !forModeling && g_NL_MODEL_ONLY && isNLpeak ) continue;


			for(curObsPeak = local_spectrum.begin(); curObsPeak != local_spectrum.end(); curObsPeak++) {
				mz = curObsPeak->first;
				intensity = curObsPeak->second;

				if( (mz >= a) && (mz <= b) ) { // record match
					bm.insert( bm_type::value_type(mz, intensity) );
					I.push_back(intensity);
				}
			}

			/*
			 * This is specific to model parameter acquisition. We don't do any
			 * of what follows for actually scoring peaks.
			 */
			if( !I.empty() ) { // you've got at least one candidate peak
				I.sort();
				intensity = I.back(); // take most intense match
				r_iter = bm.right.find( intensity );
				mz = r_iter->second;

				tmpD = (mz - theo_mz);

				if(g_usePPM) { // scale distances to PPM units instead of Da.
					          // this is an approximation but it works fairly well
					double x = ( tmpD / theo_mz ) / PPM ;
					tmpD = x;
				}


				peakPtr = new peakStruct;
				peakPtr->MZdistance = tmpD;
				peakPtr->intensity = intensity;
				peakPtr->ionType = (iter == 0 ? 'b' : 'y');
				peakPtr->hasSTY =  containsSTY(ionSeq);
				peakPtr->ionStr = ionSeq;

				matchedPeaks[ mz ] = *peakPtr;

				delete(peakPtr); peakPtr = NULL;
			}

		}
	} // end for loop over iter
}



/*********************************
// Function records the unmatched peaks and their distances from the nearest
// matched peak. This function is only used to get the modeling parameters
void MSProductClass::recordUnmatchedPeaks() {

	map<double, peakStruct>::iterator matchedPeak, mIter;
	map<double, double>::iterator curPeak;
	double mz, intensity;
	peakStruct *peakPtr = NULL;

	list<double> distance_list;
	double minDist, absDist;
	double tmpDist;

	multimap<double, double> distMultiMap;
	multimap<double, double>::iterator mm;

	unmatchedPeaks.clear();

	for(curPeak = local_spectrum.begin(); curPeak != local_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second;

		matchedPeak = matchedPeaks.find(mz);
		if(matchedPeak == matchedPeaks.end()) { // unmatched peak
			minDist = 0;
			distance_list.clear();

			// record the distances (along the m/z values) of this unmatched peak
			// from all of the ions in 'matchedPeaks' map. Store the distances into
			// the distance_list.
			for(mIter = matchedPeaks.begin(); mIter != matchedPeaks.end(); mIter++) {
				minDist = mIter->first - mz;
				if( dbl_isnan(minDist) ) minDist = BIG_NUM;
				absDist = fabs( minDist );

				distance_list.push_back( absDist );
				distMultiMap.insert(pair<double, double>(absDist, minDist));
			}


			// From here down, the function picks the minimum distance of this
			// unmatched peak to the nearest matched peak for this permutation

			distance_list.sort(); // sorted low-to-high
			mm = distMultiMap.find( distance_list.front() );
			absDist = (double) mm->first;

			distance_list.clear();
			for(mm = distMultiMap.equal_range(absDist).first; mm != distMultiMap.equal_range(absDist).second; mm++) {
				distance_list.push_back( (*mm).second );
			}
			distance_list.unique();
			distance_list.sort();

			peakPtr = new peakStruct;
			peakPtr->intensity = intensity;
			tmpDist = distance_list.back();
			peakPtr->MZdistance = mz_adjust_dist(    distance_list.back()   );  //////HW: 
			peakPtr->ionType = 'u';

			if( fabs(tmpDist) > 5 && fabs(peakPtr->MZdistance) < 0.25 ) unmatchedPeaks[ mz ] = *peakPtr;

			delete(peakPtr); peakPtr = NULL;
		}
	}
	distance_list.clear();
}
***************************/



// Function returns a map of the peaks in the passed map that are matched
// to the theoretical peaks stored in the b_ions and y_ions maps. This function is
// different from matchPeaks, since it is used for specific phospho-peptide permutations
// while recordMatchPeaks() is used for model parameter building
//void MSProductClass::getMatchedPeaks(map<double, double> *srcMapPtr, map<double, peakStruct> *Mptr) {
void MSProductClass::getMatchedPeaks(map<double, peakStruct> *Mptr) {

	map<double, double>::iterator curObsPeak;
	map<double, peakStruct>::iterator matchIter;
	map<string, double>::iterator theoPeak;
	map<string, double> *ionPtr = NULL;
	double a, b, mz, intensity, theo_mz;
	double tmpD;
	bool isNLpeak;
	peakStruct *X = NULL;
	size_t found;
	string ionSeq, ss;

	// For our BIMAP left bimap: k = mz , v = intensity
	//              right bimap: k = intensity, v = mz
	typedef boost::bimap<double,double> bm_type; // define the bimap type we will use
	bm_type bm;
	bm_type::right_const_iterator r_iter;

	list<double> I;


	for(int iter = 0; iter < 2; iter++) {

		if(iter == 0) ionPtr = &b_ions;
		else ionPtr = &y_ions;

		bm.clear();
		I.clear();
		for(theoPeak = ionPtr->begin(); theoPeak != ionPtr->end(); theoPeak++) {

			ionSeq = theoPeak->first;
			theo_mz = theoPeak->second;

			a = theo_mz - mz_err;
			b = theo_mz + mz_err;

			if(g_usePPM) {
				a = theo_mz - ppmErrMap[ ionSeq ];
				b = theo_mz + ppmErrMap[ ionSeq ];
			}

			// determine if this peak is a neutral loss peak
			isNLpeak = false;
			found = ionSeq.find("-");
			if(found != string::npos) isNLpeak = true;

			// determine if the peak should be used for scoring
			if(isNLpeak && g_NL_MODEL_ONLY) continue;


			// determine if 'ionSeq' contains a phosphorylated AA and if it should be matched or not
			found = ionSeq.find(":") + 1;
			ss = "";
			ss = ionSeq.substr(found);

			bm.clear();
			I.clear();

			for(curObsPeak = local_spectrum.begin(); curObsPeak != local_spectrum.end(); curObsPeak++) {
				mz = curObsPeak->first;
				intensity = curObsPeak->second;

				if( (mz >= a) && (mz <= b) ) { // record match into Mptr
					bm.insert( bm_type::value_type(mz, intensity) );
					I.push_back(intensity);
				}
			}

			if( !I.empty() ) {
				I.sort();
				intensity = I.back(); // take most intense match
				r_iter = bm.right.find( intensity );
				mz = r_iter->second;

				tmpD = (mz - theoPeak->second);
				if(g_usePPM) { // scale distances to PPM units instead of Da.
					           // this is an approximation but it works fairly well.
					double x = ( tmpD / theo_mz ) / PPM;
					tmpD = x;
				}

				X = new peakStruct;
				X->MZdistance = tmpD;
				X->intensity = intensity;
				X->ionType = (iter == 0 ? 'b' : 'y');
				X->ionStr = theoPeak->first;

				matchIter = Mptr->find(mz);
				if(matchIter == Mptr->end()) Mptr->insert(pair<double, peakStruct>(mz, *X));
				else {
					if(matchIter->second.intensity < intensity) matchIter->second = *X;
				}
				delete(X); X = NULL;
			}
		}
	} // end for loop over iter

	// record the fraction of the theoretical ions that were matched for this peptide
	double N = (double) Mptr->size();
	frac_matched = N / (double) totNumIons;
}



// Function APPENDS the data in 'matchedPeaks' map into the map that is passed in
void MSProductClass::addPeakData(map<double, peakStruct> *targetPtr, char whichMap) {

	map<double, peakStruct>::iterator curPeak;
	map<double, peakStruct> *ptr = NULL;

	if(whichMap == 'm') ptr = &matchedPeaks;
	else ptr = &unmatchedPeaks;

	for(curPeak = ptr->begin(); curPeak != ptr->end(); curPeak++)
		targetPtr->insert(pair<double, peakStruct>(curPeak->first, curPeak->second));
}





// Function extracts the top N peaks from the spectrum and scores them using the
// parameters previously computed. We use the BOOST BIMAP to go between m/z and intensity values
scoreStruct MSProductClass::scorePermutation(modelParamStruct *paramPtr, string specId) {

	double mz, intensity;
	double score;
	int i, numPeaksTotal;
	scoreStruct *curScore = NULL;
	scoreStruct retScore;
	map<double, double> curTopNpeaks;
	map<double, peakStruct> M;
	map<double, double>::iterator curPeak;

	if(local_spectrum.empty()) {
		cerr << "\nERROR: local_spectrum map is empty!!\n";
		cerr << "Calling function: MSProductClass::scoreTopPeaks()\n\n";
		exit(-1);
	}

	// Now identify which peaks can be matched to theoretical peaks
	M.clear(); // matched peaks
	getMatchedPeaks(&M);

	score = 0.0;
	if(g_IS_HCD) score = calcSpectrumScore_HCD(&M, paramPtr);
	else score = calcSpectrumScore(&M, paramPtr);

	curScore = new scoreStruct();
	curScore->topNpeaksConsidered = i;
	curScore->matchedPeaks = (signed) M.size();
	curScore->unmatchedPeaks = (i - curScore->matchedPeaks);
	curScore->peakFraction = (double) curScore->matchedPeaks / (double) numPeaksTotal;
	curScore->score = score;
	curScore->scoreByPeaks = score * (double) curScore->matchedPeaks;
	curScore->charge = charge;
	curScore->seq = seq;


	retScore = *curScore;
	delete(curScore); curScore = NULL;

	M.clear();
	curTopNpeaks.clear();

	return retScore;
}




// Function returns the unmatched peaks for the current spectrum object
/*
void MSProductClass::getUnmatchedPeaks(map<double, double> *srcMapPtr, map<double, peakStruct> *Mptr, map<double, peakStruct> *Uptr) {
	map<double, double>::iterator curObsPeak;
	map<double, peakStruct>::iterator m_iter;
	peakStruct *peakPtr = NULL;
	list<double> D;
	double mz, intensity, dist;

	for(curObsPeak = srcMapPtr->begin(); curObsPeak != srcMapPtr->end(); curObsPeak++) {
		mz = curObsPeak->first;
		intensity = curObsPeak->second;

		m_iter = Mptr->find(mz);
		if(m_iter == Mptr->end()) { // unmatched peak

			// this loop is for computing the distance to the nearest matched peak
			D.clear();
			for(m_iter = Mptr->begin(); m_iter != Mptr->end(); m_iter++) {

				dist = fabs(m_iter->first - mz);
				D.push_back( dist );
			}
			D.sort(); // sorted low to high
			peakPtr = new peakStruct;
			peakPtr->intensity = intensity;
			peakPtr->MZdistance = mz_adjust_dist( D.front() );      //// WHY LOG: removed log
			peakPtr->ionType= 'u';
			Uptr->insert(pair<double, peakStruct>(mz, *peakPtr));

			delete(peakPtr); peakPtr = NULL;
			D.clear();
		}
	}
}
*/

// Function returns the unmatched peaks for the current spectrum object
void MSProductClass::getUnmatchedPeaks(map<double, double> *srcMapPtr, map<double, peakStruct> *Mptr, map<double, peakStruct> *Uptr) {
	map<double, double>::iterator curObsPeak;
	map<double, peakStruct>::iterator m_iter;
	peakStruct *peakPtr = NULL;
	list<double> D;
	double mz, intensity;
	double absDist, minDist, dist;
	double tmpDist;

	multimap<double, double> distMultiMap;
	multimap<double, double>::iterator mm;

	for(curObsPeak = srcMapPtr->begin(); curObsPeak != srcMapPtr->end(); curObsPeak++) {
		mz = curObsPeak->first;
		intensity = curObsPeak->second;

		m_iter = Mptr->find(mz);
		if(m_iter == Mptr->end()) { // unmatched peak
			minDist = 0;
			D.clear();

			// record the distances (along m/z values) of this unmatched peak
			// from all of the matched peaks in Mptr
			for(m_iter = Mptr->begin(); m_iter != Mptr->end(); m_iter++) {
				minDist = m_iter->first - mz;
				if( dbl_isnan(minDist) ) minDist = BIG_NUM;
				absDist = fabs(minDist);

				D.push_back(absDist);
				distMultiMap.insert(pair<double, double>(absDist, minDist));
			}

			// this code picks the smallest distance observed
			D.sort();
			mm = distMultiMap.find(D.front());

			D.clear();
			for(mm = distMultiMap.equal_range(absDist).first; mm != distMultiMap.equal_range(absDist).second; mm++) {
				D.push_back( (*mm).second );
			}
			D.unique();
			D.sort();

			peakPtr = new peakStruct;
			peakPtr->intensity = intensity;
			tmpDist = D.back();
			// peakPtr->MZdistance = mz_adjust_dist( D.back() );
			peakPtr->MZdistance = ( D.back() );
			peakPtr->ionType = 'u';
			Uptr->insert(pair<double, peakStruct>(mz, *peakPtr));

			delete(peakPtr); peakPtr = NULL;
			D.clear();
		}
	}
}


// function actually does the scoring of the map pointed to by specPtr
double MSProductClass::calcSpectrumScore(map<double, peakStruct> *Mpeaks, modelParamStruct *paramPtr) {
	map<double, peakStruct>::iterator curPeak;
	double muM, muU, varM, varU, muM_dist, varM_dist, muU_dist, varU_dist;
	double mz, intensity, log_prob_M, log_prob_U;
	double mzDist, log_dist_M, log_dist_U;
	double score = 0.0;
	double Iscore, Dscore, x;
	int N = 0;
	bool isNLpeak;
	char peakType; // b, y
	string ionSeq, ss;
	size_t found;

	// variables for unmatched peaks
	muU = paramPtr->unMatched_mean;
	varU = paramPtr->unMatched_var;
	muU_dist = paramPtr->unMatched_dist_mean;
	varU_dist = paramPtr->unMatched_dist_var;

	N = (signed)Mpeaks->size();

	if(g_DEBUG_MODE) {
		// check to see if the debug file already exists, if so, open it for
		// appending. Otherwise, create it.
		if( fileExists("ionScores.debug") ) {
			debug_ionScores.open("ionScores.debug", ios::out | ios::app);
		}
		else { // create file
			debug_ionScores.open("ionScores.debug", ios::out);
			debug_ionScores << "specId" << "\t"
							<< "ionSeq" << "\t"
							<< "mz" << "\t"
							<< "intensity" << "\t"
							<< "mzDist" << "\t"
							<< "intense_wt" << "\t"
							<< "log_ints_M" << "\t"
							<< "log_ints_U" << "\t"
							<< "log_dist_M" << "\t"
							<< "log_dist_U" << "\t"
							<< "Iscore" << "\t"
							<< "Dscore" << "\t"
							<< "score" << endl; // final score for peak
		}
	}



	if(N == 0) score = 0.0; // the spectrum has no matched peaks
	else {

		score = 0.0;
		for(curPeak = Mpeaks->begin(); curPeak != Mpeaks->end(); curPeak++) {
			mz = curPeak->first;
			intensity = curPeak->second.intensity;
			mzDist    = curPeak->second.MZdistance;
			peakType  = curPeak->second.ionType;
			ionSeq    = curPeak->second.ionStr;

			// determine if 'ionSeq' contains a phosphorylated AA and if it should be scored or not
			found = ionSeq.find(":") + 1;
			ss = "";
			ss = ionSeq.substr(found);

			// determine if 'ionSeq' is a neutral loss peak
			isNLpeak = false;
			found = ionSeq.find("-");
			if(found != string::npos) isNLpeak = true;


			log_prob_M = log_prob_U = log_dist_M = log_dist_U = 0.0;
			muM = varM = muM_dist = varM_dist = 0.0;
			Dscore = Iscore = 0.0;

			if(peakType == 'b') {
				muM  = paramPtr->matched_mean_b;
				varM = paramPtr->matched_var_b;
				muM_dist  = paramPtr->matched_dist_mean_b;
				varM_dist = paramPtr->matched_dist_var_b;
			}
			else if(peakType == 'y') {
				muM  = paramPtr->matched_mean_y;
				varM = paramPtr->matched_var_y;
				muM_dist  = paramPtr->matched_dist_mean_y;
				varM_dist = paramPtr->matched_dist_var_y;
			}


			/*
			 * INTENSITY
			 */
			log_prob_M = log_gaussianProb(muM, varM, intensity);
			log_prob_U = log_gaussianProb(muU, varU, intensity);
			Iscore = log_prob_M - log_prob_U;

			/*
			 * DISTANCE
			 */
			log_dist_M = log_gaussianProb(muM_dist, varM_dist, mzDist);
			log_dist_U = log_gaussianProb(muU_dist, varU_dist, mzDist);
			Dscore = log_dist_M - log_dist_U;

			double intense_wt = 1.0 / ( 1.0 + exp(-Iscore) );
			
			if(dbl_isnan(Dscore) || isInfinite(Dscore)) x = 0;
			else {
				x = intense_wt * Dscore;
			}

			score += x;

			/***************************************************/
			/***************************************************/
			/*  This is where we print out scores for each ion */
			/***************************************************/
			/***************************************************/
			if(g_DEBUG_MODE) {
				debug_ionScores << specId << "\t"
								<< ionSeq << "\t"
								<< mz << "\t"
								<< intensity << "\t"
								<< mzDist << "\t"
								<< intense_wt << "\t"
								<< log_prob_M << "\t"
								<< log_prob_U <<  "\t"
								<< log_dist_M << "\t"
								<< log_dist_U << "\t"
								<< Iscore << "\t"
								<< Dscore << "\t"
								<< x << endl; // final score for peak
			}
		}
	}

	if(g_DEBUG_MODE) debug_ionScores.close();

	// if(score < 0) score = 0;
	return score;
}




// function actually does the scoring of the map pointed to by specPtr
// This function is specifically for HCD data which has a different distribution
// from CID data.
double MSProductClass::calcSpectrumScore_HCD(map<double, peakStruct> *Mpeaks, modelParamStruct *paramPtr) {
	map<double, peakStruct>::iterator curPeak;

	// intensity variables
	double muM_ints, muU_ints, varM_ints, varU_ints, log_int_M, log_int_U;

	// distance variables
	double muM_dist, varM_dist, varM_dist_IQR, log_dist_M, log_dist_U;

	double muU_dist, varU_dist;

	double mz, intensity, mzDist, score, Iscore, Dscore, x;

	double pi = 0.60; // fixed for this function call
	int N = (signed)Mpeaks->size();
	bool isNLpeak;
	char peakType; // b, y
	string ionSeq, ss;


	if(g_DEBUG_MODE) {
		// check to see if the debug file already exists, if so, open it for
		// appending. Otherwise, create it.
		if( fileExists("ionScores.debug") ) {
			debug_ionScores.open("ionScores.debug", ios::out | ios::app);
		}
		else { // create file
			debug_ionScores.open("ionScores.debug", ios::out);
			debug_ionScores << "specId" << "\t"
							<< "ionSeq" << "\t"
							<< "mz" << "\t"
							<< "intensity" << "\t"
							<< "mzDist" << "\t"
				            << "log_dist_M\t"
							<< "log_dist_U\t"
							<< "Dscore\t"
							<< "log_ints_M" << "\t"
							<< "log_ints_U" << "\t"
							<< "log_dist_M" << "\t"
							<< "log_dist_U" << "\t"
							<< "Iscore" << "\t"
							<< "Dscore" << "\t"
						 	<< "score\n";
		}
	}


	if(N == 0) score = 0.0; // the spectrum has no matched peaks
	else {

		// unmatched peak intensity  parameters are the same for all ions
		// so assign them here
		muU_ints = paramPtr->unMatched_mean;
		varU_ints = paramPtr->unMatched_var;
		muU_dist = paramPtr->unMatched_dist_mean;
		varU_dist = paramPtr->unMatched_dist_var;


		score = 0.0;
		for(curPeak = Mpeaks->begin(); curPeak != Mpeaks->end(); curPeak++) {
			mz = curPeak->first;
			intensity = curPeak->second.intensity;
			mzDist    = curPeak->second.MZdistance;
			peakType  = curPeak->second.ionType;
			ionSeq    = curPeak->second.ionStr;

			if(peakType == 'b') {
				muM_ints  = paramPtr->matched_mean_b;
				varM_ints = paramPtr->matched_var_b;
				muM_dist  = paramPtr->matched_dist_mean_b;
				varM_dist = paramPtr->matched_dist_var_b;
				varM_dist_IQR = paramPtr->matched_dist_var_IQR_b;
			}
			else if(peakType == 'y') {
				muM_ints  = paramPtr->matched_mean_y;
				varM_ints = paramPtr->matched_var_y;
				muM_dist  = paramPtr->matched_dist_mean_y;
				varM_dist = paramPtr->matched_dist_var_y;
				varM_dist_IQR = paramPtr->matched_dist_var_IQR_y;
			}


			/*
			 * INTENSITY
			 */
			if(peakType == 'b') log_int_M = getLogNPdensityInt_b(intensity, paramPtr);
			if(peakType == 'y') log_int_M = getLogNPdensityInt_y(intensity, paramPtr);

			log_int_U = getLogNPdensityInt_U(intensity, paramPtr);
			Iscore = log_int_M - log_int_U;

			/*
			 * DISTANCE
			 */
			log_dist_M = getLogNPdensityDist(mzDist, paramPtr);
			log_dist_U = getLogNPdensityDist_U(mzDist, paramPtr);

			Dscore = log_dist_M - log_dist_U;


			// scoring with both intensity and m/z distances equally weighted
			if(dbl_isnan(Iscore) || isInfinite(Iscore) ) Iscore = 0;
			if(dbl_isnan(Dscore) || isInfinite(Dscore) ) Dscore = 0;

			x = Iscore + Dscore;

			score += x;

			/******************************************************
			 *
			 * Hyungwon removed the intensity component to the HCD score in this code.
			 * It seems that the intensity component is unnessary for HCD data.
			 *
				if(dbl_isnan(Dscore) || isInfinite(Dscore) ) x = 0;
				else x = Dscore;

				score += x;
			 *
			/*************************************/





			/***************************************************/
			/***************************************************/
			/*  This is where we print out scores for each ion */
			/***************************************************/
			/***************************************************/
			if(g_DEBUG_MODE) {
				debug_ionScores << specId << "\t"
								<< ionSeq << "\t"
								<< mz << "\t"
								<< intensity << "\t"
								<< mzDist << "\t"
								<< log_dist_M << "\t"
								<< log_dist_U << "\t"
								<< Dscore << "\t"
								<< log_int_M << "\t"
								<< log_int_U <<  "\t"
								<< log_dist_M << "\t"
								<< log_dist_U << "\t"
								<< Iscore << "\t"
								<< Dscore << "\t"
								<< x << endl; // score for this ion
			}
		}
	}

	if(g_DEBUG_MODE) debug_ionScores.close();



	// if(score < 0) score = 0;
	return score;
}







// Function computes and records the optimal fragment ion tolerance for each
// theoretical peak
void MSProductClass::calc_ppm_err() {

    map<string, double>::iterator curIon;
    map<string, double> *ionPtr = NULL;
    double theo_mz, ppm_mz;

    for(int iter = 0; iter < 2; iter++) {
		if(iter == 0) ionPtr = &b_ions;
		else ionPtr = &y_ions;

		for(curIon = ionPtr->begin(); curIon != ionPtr->end(); curIon++) {
				theo_mz = curIon->second;
				ppm_mz = (theo_mz * g_MZ_ERR * PPM) * 0.5;
				ppmErrMap[ curIon->first ] = ppm_mz;
		}
    }
}



// Function for debugging that prints out ion list currently stored in
// this MSProductClass object
void MSProductClass::printIons() {

	map<string, double>::iterator iter;

	cout << seq << "/+" << charge << ",  " << seq_mass << " Da. " << endl;

	for(iter = b_ions.begin(); iter != b_ions.end(); iter++) {
		cout << iter->first << "\t" << iter->second << endl;
	}
	cout << endl;

	for(iter = y_ions.begin(); iter != y_ions.end(); iter++) {
		cout << iter->first << "\t" << iter->second << endl;
	}
	cout << endl;
}



// Function generates site-determining ions for the sequence in 'seq'
void MSProductClass::makeSiteDetermIons() {
	int N = (signed)seq.length();
	vector<int> v_b(N,0), v_y(N,0); // initalized to length N with zero's
	vector<string> B, Y;
	int ctr;
	size_t found;
	string modChars = "sty234567890@#$%&;?~";
	string b_ion_str, y_ion_str, tmp, rev_str;
	int numPhos = 0;

	boost::regex ion_regex("^([by].\\d+:\\w+)(-.*)?(\/\+\\d)?");
	boost::smatch matches;


	// compute number of phosphorylation sites in string
	for(int i = 0; i < N; i++) {
		char c = seq.at(i);
		found = modChars.find(c);
		if(found != string::npos) numPhos++;
	}


	ctr = 0;
	for(int i = 0; i < N; i++) {
		char c = seq.at(i);
		found = modChars.find(c);
		if(found != string::npos) ctr++;
		v_b.at(i) = ctr;
	}


	ctr = 0;
	for(int i = N-1; i > -1; i--) {
		char c = seq.at(i);
		found = modChars.find(c);
		if(found != string::npos) ctr++;
		v_y.at(i) = ctr;
	}


	/*
	 * Wherever the v_b and v_y vectors have a number equal to 'numPhos'
	 * for this peptide, that is the ion ladder number for a site-determining ion
	 */
	// B-ions
	for(int i = 1; i < N-1; i++) { // we do not include/keep the last ion (which is the whole peptide)
		if( v_b.at(i) == numPhos ) {
			b_ion_str = seq.substr(0, (i+1) );
			tmp = "b^" + int2string((i+1)) + ":" + b_ion_str;
			B.push_back(tmp);
		}
	}

	// Y-ions (we have to work from the C-term (right-hand-side) of the string)
	for(int i = N-1; i > 0; i--) {
		if(v_y.at(i) == numPhos) {
			int j = N - i;
			y_ion_str = seq.substr(i, j);
			tmp = "y^" + int2string(j) + ":" + y_ion_str;
			Y.push_back(tmp);
		}
	}


	if(g_DEBUG_MODE) {
		cerr << "\nBEFORE site determining ion trimming:\n";
		for(map<string, double>::iterator mm = b_ions.begin(); mm != b_ions.end(); mm++) {
			cerr << mm->first << endl;
		}
		for(map<string, double>::iterator mm = y_ions.begin(); mm != y_ions.end(); mm++) {
			cerr << mm->first << endl;
		}
		cerr << endl;
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
		cerr << "AFTER site determining ion trimming:\n";
		for(map<string, double>::iterator mm = b_ions.begin(); mm != b_ions.end(); mm++) {
			cerr << mm->first << endl;
		}
		for(map<string, double>::iterator mm = y_ions.begin(); mm != y_ions.end(); mm++) {
			cerr << mm->first << endl;
		}
		cerr << endl;
	}
}

