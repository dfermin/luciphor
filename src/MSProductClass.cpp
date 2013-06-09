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

	if(isDecoyPep(&seq) && g_IS_HCD ) {
		mz_err = g_DECOY_MZ_ERR * 0.5;
		if(g_usePPM) mz_err = g_DECOY_MZ_ERR * PPM;
	}

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
		ret += -H3PO4;
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
	string phosphoChars = "sty";
	size_t f;

	bool hasPhospho = false;
	bool loseNH3 = false;
	bool loseH2O = false;


	// In CID data there is sufficient noise in the spectra to allow for matches
	// to decoy fragment ions. This is not the case for HCD data.
	// Therefore we need to allow neutral losses in HCD data to be generated
	// using all of the decoy residues
	if( g_IS_HCD ) phosphoChars = "sty234567890@#$%&;?~";


	// compute mass of ion
	mass = 0.0;
	for(int i = 0; i < N; i++) { mass += AAmass[ ion.at(i) ]; }


	// determine if the ion contains an STY letter, if it does, then the loss
	// of a phospho-group is possible.
	int H3PO4ctr = 0;
	for(int i = 0; i < N; i++) {
		f = phosphoChars.find( ion.at(i) );
		if(f != string::npos) H3PO4ctr++;
	}
	if( H3PO4ctr > 0 ) hasPhospho = true; // has at least 1 phosphorylated AA


/*****************************************************************************
	// determine if the ion can lose ammonia or water
	int NH3ctr = 0, H2Octr = 0;
	for(int i = 0; i < N; i++) {
		if(ion.at(i) == 'R') NH3ctr++;
		if(ion.at(i) == 'K') NH3ctr++;
		if(ion.at(i) == 'Q') NH3ctr++;
		if(ion.at(i) == 'N') NH3ctr++;

		if(ion.at(i) == 'S') H2Octr++;
		if(ion.at(i) == 'T') H2Octr++;
		if(ion.at(i) == 'E') H2Octr++;
		if(ion.at(i) == 'D') H2Octr++;
	}

	if( NH3ctr > 0 ) loseNH3 = true;
	if( H2Octr > 0 ) loseH2O = true;


	// For HCD data, neutral loss peaks help, but they don't help much
	// for CID data. Therefore, don't generate the H2O and NH3 neutral loss
	// peaks for CID data
	if( !g_IS_HCD ) {
		loseNH3 = false;
		loseH2O = false;
	}
*******************************************************************************/


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
				mz_value += -H3PO4;
				mz_value /= z;

				if( (mz_value > MIN_MZ) && (N > 1) ) {
					b_ions[ new_ion ] = mz_value;
					b_ion_set.insert(new_ion);
				}
			} // end if(hasPhospho)


			if(loseNH3) {
				new_ion.clear();
				new_ion = "b^" + Nstr + ":" + ion + "-NH3";
				if(z > 1) new_ion += "/+" + int2string(z);

				mz_value = mass + H -e + extraProton;
				mz_value += -NH3;
				mz_value /= z;

				if( (mz_value > MIN_MZ) && (N > 1) ) {
					b_ions[ new_ion ] = mz_value;
					b_ion_set.insert(new_ion);
				}
			} // end if(loseNH3)


			if(loseH2O) {
				new_ion.clear();
				new_ion = "b^" + Nstr + ":" + ion + "-H2O";
				if(z > 1) new_ion += "/+" + int2string(z);

				mz_value = mass + H -e + extraProton;
				mz_value += -H2O;
				mz_value /= z;

				if( (mz_value > MIN_MZ) && (N > 1) ) {
					b_ions[ new_ion ] = mz_value;
					b_ion_set.insert(new_ion);
				}
			} // end if(loseH2O)

		} // end b-ion
		else if(ion_type == 'y') {

			// this ion contains a phospho group that can undergo neutral loss
			if(hasPhospho) {
				new_ion.clear();
				new_ion = "y^" + Nstr + ":" + ion + "-H3PO4";
				if(z > 1) new_ion += "/+" + int2string(z);

				mz_value = mass + H2O + H + extraProton;
				mz_value += -H3PO4;
				mz_value /= z;

				if( (mz_value > MIN_MZ) && (N > 1) ) {
					y_ions[ new_ion ] = mz_value;
					y_ion_set.insert(new_ion);
				}
			} // end if(hasPhospho)


			if(loseNH3) {
				new_ion.clear();
				new_ion = "y^" + Nstr + ":" + ion + "-NH3";
				if(z > 1) new_ion += "/+" + int2string(z);

				mz_value = mass + H2O + H + extraProton;
				mz_value += -NH3;
				mz_value /= z;

				if( (mz_value > MIN_MZ) && (N > 1) ) {
					y_ions[ new_ion ] = mz_value;
					y_ion_set.insert(new_ion);
				}
			} // end if(loseNH3)


			if(loseH2O) {
				new_ion.clear();
				new_ion = "y^" + Nstr + ":" + ion + "-H2O";
				if(z > 1) new_ion += "/+" + int2string(z);

				mz_value = mass + H2O + H + extraProton;
				mz_value += -H2O;
				mz_value /= z;

				if( (mz_value > MIN_MZ) && (N > 1) ) {
					y_ions[ new_ion ] = mz_value;
					y_ion_set.insert(new_ion);
				}
			} // end if(loseH2O)

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

	ofstream winF;
	if(g_DEBUG_MODE == 4) {
		if(fileExists("mzErrWinMatched.txt") ) winF.open("mzErrWinMatched.txt", ios::app);
		else {
			winF.open("mzErrWinMatched.txt", ios::out);
			winF << "specId\t"
				 << "theoMZ\t"
				 << "obsMZ\t"
				 << "diff\n";
		}
	}

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

					if(g_DEBUG_MODE == 4) {
						winF << specId << "\t"
							 << theo_mz << "\t"
							 << mz << "\t"
							 << (theo_mz - mz)
							 << endl;
					}
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
	if(g_DEBUG_MODE == 4) winF.close();
}



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



// Function returns the m/z values for all of the theoretical peaks for this peptide perumtation
void MSProductClass::assignFragmentIonsMZ(list<double> &ret) {
	map<string, double>::iterator ion;

	for(ion = b_ions.begin(); ion != b_ions.end(); ion++) ret.push_back(ion->second);
	for(ion = y_ions.begin(); ion != y_ions.end(); ion++) ret.push_back(ion->second);
}



// Function extracts the top N peaks from the spectrum and scores them using the
// parameters previously computed. We use the BOOST BIMAP to go between m/z and intensity values
scoreStruct MSProductClass::scorePermutation() {

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

	// Now identify which observed peaks can be matched to theoretical peaks
	M.clear(); // matched peaks
	getMatchedPeaks(&M);

	score = 0.0;
	if(g_IS_HCD) score = calcSpectrumScore_HCD(&M);
	else score = calcSpectrumScore(&M);

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
			peakPtr->MZdistance = ( D.back() );
			peakPtr->ionType = 'u';
			Uptr->insert(pair<double, peakStruct>(mz, *peakPtr));

			delete(peakPtr); peakPtr = NULL;
			D.clear();
		}
	}
}


// function actually does the scoring of the map pointed to by specPtr
double MSProductClass::calcSpectrumScore(map<double, peakStruct> *Mpeaks) {
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

	modelParamStruct *paramPtr = new modelParamStruct;
	*paramPtr = g_modelParamsMap_CID[ charge ];


	// variables for unmatched peaks
	muU = paramPtr->unMatched_mean;
	varU = paramPtr->unMatched_var;
	muU_dist = paramPtr->unMatched_dist_mean;
	varU_dist = paramPtr->unMatched_dist_var;

	N = (signed)Mpeaks->size();

	if(g_DEBUG_MODE == 2) {
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
			if(g_DEBUG_MODE == 2) {
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

	if(g_DEBUG_MODE == 2) debug_ionScores.close();

	delete(paramPtr); paramPtr = NULL;

	return score;
}




// function actually does the scoring of the map pointed to by specPtr
// This function is specifically for HCD data which has a different distribution
// from CID data.
double MSProductClass::calcSpectrumScore_HCD(map<double, peakStruct> *Mpeaks) {

	map<double, peakStruct>::iterator curPeak;

	// intensity variables
	double muM_ints, muU_ints, varM_ints, varU_ints, log_int_M, log_int_U;

	// distance variables
	double muM_dist, varM_dist, varM_dist_IQR, log_dist_M, log_dist_U;

//	double muU_dist, varU_dist;

	double mz, intensity, mzDist, score, Iscore, Dscore, x;

	double pi = 0.60; // fixed for this function call
	int N = (signed)Mpeaks->size();
	bool isNLpeak;
	char peakType; // b, y
	string ionSeq, ss;

	 modelParamStruct *paramPtr = &g_modelParams_HCD;


	if(g_DEBUG_MODE == 2) {
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
//		muU_dist = paramPtr->unMatched_dist_mean;
//		varU_dist = paramPtr->unMatched_dist_var;


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
//			log_dist_U = getLogNPdensityDist_U(mzDist, paramPtr); // what we had 2013 June 6
			log_dist_U = 0; // log of Uniform distribution between -1 to 1 is zero

			Dscore = log_dist_M - log_dist_U;


			// scoring with both intensity and m/z distances equally weighted
			if(dbl_isnan(Iscore) || isInfinite(Iscore) ) Iscore = 0;
			if(dbl_isnan(Dscore) || isInfinite(Dscore) ) Dscore = 0;

			x = Iscore + Dscore;  // official HCD scoring method

			score += x;



			/***************************************************/
			/***************************************************/
			/*  This is where we print out scores for each ion */
			/***************************************************/
			/***************************************************/
			if(g_DEBUG_MODE == 2) {
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

	if(g_DEBUG_MODE == 2) debug_ionScores.close();

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




// Function retains only the site-determining ions found in passed set object
// for the sequence currently assigned to this PSM
void MSProductClass::keepOnlySiteDetermIons(set<string> &ions) {
	set<string>::iterator sdIter; //site-determining iter
	map<string, double> keepers;
	map<string, double>::iterator m;
	string tmp;
	int f;

	// b-ions
	for(m = b_ions.begin(); m != b_ions.end(); m++) {
		f = m->first.find_first_of("/");

		if(f != string::npos) tmp = m->first.substr(0,f);
		else tmp = m->first;

		if( !containsSTY(tmp) ) continue;

		sdIter = ions.find(tmp);
		if(sdIter != ions.end()) keepers[ m->first ] = m->second;
	}


	// y-ions
	for(m = y_ions.begin(); m != y_ions.end(); m++) {
		f = m->first.find_first_of("/");

		if(f != string::npos) tmp = m->first.substr(0,f);
		else tmp = m->first;

		if( !containsSTY(tmp) ) continue;

		sdIter = ions.find(tmp);
		if(sdIter != ions.end()) keepers[ m->first ] = m->second;
	}

	// keep only the ions in the 'keepers' object
	b_ions.clear();
	y_ions.clear();

	for(m = keepers.begin(); m != keepers.end(); m++) {
		if(m->first.at(0) == 'b') b_ions[ m->first ] = m->second;

		if(m->first.at(0) == 'y') y_ions[ m->first ] = m->second;
	}
}



