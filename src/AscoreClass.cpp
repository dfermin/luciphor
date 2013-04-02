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

	if(g_DEBUG_MODE) {
		cerr << endl << seq << endl;
		map<string, double>::iterator m;
		for(m = b_ions.begin(); m != b_ions.end(); m++) cerr << m->first << "\t" << m->second << endl;
		for(m = y_ions.begin(); m != y_ions.end(); m++) cerr << m->first << "\t" << m->second << endl;
		cerr << endl;
	}

}

// Function returns the mass of the current peptide string recorded in 'seq'
double AscoreClass::getMass() {
	double ret = H2O;
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
		//generateNL_ionsMZ(b, 'b'); // neutral loss peaks don't help Ascore, so comment out

		generateIonsMZ(y, 'y');
		//generateNL_ionsMZ(y, 'y'); // neutral loss peaks don't help Ascore, so comment out
	}
}



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

				b_ions[ new_ion ] = round_dbl(mz_value, 2);
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

				y_ions[ new_ion ] = round_dbl(mz_value, 2);
				y_ion_set.insert(new_ion);
			}
		}
	}
}


// Function generates neutral loss fragment ions for the given sequence
void AscoreClass::generateNL_ionsMZ(string ion, char ion_type) {
	double mass, mz_value;
	int N = (signed) ion.length();
	string Nstr = int2string(N);
	string new_ion;

	bool hasPhospho = false;
	bool looseNH3 = false;
	bool looseH2O = false;


	int H3PO4ctr = 0;
	for(int i = 0; i < N; i++) {
		if(ion.at(i) == 's') H3PO4ctr++;
		if(ion.at(i) == 't') H3PO4ctr++;
		if(ion.at(i) == 'y') H3PO4ctr++;
	}
	if(H3PO4ctr > 0) hasPhospho = true;


	int NH3ctr = 0;
	int H2Octr = 0;
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

	if( NH3ctr > 0 ) looseNH3 = true;
	if( H2Octr > 0 ) looseH2O = true;



	for(int z = 1; z < (charge + 0); z++) { // was (charge + 1) 2011.12.05
			mass = 0.0;
			mz_value = 0.0;

			for(int i = 0; i < N; i++) {
				mass += AAmass[ ion.at(i) ];
			}

			if(ion_type == 'b') {

				if(hasPhospho) {
					new_ion.clear();
					new_ion = "b^" + Nstr + ":" + ion + "-H3PO4";
					if(z > 1) new_ion += "/+" + int2string(z);

					mz_value = mass + H - e + ( H * (z-1) ) + nterm_mass;
					mz_value += -H3PO4;
					mz_value /= z;

					if( (mz_value > MIN_MZ) && (N > 1) ) {
						b_ions[ new_ion ] = mz_value;
						b_ion_set.insert(new_ion);
					}
				} // end phospho NL


				if(looseNH3) {
					new_ion.clear();
					new_ion = "b^" + Nstr + ":" + ion + "-NH3";
					if(z > 1) new_ion += "/+" + int2string(z);

					mz_value = mass + H - e + ( H * (z-1) ) + nterm_mass;
					mz_value += -NH3;
					mz_value /= z;

					if( (mz_value > MIN_MZ) && (N > 1) ) {
						b_ions[ new_ion ] = mz_value;
						b_ion_set.insert(new_ion);
					}
				} // end looseNH3


				if(looseH2O) {
					new_ion.clear();
					new_ion = "b^" + Nstr + ":" + ion + "-H2O";
					if(z > 1) new_ion += "/+" + int2string(z);

					mz_value = mass + H - e + ( H * (z-1) ) + nterm_mass;
					mz_value += -H2O;
					mz_value /= z;

					if( (mz_value > MIN_MZ) && (N > 1) ) {
						b_ions[ new_ion ] = mz_value;
						b_ion_set.insert(new_ion);
					}
				} // end looseH2O
			} // end b-ion



			if(ion_type == 'y') {

				if(hasPhospho) {
					new_ion.clear();
					new_ion = "y^" + Nstr + ":" + ion + "-H3PO4";
					if(z > 1) new_ion += "/+" + int2string(z);

					mz_value = mass + H2O + H + ( H * (z-1) );
					mz_value += -H3PO4;
					mz_value /= z;

					if( (mz_value > MIN_MZ) && (N > 1)) {
						y_ions[ new_ion ] = round_dbl(mz_value, 2);
						y_ion_set.insert(new_ion);
					}
				} // end hasPhospho


				if(looseNH3) {
					new_ion.clear();
					new_ion = "y^" + Nstr + ":" + ion + "-NH3";
					if(z > 1) new_ion += "/+" + int2string(z);

					mz_value = mass + H2O + H + ( H * (z-1) );
					mz_value += -NH3;
					mz_value /= z;

					if( (mz_value > MIN_MZ) && (N > 1)) {
						y_ions[ new_ion ] = round_dbl(mz_value, 2);
						y_ion_set.insert(new_ion);
					}
				} // end looseNH3


				if(looseH2O) {
					new_ion.clear();
					new_ion = "y^" + Nstr + ":" + ion + "-H2O";
					if(z > 1) new_ion += "/+" + int2string(z);

					mz_value = mass + H2O + H + ( H * (z-1) );
					mz_value += -H2O;
					mz_value /= z;

					if( (mz_value > MIN_MZ) && (N > 1)) {
						y_ions[ new_ion ] = round_dbl(mz_value, 2);
						y_ion_set.insert(new_ion);
					}
				} // end looseH2O
			} //end y-ion
		}
}





// Function records the passe spectrum into the local spectrum map
void AscoreClass::assignSpectrumMap(map<double, vector<double> > *ptr) {

	map<double, vector<double> >::iterator curPeak;
	map<double, double>::iterator m_iter;
	list<double> mzL;
	double mz, intensity, E;
	double x1, x2, remainder;

	minMZ = 0;
	maxMZ = 0;
	local_spectrum.clear();
	maxIntensity = 0;

	for(curPeak = ptr->begin(); curPeak != ptr->end(); curPeak++) {
		mz = round_dbl( curPeak->first, 2 );
		intensity = curPeak->second.at(0);

		if(intensity > maxIntensity) maxIntensity = intensity;

		m_iter = local_spectrum.find(mz);
		if(m_iter != local_spectrum.end()) {
			if(intensity > m_iter->second) local_spectrum[ mz ] = intensity;
		}

		local_spectrum[ mz ] = intensity;

		mzL.push_back(mz);
	}

	mzL.sort();

	x1 = mzL.front();
	remainder = fmod(x1,100.0);
	minMZ = x1 - remainder;

	x2 = mzL.back();
	remainder = fmod(x2,100.0);
	maxMZ = x2 - remainder + 100.0;

	// remove precursor neutral loss peaks
	removePrecursorNL();

	// normalize intensities to range of 0-100
	for(m_iter = local_spectrum.begin(); m_iter != local_spectrum.end(); m_iter++) {
		E = runif() * 1/10000; // random noise to make unique peak intensities
		intensity = ( (m_iter->second / maxIntensity) * 100.0 ) + E;
		m_iter->second = intensity;
	}

}



// Function removes the precursor neutral loss peaks (if any are present)
void AscoreClass::removePrecursorNL() {
	map<double, double> candPeaks;
	map<double, double>::iterator curPeak;

	list<double> L;

	double mz, intensity, a, b;
	double maxI = 0;
	double Z = (double)charge;

	double precursor_mz = (seq_mass + (H*Z)) / Z;
	double precursor_H3PO4_mz = precursor_mz - (H3PO4/Z);
	double precursor_H2O_mz   = precursor_mz - (H2O/Z);

	double *massPtr = NULL;

	for(int i = 0; i < 2; i++) {

		if(i == 0) massPtr = &precursor_H2O_mz;
		else massPtr = &precursor_H3PO4_mz;

		// prepare for next iteration
		candPeaks.clear();
		L.clear();
		maxI = 0;

		for(curPeak = local_spectrum.begin(); curPeak != local_spectrum.end(); curPeak++) {

			mz = curPeak->first;
			intensity = curPeak->second;

			a = *massPtr - mz_err;
			b = *massPtr + mz_err;

			if( (mz >= a) && (mz <= b) ) { // candidate peak for neutral loss
				candPeaks[ intensity ] = mz;
				L.push_back(intensity);
			}
		}

		if( !candPeaks.empty() ) {
			L.sort(); // sorted low to high
			maxI = L.back(); // most intense peak in candPeaks

			mz = candPeaks[ maxI ];
			local_spectrum.erase( mz );
		}
	}


	// need to record the maximum peak intensity in spectrum
	maxIntensity = 0;
	for(curPeak = local_spectrum.begin(); curPeak != local_spectrum.end(); curPeak++) {
		if(curPeak->second > maxIntensity) maxIntensity = curPeak->second;
	}

}




// Function creates an artificial spectrum from local_spectrum consisting of only the
// top 'numPeaks' found within each MZ_WIN of local_spectrum.
// This returned map is what will be used for calculating Ascore
void AscoreClass::getSpectrumAt_X_depth(int numPeaks, map<double, double> *Mptr) {
	double curMZ, nextMZ;
	double mz, intensity;
	map<double, double>::iterator curPeak;
	map<double, double> bestPeaksMap;
	map<double, double> cp_I2M, cp_M2I;  // candidate peaks Intensity-to-M/Z or M/Z-to-Intensity
	map<double, double>::iterator candIter;
	list<double> *Lptr = NULL;
	list<double>::iterator iterL;

	for(curMZ = minMZ; curMZ < (maxMZ + MZ_WIN); curMZ += MZ_WIN) {
		nextMZ = curMZ + MZ_WIN;

		mz = 0;
		intensity = 0;
		cp_M2I.clear();
		cp_I2M.clear();

		Lptr = new list<double>;

		for(curPeak = local_spectrum.begin(); curPeak != local_spectrum.end(); curPeak++) {
			mz = curPeak->first;
			intensity = curPeak->second;

			if( (mz >= curMZ) && (mz < nextMZ) ) {
				Lptr->push_back(intensity);
				cp_M2I[ mz ] = intensity;
				cp_I2M[ intensity ] = mz;
			}
		}

		if( ((signed)cp_M2I.size()) != ((signed)cp_I2M.size()) ) {
			cerr << "ERROR cp_M2I != cp_I2M:  cannot score with Ascore\n";
			exit;
		}

		// pick the top 'i' most intense peaks in this window out of bimap object
		Lptr->sort();
		Lptr->reverse();
		int i = numPeaks;
		iterL = Lptr->begin();
		while(i > 0) {
			intensity = *iterL;

			candIter = cp_I2M.find(intensity);

			if(candIter != cp_I2M.end()) {
				mz = candIter->second;
				Mptr->insert(pair<double, double>(mz, intensity));
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
		expected_mz = round_dbl(curTheoPeak->second, 0);
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
		expected_mz = round_dbl(curTheoPeak->second, 0);
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



// Function returns a deque of the peptide scores for the currently assigned
// phospho permutation at each peak depth
deque<double> AscoreClass::getPeptideScore() {

	map<double, double>::iterator curPeak;
	map<string, double>::iterator curTheoPeak;
	map<double, double> *curPeakMap = NULL;

	deque<double> ret(10,0); // deque initialized to 10 zeros

	double maxScore = 0;
	double expected_mz, a, b, mz;
	set<double> matched_ions;
	double negLogProb = 0.0;
	double prob = 0.0;
	int numMatchedPeaks = 0;
	int bestDepth = 0;

	double depthWt[] = { 0.5, 0.75, 1, 1, 1, 1, 0.75, 0.5, 0.25, 0.25 };

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

			double pr = getPeakProb( (double)numPeaks );

			prob = cum_binomial_prob(N, k, pr);
			negLogProb = ( -10.0 * log(prob) ) * depthWt[ (numPeaks - 1) ];
		}

		ret.at( (numPeaks-1) ) = negLogProb;

		if(negLogProb > maxScore) maxScore = negLogProb;

		delete(curPeakMap); curPeakMap = NULL;
	}


	return(ret);
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

	matched_spectrum.clear();

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

	//////////////////////////////////
	// b-ions
	//////////////////////////////////
	candMatchedPeaks.clear();
	for(curTheoPeak = b_ions.begin(); curTheoPeak != b_ions.end(); curTheoPeak++) {
		expected_mz = round_dbl( curTheoPeak->second, 0 );
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


	///////////////////////////
	// y-ions
	//////////////////////////
	candMatchedPeaks.clear();
	for(curTheoPeak = y_ions.begin(); curTheoPeak != y_ions.end(); curTheoPeak++) {
		expected_mz = round_dbl( curTheoPeak->second, 0 );
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


	matched_ions = (signed) matched_spectrum.size();

	if(matched_ions == 0) ret = 0.0;
	else {
		double k = matched_ions;
		int N = (signed)b_ion_set.size() + (signed)y_ion_set.size();

		double pr = getPeakProb( (double)optimalPeakDepth );

		double prob = cum_binomial_prob(N, k, pr);
		ret = -10.0 * log(prob);
	}

	delete(curPeakMap);
	return ret;
}



// Function retains only the site-determining ions for this peptide permutation
void AscoreClass::recordSiteDetermIons(set<string> &ions) {

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

		if(g_DEBUG_MODE) {
			cerr << "|" << seq << "|\t" << m->first << "\t" << m->second << endl;
		}
	}

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



// Function just returns the probability value to be used in Ascore calculation
// We are using this function because there are 3 competing methods for how to
// compute the probability. It's a lot easier to simply change this function
// than it is to change all instances of the formula throughout the code
double AscoreClass::getPeakProb(double i) {
	double ret = 0;


	if(i == -1) ret = getPhosphoRS_pr();
	else {

		//ret = (i / MZ_WIN) * g_MZ_ERR; // this equation is from Alexey's old python code

		ret = i / MZ_WIN;  // used in original Ascore manuscript
	}

	return ret;
}


