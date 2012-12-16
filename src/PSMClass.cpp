/*
 * PSMClass.cpp
 *
 *  Created on: Mar 18, 2011
 *      Author: dfermin
 */


#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <set>
#include <boost/regex.hpp>
//#include <boost/thread/mutex.hpp>
#include <boost/filesystem/operations.hpp>
#include "PSMClass.hpp"
#include "statsFunctions.hpp"
#include "MSProductClass.hpp"
#include "AscoreClass.hpp"
#include "structs.hpp"
#include "globals.hpp"

using namespace std;
using namespace boost;

//boost::mutex mutex1;


// class ctor that takes in a struct of data values
PSMClass::PSMClass(matchDataStruct *ptr) {

	// this code assigns the modifications to the positions held  in ptr->mods map
	map<int,double>::iterator mod_pos;
	map<char, string>::iterator modAAmap_iter;
	char modLetter = '*';
	char AAchar = 'X';
	double curModMass = 0.0;
	int phosphoCtr = 0; // for counting actual phospho-sites currently reported
	int STYCtr = 0; // for counting potential phospho-sites
	int idx = 0; // position of modification

	is_unambiguous = false;
	use_for_model = false;
	numPhosphoSites = 0;
	numSupportingSpectra = 0;
	nterm_mass = 0.0;
	min_intensity = 0;
	numDecoyPermutations = 0;


	mz_err = g_MZ_ERR * 0.5;
	if(g_usePPM) mz_err *= PPM;


	modPeptide = ptr->peptide;

	matchedPeakMap.clear();

	if( !ptr->mods.empty() ) { // this sequence contains at least 1 modification

		for(mod_pos = ptr->mods.begin(); mod_pos != ptr->mods.end(); mod_pos++) {

			idx = mod_pos->first; // position of the modification in peptide string
			curModMass = round_dbl( mod_pos->second, 2 );

			if(idx != -1) {
				AAchar = tolower( ptr->peptide.at(idx) ); // residue that is reportedly modified
				modAAmap_iter = modAAmap.find(AAchar);

				if(modAAmap_iter != modAAmap.end()) { // you found the modification
					modPeptide[ idx ] = AAchar;
				}
			}
			else { // n-terminal modification found
				nterm_mass = curModMass;
			}
		}
	}


	// check to see if this modPeptide contains at least 2 potential
	// phopsho-sites and at least 1 of them is phosphorylated
	for(int i = 0; i < (signed)modPeptide.length(); i++) {
		if(modPeptide.at(i) == 's') phosphoCtr++;
		if(modPeptide.at(i) == 't') phosphoCtr++;
		if(modPeptide.at(i) == 'y') phosphoCtr++;
	}
	for(int i = 0; i < (signed)ptr->peptide.length(); i++) {
		if(ptr->peptide.at(i) == 'S') STYCtr++;
		if(ptr->peptide.at(i) == 'T') STYCtr++;
		if(ptr->peptide.at(i) == 'Y') STYCtr++;
	}

	// we only keep PSM's where there are multiple potential phospho-sites
	// and at least one of them is annotated as being phosphorylated.
	is_valid_phosphoPSM = false;


	//if( (phosphoCtr > 0) && (STYCtr > 1) && (phosphoDelta > 0) ) { // this is a keeper
	if( (phosphoCtr > 0) && (STYCtr > 0) ) { // this is a keeper
		specId = ptr->specId;
		peptide = ptr->peptide;
		origModPeptide = modPeptide;
		mass = ptr->mass;
		charge = ptr->charge;
		iniProb = ptr->iniProb;
		numPotentialSites = STYCtr;
		numPhosphoSites = phosphoCtr;
		is_unambiguous = false;
		scanNum = ptr->scanNum;

		identify_STY_sites(); // record the positions of the STY characters

		is_valid_phosphoPSM = true;
		if(iniProb >= g_model_prob) use_for_model = true;
	}

	// This means all potential phosphorylation sites are phosphorylated in the peptide
	// there is no actual need to run luciphor on these cases but we score them anyways
	if( (phosphoCtr == STYCtr) ) {
		is_unambiguous = true;
		numPermutations = 1;

	}
	else numPermutations = combinatorial(numPotentialSites, numPhosphoSites);

	if(g_randDecoyAA) calcNumDecoyPermutations(); // record decoy permutation count
}




// Function records the positions of STY on peptide sequence (0-based counting)
void PSMClass::identify_STY_sites() {

	styPos.clear();
	styPos.reserve(numPotentialSites);

	for(int i = 0; i < (signed)peptide.length(); i++) {

		if( (peptide.at(i) == 'S') ||
			(peptide.at(i) == 'T') ||
			(peptide.at(i) == 'Y')
		   ) {
			styPos.push_back(i);
		}
	}
}


// Function returns the name of the spectrum file (mzXML, mzML) for the PSMClass
string PSMClass::getSpectrumFileName() {

	boost::smatch matches; // used to capture REGEX matches, element[0] is the whole string
	boost::regex spectrum_regex("(.*\/)?(.+)\\.\\d+\\.\\d+\\.\\d+$");
	string ret;

	if( boost::regex_match(specId, matches, spectrum_regex) ) {
		ret.assign( matches[2].first, matches[2].second);
		ret += "." + g_ext;
	}
	return ret;
}



// Function records the maximum intensity observed in the raw spectrum
void PSMClass::recordMaxIntensity() {
	map<double, vector<double> >::iterator curPeak;
	list<double> *allI = NULL;

	allI = new list<double>;
	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		allI->push_back( curPeak->second.at(0) ); // record intensities
	}
	allI->sort(); // sorted low to high
	max_intensity = allI->back();

	delete(allI); allI = NULL;
}



// Function records spectrum passed to it from struct.
// It also records the maximum intensity observed from the raw spectrum object
void PSMClass::recordSpectrum(SpecStruct spec) {
	double curMZ, curI;
	int N = (signed) spec.mz.size();
	list<double> *all_intensities = NULL;
	vector<double> *v = NULL;
	max_intensity = 0;
	max_mz = 0;
	//all_intensities = new list<double>();

	for(int i = 0; i < N; i++) {
		curMZ = spec.mz.at(i);
		curI  = spec.intensity.at(i);

		if(curMZ == 0) continue;
		if(curI == 0) continue;
		//all_intensities->push_back( curI );

		if(curMZ > max_mz) max_mz = curMZ;
		if(curI > max_intensity) max_intensity = curI;

		v = new vector<double>;
		v->reserve(3);
		v->assign(3,0);
		v->at(0) = curI; // raw peak intensity
		raw_spectrum[ curMZ ] = *v;
		delete(v); v = NULL;
	}

//	all_intensities->sort(); // sorted low to high
//	max_intensity = all_intensities->back(); // most intense value is last
//
//	delete(all_intensities); all_intensities = NULL;

	deisotopeSpectrum();
	//binPeaks();
	normalizeSpectrum(); // scale spectrum to be in the range of 0-100
	//reduceNeutralLossPeak(); // now that you normalized the spectrum, reduce the impact of the NL peak
	medianNormalizeIntensities(); // divide the intensities by their median value
}




// Function normalizes the data in spectrum map to be in the range of 0-10
void PSMClass::normalizeSpectrum() {
	map<double, vector<double> >::iterator curPeak;
	double mz, curI;

	// now normalize the values in spectrum to max_intensity
	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		curI = curPeak->second.at(0);
		curPeak->second.at(1) = ( curI / max_intensity ) * 100.00;
	}
}



// Function normalizes peak intensities to the median peak intensity in the spectrum
void PSMClass::medianNormalizeIntensities() {
	map<double, vector<double> >::iterator curPeak;
	double mz, intensity;
	list<double> intensity_list;
	list<double>::iterator curI;
	int i, N, idx;
	double part1, part2;

	// get the scaled intensities and store them in a list
	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second.at(1);
		intensity_list.push_back(intensity);
	}

	// sort the intensities low to high
	intensity_list.sort();
	N = (signed) intensity_list.size();

	idx = 0;
	if( (N % 2) == 0 ) { // even number of elements, median value is average of middle 2 values
		idx = (N / 2) - 1;
		curI = intensity_list.begin();
		for(i = 0; i < idx; i++) curI++;

		part1 = *curI; // lower half
		curI++;
		part2 = *curI; // upper half

		median_intensity = (part1 + part2) / 2.0;
	}
	else { // odd number of elements
		idx = (N / 2); // in zero-based coordinates, this is the middle value
		curI = intensity_list.begin();
		for(i = 0; i < idx; i++) curI++;
		median_intensity = *curI;
	}

	// now normalize the peak intensities to the median intensity for the spectrum
	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second.at(1);
		curPeak->second.at(2) = log( intensity / median_intensity );
	}
}




// Function to increase peak intensities by reducing the intensity of the neutral loss peak
void PSMClass::reduceNeutralLossPeak() {

	map<double, vector<double> >::iterator curPeak;

	double mz, intensity;
	double NL_mass, NL_mz, NL_a, NL_b, NL_peak_mz, NL_peak_I;
	double orig_maxIntensity; // used to hold the "official" max. intensity value for this spectrum
	double raw_intensity;

	list<double> *NL_intensity_list = NULL;
	list<double>::iterator curI;

    // For our BIMAP left bimap: k = mz , v = intensity
    //              right bimap: k = intensity, v = mz
    typedef boost::bimap<double,double> bm_type; // define the bimap type we will use
    bm_type bm;
    bm_type::right_const_iterator r_iter;


	NL_mass = mass - 97.976895;
	NL_mz = NL_mass / (double) charge;
	NL_a = NL_mz - mz_err;
	NL_b = NL_mz + mz_err;


	NL_intensity_list = new list<double>;
	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second.at(1); // 0-100 scaled intensities

		if( (mz >= NL_a) && (mz <= NL_b) ) {
			bm.insert( bm_type::value_type(mz, intensity) );
			NL_intensity_list->push_back(intensity);
		}
	}

	if( !NL_intensity_list->empty() ) { // a potential neutral loss peak was found
		NL_intensity_list->sort(); // sorted low to high
		NL_peak_I = NL_intensity_list->back(); // most intense candidate NL peak in list

		r_iter = bm.right.find( NL_peak_I );
		NL_peak_mz = r_iter->second; // get the m/z for this peak

		curPeak = raw_spectrum.find(NL_peak_mz);
		raw_intensity = curPeak->second.at(0);

		//raw_spectrum.erase(curPeak);
		raw_spectrum[NL_peak_mz].at(0) = 0;
		raw_spectrum[NL_peak_mz].at(1) = 0;
		raw_spectrum[NL_peak_mz].at(2) = 0;



		// since we deleted this peak we need to recompute the new max_intensity
		// for the spectrum.
		orig_maxIntensity = max_intensity;
		recordMaxIntensity(); // record the NEW max intensity value for the spectrum

		normalizeSpectrum(); // scale spectrum again to be in the range of 0-100

		//re-insert the NL peak
		raw_spectrum[NL_peak_mz].at(0) = raw_intensity;
		raw_spectrum[NL_peak_mz].at(1) = NL_peak_I;
		raw_spectrum[NL_peak_mz].at(2) = 0;

		max_intensity = orig_maxIntensity; // restore the original max_intensity value
	}
	bm.clear();
	delete(NL_intensity_list); NL_intensity_list = NULL;
}




// Function combines peak intensities into bins fragmented based upon the
// fragment ion tolerance given by the user
void PSMClass::binPeaks() {
	map<double, vector<double> >::iterator curPeak;
	map<double, vector<double> > binMap;
	vector<double> *v = NULL;
	double mz, intensity;
	double i, j;
	double startMZ = MIN_MZ - g_MZ_ERR;
	double endMZ   = max_mz + g_MZ_ERR;
	double C = round_dbl((g_MZ_ERR * 5.0), 0 );


	// initialize new map
	for(mz = startMZ; mz <= endMZ; mz += g_MZ_ERR) {
		v = new vector<double>(3,0);
		binMap[ mz ] = *v;
		delete(v);
	}

	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second.at(0);

		startMZ = round_dbl(mz, 0) - C;
		endMZ   = round_dbl(mz, 0) + C;

		for(i = startMZ; i < endMZ; i += g_MZ_ERR) {
			j = i + g_MZ_ERR;

			if( (mz >= i) && (mz < j) ) {
				binMap[ i ].at(0) += intensity;
				break;
			}
		}
	}

	raw_spectrum.clear();
	raw_spectrum = binMap;
}





// Function identifies the threshold for noisy peaks in the spectrum
void PSMClass::identifyNoisyPeaks() {
	list<double> *I_list = NULL;
	map<double, vector<double> >::iterator curPeak;
	int N = 0;

	peakType = 2; // median normalized

	N = (signed) raw_spectrum.size();

	I_list = new list<double>(0,N);

	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		I_list->push_back(curPeak->second.at(peakType));
	}

	double mu = getMean(I_list);
	double sigma2 = getVar(I_list);
	double sigma = sqrt(sigma2);

	delete(I_list); I_list = NULL;

	// log-scaled QN intensity values that are 3 standard deviations from mean intensity
	// will be considered as noise
	min_intensity = mu - ( 3.0 * sigma );
}




// Function sets which map the 'spectrum' pointer should point to
void PSMClass::setSpectrumPtr(string whichMap) {

	if(whichMap.compare("raw") == 0) peakType = 0;
	else if(whichMap.compare("scaled") == 0) peakType = 1;
	else if(whichMap.compare("median") == 0) peakType = 2;
}




// Function writes the final spectrum scored by Ascore to disk
void PSMClass::writeAscoreSpectrum() {

	// this part figures out what to name the output directory
	string fldr;
	if( !g_userDefinedOutput ) {
		// set the output directory's name based upon input file's name
		boost::filesystem::path curFile( g_srcXMLfile.c_str() );
		string basename = curFile.leaf();
		size_t found = basename.find('.');
		if(found != string::npos) {
			fldr = basename.substr(0, found);
		}
	}
	else {
		size_t found = g_outputName.find('.');
		if(found != string::npos) fldr = g_outputName.substr(0,found);
		else fldr = g_outputName;
	}


	//construct a output directory
	fldr += ".spectra_files_ascore." + getTimeStamp();
	filesystem::path f( fldr.c_str() );
	if( !filesystem::exists(f) ) filesystem::create_directory(f);

	string fileName = specId + ".tsv";
	filesystem::path filePath(f/fileName);

	ofstream out;
	string opf = filePath.file_string();
	out.open( opf.c_str(), ios::out);
	if(!out) {
		cerr << "\nUnable to create '" << opf << "'\n\n";
		exit(-1);
	}

	string ionSeq;
	double mz, intensity;
	double norm_intensity;
	map<double, vector<double> >::iterator curPeak;
	map<double, peakStruct>::iterator ascoreIter;
	peakType = g_intensityType - 1;


	out << "mz\tI\tion\n";

	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second.at(peakType);

		out << mz << "\t" << intensity << "\t";

		ascoreIter = ascoreMatchedSpectrum1.find( round_dbl(mz, 0) );
		if(ascoreIter != ascoreMatchedSpectrum1.end()) {

			if(g_intensityType == 2) {
				if(ascoreIter->second.norm_intensity == intensity) out << ascoreIter->second.ionStr;
			}
			else if(ascoreIter->second.intensity == intensity) out << ascoreIter->second.ionStr;
			//out << ascoreIter->second.ionStr;
		}
		out << endl;
	}
	out.close();
}



// Write spectrum to file
void PSMClass::writeSpectrumToDisk() {

	// this part figures out what to name the output directory
	string fldr;
	if( g_userDefinedOutput ) {
		size_t found = g_outputName.find('.');
		if(found != string::npos) fldr = g_outputName.substr(0,found);
		else fldr = g_outputName;
	}
	else {
		// set the output directory's name based upon input file's name
		boost::filesystem::path curFile( g_srcXMLfile.c_str() );
		string basename = curFile.leaf();
		size_t found = basename.find('.');
		if(found != string::npos) {
			fldr = basename.substr(0, found);
		}
	}


	//construct a output directory
	fldr += ".spectra_files." + getTimeStamp();
	filesystem::path f( fldr.c_str() );
	if( !filesystem::exists(f) ) filesystem::create_directory(f);

	string fileName = specId + ".tsv";
	filesystem::path filePath(f/fileName);

	ofstream out;
	string opf = filePath.file_string();
	out.open( opf.c_str(), ios::out);
	if(!out) {
		cerr << "\nUnable to create '" << opf << "'\n\n";
		exit(-1);
	}

	string ionSeq;
	double mz, intensity;
	double dist, finalScore;
	map<double, vector<double> >::iterator curPeak;
	map<double, peakStruct>::iterator matchIter;
	size_t found;
	string ss;

	peakType = g_intensityType - 1;

	// Need one column header because user elected to print features for both of
	// the predicted peptide permutations
	if(g_WRITE_TOP_TWO) out << "num\t";
	out << "mz\tI\tion\tabs_dist\tIscore\tDscore\tFinalPeakScore\n";

	int stopAt = 1;
	if(g_WRITE_TOP_TWO) stopAt = 2;

	matchedSpectrumStruct *ptr = NULL;
	ptr = &bestSpectrum;
	for(int iter = 0; iter < stopAt; iter++) {

		if(iter > 0) ptr = &nextBestSpectrum;


		for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
			mz = curPeak->first;
			intensity = curPeak->second.at(peakType);

			if(g_WRITE_TOP_TWO) out << iter << "\t";

			out << mz << "\t" << intensity; // default output

			matchIter = ptr->matched_ions.find(mz);
			if(matchIter != ptr->matched_ions.end()) {


				dist = fabs(matchIter->second.MZdistance);
				ionSeq = matchIter->second.ionStr;

				// determine if 'ionSeq' contains a phosphorylated AA and if it should be reported or not
				found = ionSeq.find(":") + 1;
				ss = "";
				ss = ionSeq.substr(found);

				out << "\t" << ionSeq << "\t"
					<< fabs(dist) << "\t"
					<< ptr->IscoreMap[ ionSeq ] << "\t"
					<< ptr->DscoreMap[ ionSeq ] << "\t"
					<< ptr->FinalScoreMap[ ionSeq ] << "\t";
			}
			else out << "\t\t\t\t\t"; // not technically necessary but makes for easier parsing

			out << endl;
		}
	} // end for loop over iter


	out.close();
}



// Function to generate all the possible phospho-peptide permutations from the peptide
// sequence currently held in 'peptide' and the value of 'numSites'
void PSMClass::generatePermutations() {

	/*
	 * Generate each phosphoPeptide permutation and store it into phosphoVersionSet.
	 * The value of each map element will be the score for that sequence using the
	 * stored spectrum.
	 */
	set<vector<int> > bitMat = generateBinaryGrid(numPotentialSites, (double)numPhosphoSites);

	/*
	 * Use the bit matrix (bitMax) to change the case of the phospho-sites in
	 * the peptide sequence. Where ever a bitMat entry is 1,
	 * that site needs to be "phosphorylated" in that peptide model. We store
	 * each version of the phosphroylated peptide in a map.
	 */
	set< vector<int> >::iterator s_iter;

	string *curPepVersion = NULL;
	vector< int > *curBitMap = NULL;

	for(s_iter = bitMat.begin(); s_iter != bitMat.end(); s_iter++) {
		curBitMap = new vector<int>;
		curPepVersion = new string( peptide );

		*curBitMap = *s_iter;
		for(int i = 0; i < (signed) curBitMap->size(); i++) {

			if(curBitMap->at(i) == 1) { // flip this letter's case
				int pos = styPos[i]; // get letter's position in string
				char newLtr = tolower( curPepVersion->at( pos ) ); // flip case
				curPepVersion->at( pos ) = newLtr; // replace CAPS with lower case
			}
		}

		// The above code converts all non-phosphorylated PTMs to CAPS.
		// We need to change the case of all non-phosphorylated PTMs.
		for(int i = 0; i < (signed) origModPeptide.length(); i++) {
			char x = origModPeptide.at(i);
			if( (x == 's') || (x == 'S') ||
				(x == 't') || (x == 'T') ||
				(x == 'y') || (x == 'Y')
			) continue;

			// retain the case of the current non-phospho AA
			curPepVersion->at(i) = x;
		}

		phosphoVersionSet.insert( *curPepVersion );

		delete(curPepVersion);
		delete(curBitMap);
	}

}



// Function classifies all the peaks in the spectrum as being matched or unmatched
void PSMClass::classifyPeaks() {

	set<string>::iterator curPermutation;

	MSProductClass *curMSProduct = NULL;
	map<double, double> *spectrum = NULL;
	map<double, vector<double> >::iterator curPeak;
	double mz, intensity;

	matchedPeakMap.clear();
	unmatchedPeakMap.clear();

	if(raw_spectrum.empty()) return; // If this PSM has no spectrum associated with it skip it.

	// store the relevant peak type into a new map that is passed to the curMSProduct object
	spectrum = new map<double, double>;
	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second.at(peakType);
		spectrum->insert( pair<double, double>(mz, intensity) );
	}

	for(curPermutation = phosphoVersionSet.begin(); curPermutation != phosphoVersionSet.end(); curPermutation++) {

		curMSProduct = new MSProductClass(specId, *curPermutation, charge, nterm_mass);

		curMSProduct->assignSpectrumMap( *spectrum );
		if(g_usePPM) curMSProduct->calc_ppm_err();

		curMSProduct->recordMatchPeaks( use_for_model );

		//curMSProduct->recordMatchPeaks( use_for_model );
		if( use_for_model ) curMSProduct->addPeakData(&matchedPeakMap, 'm');

		delete(curMSProduct); curMSProduct = NULL;
	}


	// This can happen if there were no matched or unmatched peaks
	// identified for the current PSM due to the exclusion of neutral loss peaks.
	// When this happens we will change that modeling status of the PSM to have it
	// discarded in the 'threaded_recordModelingParameters_*' functions
	if( matchedPeakMap.empty() ) use_for_model = false;
	else {
		/*
		 * All of the peaks that *can* be matched for this spectrum in every permutation
		 * are stored in matchedPeakMap. Now we record all of the remaining peaks as unmatched
		 */
		map<double, peakStruct>::iterator curMatchedPeakIter;
		for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
			mz = curPeak->first;
			intensity = curPeak->second.at(peakType);

			curMatchedPeakIter = matchedPeakMap.find(mz);
			if(curMatchedPeakIter == matchedPeakMap.end()) { // this is an unmatched peak in all permutations
				getMinDistance_unmatched(mz, intensity);
			}
		}
	}
	delete(spectrum); spectrum = NULL;

}




// Function records the distance of an unmatched peak to the nearest
// matched peak for the current spectrum based upon all the permutations tried
// on this spectrum
void PSMClass::getMinDistance_unmatched(double mz, double intensity) {

	map<double, peakStruct>::iterator matchedPeak;
	peakStruct *peakPtr = NULL;

	list<double> distance_list;
	double absDist, minDist;

	multimap<double, double> distMultiMap;
	multimap<double, double>::iterator mm;

	for(matchedPeak = matchedPeakMap.begin(); matchedPeak != matchedPeakMap.end(); matchedPeak++) {
		minDist = matchedPeak->first - mz;
		if( dbl_isnan(minDist) ) minDist = TINY_NUM;
		absDist = fabs( minDist );

		distance_list.push_back( absDist );
		distMultiMap.insert(pair<double, double>(absDist, minDist));
	}

	// From here down, the function picks the minimum distance of this
	// unmatched peak to the nearest matched peak for this permutation

	distance_list.sort(); // low to high
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
	peakPtr->MZdistance = distance_list.front();
	peakPtr->ionType = 'u';

	unmatchedPeakMap[ mz ] = *peakPtr;

	delete(peakPtr); peakPtr = NULL;
}





// Function to score each phospho-peptide permutation that is associated with this
// PSMClass object
void PSMClass::calcScore() {
	set<string>::iterator curPermutation;
	set<string> *curSetPtr = NULL;
	time_t start_t, end_t;

	MSProductClass *curMSProduct = NULL;
	scoreStruct *curScore = NULL;
	deque<scoreStruct> scoreDeq; // hold scores for this PSMClass object

	map<double, double> *spectrum = NULL;
	map<double, vector<double> >::iterator curPeak;
	double mz, intensity;
	bool status;


	// store the relevant peak type into a new map that is passed to the curMSProduct object
	spectrum = new map<double, double>;
	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second.at(peakType);
		spectrum->insert( pair<double, double>(mz, intensity) );
	}

	time(&start_t);
	// this loop first scores the forward permutations, then the decoys
	for(int iter = 0; iter < 2; iter++) {

		if(iter == 0) curSetPtr = &phosphoVersionSet;
		else curSetPtr = &decoySet;

		// now score the data pointed to by curSetPtr
		for(curPermutation = curSetPtr->begin(); curPermutation != curSetPtr->end(); curPermutation++) {

			curMSProduct = new MSProductClass(specId, *curPermutation, charge, nterm_mass);
			curMSProduct->updateMZerr(mz_err);

			status = curMSProduct->assignSpectrumMap( *spectrum ); // spectrum should point to QN_spectrum map
			if(g_usePPM) curMSProduct->calc_ppm_err();

			if(g_useOnlySiteDetermIons) curMSProduct->makeSiteDetermIons();

			if(status) {
				cerr << "ERROR: MSProductClass.assignSpectrumMap() spectrum assignment failed.\n"
					 << specId << " has no spectrum map assigned to it.\n"
					 << "Unable to continue. Exiting now.\n\n";
				exit(-1);
			}

			curScore = new scoreStruct();
			*curScore = curMSProduct->scorePermutation(&local_params, specId);

			scoreDeq.push_back( *curScore );

			delete(curScore); curScore = NULL;
			delete(curMSProduct); curMSProduct = NULL;
		}
	}

	delete(spectrum); spectrum = NULL;

	// pick the best score and next-best score that is stored in the scoreVec vector
	pickScores(scoreDeq);
	scoreDeq.clear();
	time(&end_t);

	scoreTime = difftime(end_t, start_t); // how many seconds it took to score this PSM

	// this code is only to provide the user some kind of feed back while the program is running
	//boost::mutex::scoped_lock mutex_locker(mutex1, defer_lock); // defer_lock is initially unlocked
	//mutex_locker.lock();
	g_progressCtr++;
	if(g_progressCtr % 100 == 0 ) cerr << g_progressCtr << " ";
	if(g_progressCtr % 1000 == 0 ) cerr << endl;
	//mutex_locker.unlock();
}



// Function to score each phospho-peptide permutation that is associated with this
// PSMClass object
//void PSMClass::get_rawScores() {
//	set<string>::iterator curPermutation;
//	MSProductClass *curMSProduct = NULL;
//	map<double, double> *spectrum = NULL;
//	map<double, vector<double> >::iterator curPeak;
//	double mz, intensity;
//
//
//	// store the relevant peak type into a new map that is passed to the curMSProduct object
//	spectrum = new map<double, double>;
//	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
//		mz = curPeak->first;
//		intensity = curPeak->second.at(peakType);
//		spectrum->insert( pair<double, double>(mz, intensity) );
//	}
//
//
//	for(curPermutation = phosphoVersionSet.begin(); curPermutation != phosphoVersionSet.end(); curPermutation++) {
//
//		curMSProduct = new MSProductClass(*curPermutation, charge, nterm_mass);
//
//		curMSProduct->assignSpectrumMap( *spectrum ); // spectrum should point to QN_spectrum map
//
//		curMSProduct->scoreTopPeaks(&local_params, specId);
//
//		delete(curMSProduct); curMSProduct = NULL;
//	}
//	delete(spectrum); spectrum = NULL;
//
//	// this code is only to provide the user some kind of feed back while the program is running
//	//boost::mutex::scoped_lock mutex_locker(mutex1, defer_lock); // defer_lock is initially unlocked
//	//mutex_locker.lock();
//	g_progressCtr++;
//	if(g_progressCtr % 100 == 0 ) cerr << g_progressCtr << " ";
//	if(g_progressCtr % 1000 == 0 ) cerr << endl;
//	//mutex_locker.unlock();
//}





// Function picks out the best score and the next-best score for the data
// contained in the scoreStruct vector passed to it
void PSMClass::pickScores(deque<scoreStruct> &v) {
	deque<scoreStruct>::iterator curScore;
	list<double> score_list;
	list<double>::iterator L_iter;
	scoreStruct *bestScoreStruct = NULL;
	scoreStruct *nextBestScoreStruct = NULL;
	double best_score = 0.0;
	double next_bestScore = 0.0;
	double mz, intensity;
	bool score1_set = false;
	bool score2_set = false;

	MSProductClass *bestMatch = NULL;
	map<double, peakStruct> *Mptr = NULL;
	map<double, double> *spectrum = NULL;
	map<double, vector<double> >::iterator curPeak;


	// store the relevant peak type into a new map that is passed to the curMSProduct object
	spectrum = new map<double, double>;
	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second.at(peakType);
		spectrum->insert( pair<double, double>(mz, intensity) );
	}

	// prep for assignment
	bestScoreStruct = new scoreStruct();
	nextBestScoreStruct = new scoreStruct();

	// put all the scores into a list
	score_list.clear();
	for(curScore = v.begin(); curScore != v.end(); curScore++) {
		score_list.push_back(curScore->score);
	}
	score_list.unique();
	score_list.sort();
	score_list.reverse();


	if((signed)score_list.size() == 1) {
		// This means all permutations had the same score.
		// In this case, we will take the first two scores off the list
		*bestScoreStruct = v.at(0);

		if((signed) v.size() > 1) *nextBestScoreStruct = v.at(1);
		else *nextBestScoreStruct = v.at(0);
	}
	else { // you had at least 2 distinct scores for this PSM
		L_iter = score_list.begin();
		best_score = *L_iter;

		L_iter++;
		next_bestScore = *L_iter;

		for(curScore = v.begin(); curScore != v.end(); curScore++) {
			if( (curScore->score == best_score) && (score1_set == false) ) {
				*bestScoreStruct = *curScore;
				score1_set = true;
			}
			if( (curScore->score == next_bestScore) && (score2_set == false) ) {
				*nextBestScoreStruct = *curScore;
				score2_set = true;
			}

			// this only happens *AFTER* both variables have been assigned a value
			if( score1_set && score2_set ) break;
		}
	}

	score_list.clear(); // cleanup


	// One of the methods above should have produced a bestScore and a 2nd bestScore
	// now record them into the PSMClass variables
	bestScore_final = *bestScoreStruct;
	nextBestScore_final = *nextBestScoreStruct;

	delete(bestScoreStruct); bestScoreStruct = NULL;
	delete(nextBestScoreStruct); nextBestScoreStruct = NULL;


	/*
	 * Best scoring Permutation
	 */
	bestMatch = new MSProductClass( specId, bestScore_final.seq, bestScore_final.charge, nterm_mass );
	bestMatch->updateMZerr(mz_err);
	bestMatch->assignSpectrumMap( *spectrum );

	if(g_useOnlySiteDetermIons) bestMatch->makeSiteDetermIons();

	Mptr = new map<double, peakStruct>();
	bestMatch->getMatchedPeaks(Mptr);

	bestSpectrum.matched_ions = *Mptr;
	bestScore_final.fracMatched = bestMatch->getFractionMatched();

	recordBestSpectrumScores(1);

	delete(bestMatch); bestMatch = NULL;
	delete(Mptr); Mptr = NULL;


	/*
	 * 2nd Best scoring Permutation
	 */
	bestMatch = new MSProductClass( specId, nextBestScore_final.seq, nextBestScore_final.charge, nterm_mass );
	bestMatch->updateMZerr(mz_err);
	bestMatch->assignSpectrumMap( *spectrum );

	if(g_useOnlySiteDetermIons) bestMatch->makeSiteDetermIons();

	Mptr = new map<double, peakStruct>();
	bestMatch->getMatchedPeaks(Mptr);

	nextBestSpectrum.matched_ions = *Mptr;
	nextBestScore_final.fracMatched = bestMatch->getFractionMatched();

	recordBestSpectrumScores(2);

	delete(bestMatch); bestMatch = NULL;
	delete(spectrum); spectrum = NULL;
	delete(Mptr); Mptr = NULL;


	// compute and record the Luciphor delta score
	luciphor_deltaScore = bestScore_final.score - nextBestScore_final.score;
}




// Function generates all possible decoy-phospho peptide permutations
// for the sequence of the current PSM.
void PSMClass::genDecoys() {
	string seq = upperCaseSTY( origModPeptide );

	set<string> *localDecoySet = NULL;
	set<string>::iterator s;

	localDecoySet = new set<string>();

	while( (signed)localDecoySet->size() < (int)numDecoyPermutations ) {
		localDecoySet->insert( genRandDecoyPeptide(seq, numPhosphoSites) );

		//if( (signed)localDecoySet->size() >= decoyLimit ) break;
	}

	for(s = localDecoySet->begin(); s != localDecoySet->end(); s++) {
		decoySet.insert(*s);
	}

	delete(localDecoySet); localDecoySet = NULL;
}




// Function writes final results out
void PSMClass::write_results(ofstream &outf) {

		int isDecoy1 = 0, isDecoy2 = 0;
		double delta_score = 0;
		int N = (signed) raw_spectrum.size();

		delta_score = bestScore_final.score - nextBestScore_final.score;

		// If this PSM represents a peptide where all potential sites are
		// phosphorylated then we set the delta_score to be the best score.
		// We reset 'nextBestScore_final.score' to be zero.
		if(is_unambiguous && (g_randDecoyAA == false)) {
			delta_score = bestScore_final.score;
			nextBestScore_final.score = 0;
		}

		if( isDecoyPep( &bestScore_final.seq ) ) isDecoy1 = 1;
		if( isDecoyPep( &nextBestScore_final.seq ) ) isDecoy2 = 1;

		string tmp;
		if(nterm_mass > 0.0) { // nterminal modification is present
			tmp = "n[" + origModPeptide;
			origModPeptide = tmp;
			tmp.clear();

			tmp = "n[" + bestScore_final.seq;
			bestScore_final.seq = tmp;
			tmp.clear();

			tmp = "n[" + nextBestScore_final.seq;
			nextBestScore_final.seq = tmp;
			tmp.clear();
		}


		if(g_FULL_MONTY) {
			outf << specId << "\t"
				 << peptide << "/+" << charge << "\t"
				 << repModAAchar( &bestScore_final.seq ) << "\t"
				 << repModAAchar( &nextBestScore_final.seq ) << "\t"
				 << iniProb << "\t"
				 << numPhosphoSites << "\t"
				 << numPotentialSites << "\t";

			if(isDecoy1 == 1) outf << "NA\tNA\tNA\t";
			else {
				outf << luciphorProb << "\t"
					 << localFLR << "\t"
					 << globalFLR << "\t";
			}

			outf << delta_score << "\t"
				 << bestScore_final.score << "\t"
				 << nextBestScore_final.score << "\t"
				 << isDecoy1 << "\t"
				 << isDecoy2 << "\t"
				 << N << "\t"
				 << bestScore_final.matchedPeaks << "\t"
				 << nextBestScore_final.matchedPeaks << "\t"
				 << bestScore_final.fracMatched << "\t"
				 << scoreTime;

			outf << endl;
		}
		else { // default output
			outf << specId << "\t"
				 << repModAAchar( &origModPeptide ) << "\t"
				 << repModAAchar( &bestScore_final.seq ) << "\t"
				 << repModAAchar( &nextBestScore_final.seq ) << "\t"
				 << iniProb << "\t"
				 << numSupportingSpectra << "\t"
				 << numPhosphoSites << "\t"
				 << numPotentialSites << "\t";

			if(isDecoy1 == 1) outf << "NA\tNA\tNA\t";
			else {
				outf << luciphorProb << "\t"
					 << localFLR << "\t"
					 << globalFLR << "\t";
			}

			outf << delta_score << "\t"
				 << bestScore_final.score << "\t"
				 << nextBestScore_final.score << "\t"
				 << isDecoy1 << "\t"
				 << isDecoy2 << "\t"
				 << N << "\t"
				 << bestScore_final.matchedPeaks << "\t"
				 << nextBestScore_final.matchedPeaks;

			outf << endl;
		}
}



// Function returns the intensity values stored in 'matchedPeaksMap' or
// 'unmatchedPeaksMap'
list<double> PSMClass::getIntensities(char x, char ionType) {
	list<double> ret;
	map<double, peakStruct>::iterator curPeak;
	map<double, peakStruct> *mapPtr = NULL;
	double intensity = 0.0;

	if(x == 'm') mapPtr = &matchedPeakMap;
	else mapPtr = &unmatchedPeakMap;

	for(curPeak = mapPtr->begin(); curPeak != mapPtr->end(); curPeak++) {
		if(curPeak->second.ionType == ionType) {
			intensity = curPeak->second.intensity;
			ret.push_back( intensity );
		}
	}
	return ret;
}



// Function returns the distance values stored in 'matchedPeaksMap' or
// 'unmatchedPeaksMap'
list<double> PSMClass::getDistances(char x, char ionType) {
	list<double> ret;
	map<double, peakStruct>::iterator curPeak;
	map<double, peakStruct> *mapPtr = NULL;
	double mzDist = 0.0;

	if(x == 'm') mapPtr = &matchedPeakMap;
	else mapPtr = &unmatchedPeakMap;

	for(curPeak = mapPtr->begin(); curPeak != mapPtr->end(); curPeak++) {
		if(curPeak->second.ionType == ionType) {
			mzDist = curPeak->second.MZdistance;
			ret.push_back( mzDist );
		}
	}

	return ret;
}




// Function runs Ascore on the current PSMClass object
void PSMClass::runAscore() {

	map<string,double>::iterator curAscorePermutation;
	set<string>::iterator curPermutation;
	AscoreClass *curAS = NULL;
	string curSeq, pepSeq1, pepSeq2;
	int N = 0;
	double maxScore = 0.0, negLogProb = 0.0;
	int peakDepth = 0;
	int optimalPeakDepth = 0;
	int numIons = 0;
	double pepScore1, pepScore2;
	double delta, mz, intensity, ps;

	vector<double> pep1vec(10,0);
	vector<double> pep2vec(10,0);
	vector<double> scoreDiff(10,0);
	list<double> pepScoreList;
	map<double, int> deltaMap;
	map<double, int>::iterator deltaIter;
	map<string, double> ascoreFinalMap;
	map<double, double> *spectrum = NULL;
	map<double, vector<double> >::iterator curPeak;
	map<double, deque<string> > scoreMap;
	map<double, deque<string> >::iterator scoreMapIter;


	// initialize afs struct of this object for receiving final output
	afs.maxScoreDiff = 0;
	afs.negLogProb = 0;
	afs.nextSeq = "";
	afs.numMatchedPeaks1 = 0;
	afs.numMatchedPeaks2 = 0;
	afs.numPeaksPerBin = 0;
	afs.seqBest = "";



//	if(g_DEBUG_MODE) {
//		if(specId.compare("Pepmix1_HCD.1499.1499.2") != 0) return;
//	}

	// store the relevant peak type into a new map that is passed to the curMSProduct object
	spectrum = new map<double, double>;
	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		mz = curPeak->first;
		intensity = curPeak->second.at(peakType);
		spectrum->insert( pair<double, double>(mz, intensity) );
	}

	/*
	 * Calculate the Peptide Score as described in Gygi's paper.
	 * The PeptideScore tells you which two phospho-permutations are the top
	 * candidates to be subjected to Ascoring.
	 */
	for(curPermutation = phosphoVersionSet.begin(); curPermutation != phosphoVersionSet.end(); curPermutation++) {
		curSeq = *curPermutation;
		curAS = new AscoreClass(curSeq, charge, nterm_mass);
		curAS->assignSpectrumMap( *spectrum );
		curAS->recordMZrange();
		ps = curAS->getPeptideScore();

		scoreMap[ ps ].push_back(curSeq);
		pepScoreList.push_back(ps);

		delete(curAS); curAS = NULL;
	}


	N = (signed) pepScoreList.size(); // get number of scores observed

	if(N == 1) { // this is an unambiguous case and there is no second best permutation
		         // take the value of 'pepScore1' to be the ascore for this PSM

		pepScore1 = pepScoreList.front();

		scoreMapIter = scoreMap.find(pepScore1);
		pepSeq1 = scoreMapIter->second.at(0);

		curAS = new AscoreClass(pepSeq1, charge, nterm_mass);
		curAS->assignSpectrumMap( *spectrum );
		curAS->recordMZrange();

		// identify the max peak depth
		double maxScore = 0;
		int maxScoreIdx = 0;
		for(int i = 1; i <= 10; i++) {
			delta = curAS->getFinalAscore(i);
			if( delta > maxScore ) {
				maxScore = delta;
				maxScoreIdx = i;
			}
		}

		delta = curAS->getFinalAscore( maxScoreIdx );
		afs.negLogProb = delta;
		afs.numPeaksPerBin = maxScoreIdx;
		ascoreMatchedSpectrum1 = curAS->getMatchedSpectrumMap();
		afs.numMatchedPeaks1 = (signed) ascoreMatchedSpectrum1.size();

		afs.maxScoreDiff = afs.negLogProb;
		afs.seqBest = pepSeq1;
		afs.nextSeq = "";
		afs.numMatchedPeaks2 = 0;

		delete(curAS);
	}
	else { // identify the top two permutations that have the highest PeptideScores

		pepScoreList.sort(); // sorted low to high
		pepScoreList.reverse(); // sorted low to high

		pepScore1 = pepScoreList.front();
		pepScoreList.pop_front();
		pepScore2 = pepScoreList.front();

		afs.negLogProb = pepScore1; // the "PeptideScore" as described by Gygi

		/*
		 * Because of the way Ascore works, two different phospho-permutations
		 * can have the exact same score. For this reason we use the 'scoreMap'
		 * object to find all of the candidate permutations with the top 2 scores.
		 */
		pepSeq1.clear();
		pepSeq2.clear();

		scoreMapIter = scoreMap.find(pepScore1);
		int numSeqs  = (signed) scoreMapIter->second.size();

		if(numSeqs > 1) { // more than 1 sequence go the same top score
			pepSeq1 = scoreMapIter->second.at(0);
			pepSeq2 = scoreMapIter->second.at(1);
		}
		else { // at least two distinct scores can be extracted from the map
			pepSeq1 = scoreMapIter->second.at(0);

			scoreMapIter = scoreMap.find(pepScore2);
			pepSeq2 = scoreMapIter->second.at(0);
		}

		/*
		 * At this point you only have to deal with 2 phospho-permutations, the two
		 * with the highest peptide scores.
		 */
		// pepSeq 1
		curAS = new AscoreClass(pepSeq1, charge, nterm_mass);
		curAS->assignSpectrumMap( *spectrum );
		curAS->recordMZrange();
		pep1vec = curAS->getPeptideScoreVec();
		delete(curAS);

		// pepSeq 2
		curAS = new AscoreClass(pepSeq2, charge, nterm_mass);
		curAS->assignSpectrumMap( *spectrum );
		curAS->recordMZrange();
		pep2vec = curAS->getPeptideScoreVec();
		delete(curAS);

		// identify the maximum absolute difference between elements of pep1vec and pep2vec
		for(int i = 0; i < 10; i++) {
			delta = fabs( pep1vec.at(i) - pep2vec.at(i) );
			deltaMap[ delta ] = i;  // k = diff, v = peak depth
			scoreDiff.at(i) = delta;
		}

		sort(scoreDiff.begin(), scoreDiff.end());  //sorted low to high
		delta = scoreDiff.at(9); // last element
		deltaIter = deltaMap.find(delta);
		optimalPeakDepth = deltaIter->second + 1;


		// pepSeq 1 final ascore
		curAS = new AscoreClass(pepSeq1, charge, nterm_mass);
		curAS->assignSpectrumMap( *spectrum );
		curAS->recordMZrange();
		numIons = curAS->getSiteDetermIons();
		pepScore1 = curAS->getFinalAscore(optimalPeakDepth);

		ascoreMatchedSpectrum1 = curAS->getMatchedSpectrumMap();
		afs.numMatchedPeaks1 = (signed) ascoreMatchedSpectrum1.size();
		delete(curAS); curAS = NULL;


		// pepSeq 2 final ascore
		curAS = new AscoreClass(pepSeq2, charge, nterm_mass);
		curAS->assignSpectrumMap( *spectrum );
		curAS->recordMZrange();
		curAS->getSiteDetermIons();
		pepScore2 = curAS->getFinalAscore(optimalPeakDepth);
		ascoreMatchedSpectrum2 = curAS->getMatchedSpectrumMap();
		afs.numMatchedPeaks2 = (signed) ascoreMatchedSpectrum2.size();
		delete(curAS); curAS = NULL;

		// Get final Ascore values
		afs.numPeaksPerBin = optimalPeakDepth;
		afs.seqBest = pepSeq1;
		afs.nextSeq = pepSeq2;
		afs.maxScoreDiff = fabs(pepScore1 - pepScore2);


		if(g_DEBUG_MODE) {
			cerr << "\n## PSMClass::runAscore():\n"
				 << "PeptideScore:   " << afs.negLogProb << endl
				 << "Peptide 1:      " << afs.seqBest << endl
				 << "Peptide 2:      " << afs.nextSeq << endl
				 << "Score 1:        " << pepScore1 << endl
				 << "Score 2:        " << pepScore2 << endl
				 << "Peak Depth:     " << afs.numPeaksPerBin << endl
				 << "Matched Pks 1:  " << afs.numMatchedPeaks1 << endl
				 << "Mathced Pks 2:  " << afs.numMatchedPeaks2 << endl
				 << "Ascore:         " << afs.maxScoreDiff << endl;
		}

	} //end else

	g_progressCtr++;
	if(g_progressCtr % 100 == 0 ) cerr << g_progressCtr << " ";
	if(g_progressCtr % 1000 == 0 ) cerr << endl;
}




// Function writes Ascore results to disk
void PSMClass::write_ascore_results(ofstream &outf) {
	string tmp;

	// this line handles unambiguous cases.
	if(afs.nextSeq.length() == 0) afs.nextSeq = afs.seqBest;

	if(nterm_mass > 0.0) { // nterminal modification is pressent
		tmp = "n[" + origModPeptide;
		origModPeptide = tmp;
		tmp.clear();

		tmp = "n[" + afs.seqBest;
		afs.seqBest = tmp;
		tmp.clear();

		if(afs.nextSeq.length() > 0) {
			tmp = "n[" + afs.nextSeq;
			afs.nextSeq = tmp;
			tmp.clear();
		}
	}

//	if( isDecoyPep( &afs.seqBest ) ) isDecoy1 = 1;
//	if( isDecoyPep( &afs.nextSeq ) ) isDecoy2 = 1;

	outf << specId << "\t"
		 << repModAAchar( &origModPeptide ) << "/+" << charge << "\t"
		 << repModAAchar( &afs.seqBest ) << "\t"
		 << repModAAchar( &afs.nextSeq ) << "\t"
		 << iniProb << "\t"
		 << numPhosphoSites << "\t"
		 << numPotentialSites << "\t"
		 << afs.numPeaksPerBin << "\t"
		 << afs.negLogProb << "\t"
		 << afs.maxScoreDiff << "\t"
//		 << isDecoy1 << "\t"
//		 << isDecoy2 << "\t"
		 << (signed) raw_spectrum.size() << "\t"
		 << afs.numMatchedPeaks1 << "\t"
		 << afs.numMatchedPeaks2 << "\n";

}


// Function records the scores for all matched peaks in a spectrum and stores them
// into the 'bestSpectrum' struct
void PSMClass::recordBestSpectrumScores(int J) {

	map<double, peakStruct>::iterator curPeak;
	matchedSpectrumStruct *ptr = NULL;

	// determine which permutation will be printed: best or 2nd best
	if(J == 1) ptr = &bestSpectrum;
	else ptr = &nextBestSpectrum;

	double mz = 0.0, intensity = 0.0, mzDist = 0.0;
	double Iscore = 0.0, Dscore = 0.0, x = 0.0;
	string ionSeq;
	char ionType = 'X';
	size_t found;

	double muM = 0.0, varM = 0.0; // matched params
	double muM_d = 0.0, varM_d = 0.0;
	double log_prob_M = 0.0, log_dist_M = 0.0;

	double muU = 0.0, varU = 0.0; // unmatched params
	double muU_d = 0.0, varU_d = 0.0;
	double log_prob_U = 0.0, log_dist_U = 0.0;

	// unmatched parameters are the same regardless of the ion type
	muU = local_params.unMatched_mean;
	varU = local_params.unMatched_var;
	muU_d = local_params.unMatched_dist_mean;
	varU_d = local_params.unMatched_dist_var;

	//for(curPeak = bestSpectrum.matched_ions.begin(); curPeak != bestSpectrum.matched_ions.end(); curPeak++) {
	for(curPeak = ptr->matched_ions.begin(); curPeak != ptr->matched_ions.end(); curPeak++) {


		mz = curPeak->first;
		ionSeq = curPeak->second.ionStr;
		intensity = curPeak->second.intensity;
		mzDist =curPeak->second.MZdistance;
		ionType = curPeak->second.ionType;


		if(ionType == 'b') {
			muM = local_params.matched_mean_b;
			varM = local_params.matched_var_b;
			muM_d = local_params.matched_dist_mean_b;
			varM_d = local_params.matched_dist_var_b;
		}
		else { // y-ion
			muM = local_params.matched_mean_y;
			varM = local_params.matched_var_y;
			muM_d = local_params.matched_dist_mean_y;
			varM_d = local_params.matched_dist_var_y;
		}


		if(g_IS_HCD) {
			/*
			 * INTENSITY
			 */
			log_prob_M = log_gaussianProb(muM, varM, intensity);
			log_prob_U = log_gaussianProb(muU, varU, intensity);
			Iscore = log_prob_M - log_prob_U;
			ptr->IscoreMap[ ionSeq ] = Iscore;

			/*
			 * DISTANCE
			 */
			log_dist_M = log_laplaceProb(muM_d, varM_d, mzDist);
			log_dist_U = log_uniformProb(-mz_err, mz_err);
			Dscore = log_dist_M - log_dist_U;
			ptr->DscoreMap[ ionSeq ] = Dscore;

			double intense_wt = 1.0 / ( 1.0 + exp(-Iscore) );
			x = intense_wt * Dscore;
		}
		else { // for CID data
			/*
			 * INTENSITY
			 */
			log_prob_M = log_gaussianProb(muM, varM, intensity);
			log_prob_U = log_gaussianProb(muU, varU, intensity);
			Iscore = log_prob_M - log_prob_U;
			ptr->IscoreMap[ ionSeq ] = Iscore;

			/*
			 * DISTANCE
			 */
			log_dist_M = log_gaussianProb(muM_d, varM_d, mzDist);
			log_dist_U = log_gaussianProb(muU_d, varU_d, mzDist);
			Dscore = log_dist_M - log_dist_U;
			ptr->DscoreMap[ ionSeq ] = Dscore;

			double intense_wt = 1.0 / ( 1.0 + exp(-Iscore) );
			x = intense_wt * Dscore;
		}

		ptr->FinalScoreMap[ ionSeq ] = x;

	}
}






// Function to be passed to a threadpool and executed for each charge state
void PSMClass::threaded_recordModelingParameters_matched() {
	list<double>::iterator L;
	list<double> X;

	setSpectrumPtr("median");
	generatePermutations();
	classifyPeaks();

/**************
	if( use_for_model ) {

		X.clear();
		X = getIntensities('m', 'y');
		X.unique();
		for(L = X.begin(); L != X.end(); L++) M_ints_y.push_back(*L);

		X.clear();
		X = getIntensities('m', 'b');
		X.unique();
		for(L = X.begin(); L != X.end(); L++) M_ints_b.push_back(*L);

		X.clear();
		X = getDistances('m', 'y');
		X.unique();
		for(L = X.begin(); L != X.end(); L++) M_dist_y.push_back(*L);

		X.clear();
		X = getDistances('m', 'b');
		X.unique();
		for(L = X.begin(); L != X.end(); L++) M_dist_b.push_back(*L);

	}
*****************/
}




/*********************************
// Function to be passed to a threadpool and executed for each charge state
void PSMClass::threaded_recordModelingParameters_UNmatched() {
	list<double> X;
	list<double>::iterator L;

	setSpectrumPtr("median");

	if( use_for_model ) {
		X.clear();
		X = getIntensities('u', 'u');
		X.unique();
		for(L = X.begin(); L != X.end(); L++) U_ints_local.push_back(*L);


		X.clear();
		X = getDistances('u', 'u');
		X.unique();
		for(L = X.begin(); L != X.end(); L++) U_dist_local.push_back(*L);
	}

}
*******************************/




// Function returns pointer to the requested parameter list
// The pointer is constant so the data cannot be written to, only read from
list<double>* PSMClass::getParamList(char matchType, char ionType, char dataType) {
	list<double> *ret = new list<double>;
	map<double, peakStruct>::iterator curPeak;
	map<double, peakStruct> *mapPtr = NULL;
	double intensity, dist;
	int score;

	if(matchType == 'm') mapPtr = &matchedPeakMap;
	else mapPtr = &unmatchedPeakMap;

	for(curPeak = mapPtr->begin(); curPeak != mapPtr->end(); curPeak++) {
		intensity = curPeak->second.intensity;
		dist = curPeak->second.MZdistance;

		score = 0;
		if( !isInfinite(intensity) && !dbl_isnan(intensity) ) score++;

		if( (curPeak->second.ionType == ionType) && (score == 1) ) {
			if(dataType == 'i') ret->push_back(intensity);
			if(dataType == 'd') ret->push_back(dist);
		}
	}

	ret->unique();

	return ret;
}



// Function returns pointer to the requested parameter list
// The pointer is constant so the data cannot be written to, only read from
//list<double>* PSMClass::getParamList(char matchType, char ionType, char dataType) {
//	list<double> *retPtr = NULL;
//
//	if(dataType == 'i') { // intensity
//		if(matchType == 'm') {
//			if(ionType == 'b') retPtr = &M_ints_b;
//			else if(ionType == 'y') retPtr = &M_ints_y;
//		}
//		else retPtr = &U_ints_local;
//	}
//	else if(dataType == 'd') { // distance
//		if(matchType == 'm') {
//			if(ionType == 'b') retPtr = &M_dist_b;
//			else if(ionType == 'y') retPtr = &M_dist_y;
//		}
//		else retPtr = &U_dist_local;
//	}
//
//	return retPtr;
//}




// Function to be passed to a threadpool and executed for each charge state
void PSMClass::threaded_scorePSM() {
	setSpectrumPtr("median");
	generatePermutations();
	if(g_randDecoyAA)  genDecoys();
	calcScore();
}





// Function computes the number of decoy permutations that can be created for
// this PSM
void PSMClass::calcNumDecoyPermutations() {
	string seq = upperCaseSTY( origModPeptide );
	int seqLen = (signed) seq.length();
	int N = 0;
	char c;


	// find out how many candidate residues you have for making decoys
	for(int i = 0; i < seqLen; i++) {
		c = seq.at(i);
		if( (c != 'S') && (c != 'T') && (c != 'Y') && (c != 'X') && ( !islower(c) ) ) N++;
		//if( (c != 'S') && (c != 'T') && (c != 'Y') && (c != 'X') ) N++;
	}

	if( N >= numPotentialSites ) { // you have enough non-STY residues to make a decoy peptide
		numDecoyPermutations = combinatorial( (double)N, (double)numPhosphoSites );
	}
	else numDecoyPermutations = 0;


}



// Function returns the number of distinct sequence permutations for the
// sequence associated with this PSM
double PSMClass::getNumPerms() {
	double ret = 0;
	ret = numPermutations + numDecoyPermutations;
	return ret;
}



// Function records all of the necessary data for computing the FLR
// for this PSM into the passed struct
flrStruct PSMClass::getFLRdata() {
	double deltaScore = bestScore_final.score - nextBestScore_final.score;
	flrStruct ret;
	ret.specId = specId;
	ret.deltaScore = deltaScore;
	ret.isDecoy = isDecoyPep( &bestScore_final.seq );
	ret.globalFLR = -1000;
	ret.localFLR = -1000;
	ret.prob = -1000;
	return ret;
}


// Assign FLR data to this PSMClass object from the given "filled-in" FLR object
void PSMClass::setFLR( flrStruct *ptr ) {
	globalFLR = ptr->globalFLR;
	localFLR = ptr->localFLR;
	luciphorProb = ptr->prob;
}



// Function removes all collected results from the PSMclass object
void PSMClass::clear() {
	phosphoVersionSet.clear();
	matchedPeakMap.clear();
	unmatchedPeakMap.clear();
	M_ints_y.clear();
	M_ints_b.clear();
	M_dist_y.clear();
	M_dist_b.clear();
	U_ints_local.clear();
	U_dist_local.clear();
	mz_err = 0;

}



// Function generates a random peptide string assigned to this PSM
void PSMClass::randomizeSeq() {

	int N = (signed) peptide.length();
	int n = 0;
	char X;
	string residues = "ACDEFGHIKLMNPQRSTVWY"; // normal amino acid residues
	string newPep = "";

	deque<char> aa;
	// add all the normal amino acid residues to 'aa'
	for(int i = 0; i < 20; i++) aa.push_back( residues.at(i) );

	n = (signed) aa.size();


	// figure out how many modified residues (ie: lower case characters)
	// the original peptide for this PSM had and limit the random peptide
	// to this number of modified residues
	int modCount = 0;
	for(int i = 0; i <N; i++) {
		X = modPeptide.at(i);
		if(islower(X)) modCount++;
	}


	for(int i = 0; i < N; i++) {

		int p = rand() % n; // get index of a random element from 'aa' deque
		X = aa.at(p);

		// adding STY+80 residues 50% of the time
		if( (X == 'S') || (X == 'T') || (X == 'Y') ) {
			double r = ((double) rand() / RAND_MAX); // random value btwn 0 - 1
			if(r > 0.5) {
				char x = tolower(X);
				X = x;
			}
		}

		if( islower(X) ) {
			if(modCount > 0) { // we can append this character as is
				newPep += X;
				modCount--;
			}
			else { // we have maxed out the number of modified residues for this peptide
				newPep += toupper(X);
			}
		}
		else newPep += X; // uppercase residue
	}

	modPeptide = newPep;
	peptide = uc(newPep);
	is_randomized_seq = true;
}




// Perform deisotoping of spectrum
void PSMClass::deisotopeSpectrum() {
	map<double, double> peakClass; // k = an observed peak, v = it's monoisotopic parent (if any)
	map<double, vector<double> >::iterator curPeak;
	deque<double> mzValues;
	double curMZ, obs_mz, delta, theo_mz;
	const double DA_ERR = 0.01; // about 10ppm
	int N = 0;

	for(curPeak = raw_spectrum.begin(); curPeak != raw_spectrum.end(); curPeak++) {
		curMZ = curPeak->first;
		mzValues.push_back(curMZ);
		peakClass[ curMZ ] = 0; // a zero value here means the peak is not assigned to a monoisotopic parent
	}

	N = (signed) mzValues.size();
	sort( mzValues.begin(), mzValues.end() ); // sort m/z values from low to high

	for(double z = 1; z < 3; z++) { // assume isotopic peaks up to charge states of 2;

		for(int i = 0; i < N; i++) {
			curMZ = mzValues.at(i);

			for(double j = 1; j < 4; j++) { // up to 3 isotopic peak
				theo_mz = curMZ + ((1.0/z) * j); // current theoretical isotopic peak for 'curMZ'

				for(int k = (i+1); k < N; k++) {
					obs_mz = mzValues.at(k);
					delta = fabs(obs_mz - theo_mz);

					if(delta < DA_ERR) peakClass[ obs_mz ] = curMZ;
					else if(obs_mz > theo_mz) break;
				}
			}
		}
	}

	/*
	 * At this point, all peaks have been classified. Peaks where the map value
	 * is zero in peakClass will be retained and all others removed from spectrum.
	 */
	for(int i = 0; i < N; i++) {
		curMZ = mzValues.at(i);
		int pClass = peakClass[ curMZ ];
		if(pClass != 0) {
			curPeak = raw_spectrum.find(curMZ);
			raw_spectrum.erase(curPeak);
		}
	}
}

