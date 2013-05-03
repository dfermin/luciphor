/*
 * PepXMLClass.cpp
 *
 *  Created on: Mar 18, 2011
 *      Author: dfermin
 */

#include <iostream>
#include <fstream>
#include <utility>
#include <iomanip> // for pretty formatting
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/multi_array.hpp>
#include <boost/regex.hpp>
#include <boost/thread.hpp>  // for multi-threading
#include <boost/bind.hpp>    // for multi-threading
#include "threadpool.hpp"  // for multi-threading
#include <algorithm>
#include <map>
#include <deque>
//#include "PepClass.hpp"
#include "FLRClass.hpp"
#include "PepXMLClass.hpp"
#include "PSMClass.hpp"
#include "statsFunctions.hpp"
#include "structs.hpp"
#include "nonParamDist.hpp"



// for reading spectra files using pwiz library
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/data/msdata/MSData.hpp"

using namespace std;
using namespace pwiz::msdata;
using namespace pwiz::cv;
using namespace boost;


int g_progressCtr = 0;
int g_totalNumPSM = 0;
double g_MIN_DIST = -20.00;


const double min_deltaCN = 0.1;


// Function opens the user-given pepXML file and parses out the relevant data
void PepXMLClass::parsePepXMLfile() {

	matchDataStruct *mds = NULL; // holds information about the current spectrum_query record
	string line, AA, tmp;
	char isVar = 'X';
	boost::smatch matches; // used to capture REGEX matches, element[0] is the whole string
	double mass;
	int modPos;
	int score = 0;  // must be >= 4 or the PSM doesn't get recorded, some PSMs don't have any modifications
	PSMClass *curPSM = NULL;

	numPSMs_forModeling = 0;

	PSMvec = new deque<PSMClass>;

	fileName = g_srcXMLfile;

	/*
	 * These are the regex pattern objects to parse the pepXML file
	 */

	boost::regex aminoacid_mod_regex("^.*<aminoacid_modification aminoacid=\"([^\"])\".*mass=\"([^\"]+)\".*variable=\"(Y|N)\".*");
	boost::regex spectrum_query_regex("^.*<spectrum_query spectrum=\"([^\"]+)\".+start_scan=\"([^\"]+)\".+ assumed_charge=\"(\\d)\".*");
	boost::regex search_hit_regex("^.*<search_hit hit_rank=\"1\" peptide=\"(\\w+)\" peptide_prev_aa.* calc_neutral_pep_mass=\"([^\"]+)\".*");
	boost::regex modification_mass_regex("^.*<mod_aminoacid_mass position=\"(\\d+)\" mass=\"([^\"]+)\".*");
	boost::regex nterm_mod_regex("^.*<modification_info mod_nterm_mass=\"([^\"]+)\".*");
	boost::regex end_spectrum_query("^.*</search_hit>.*");

	// depending upon what the user chooses as a PSM selection criterion, the scoring_regex_ptr
	// is set to point to the right regex pattern
	boost::regex *scoring_regex_ptr = NULL;
	boost::regex peptideprophet_result_regex("^.*<peptideprophet_result probability=\"([^\"]+)\" all_ntt.*");
	boost::regex sequest_xcorr_regex("^.*<search_score name=\"xcorr\" value=\"([^\"]+)\".*");
	boost::regex sequest_deltacn_regex("^.*<search_score name=\"deltacn\" value=\"([^\"]+)\".*");
	boost::regex xtandem_expect_regex("^.*<search_score name=\"expect\" value=\"([^\"]+)\".*");
	boost::regex mascot_ionscore_regex("^.*<search_score name=\"ionscore\" value=\"([^\"]+)\".*");

	scoring_regex_ptr = &peptideprophet_result_regex; // the default
	if(g_scoringMethod == 1) scoring_regex_ptr = &sequest_xcorr_regex;
	if(g_scoringMethod == 2) scoring_regex_ptr = &xtandem_expect_regex;
	if(g_scoringMethod == 3) scoring_regex_ptr = &mascot_ionscore_regex;


	ifstream in;
	in.open(g_srcXMLfile.c_str());

	if(!in) {
		cerr << "\nUnable to open '" << g_srcXMLfile << "'\nExiting...\n\n";
		exit(0);
	}

	cerr << "\nParsing: " << g_srcXMLfile << endl;

	while( !in.eof() ) {

		getline(in, line);

		if( boost::regex_match(line, matches, aminoacid_mod_regex) ) {
			// record this AA modification into the AAmass map object
			AA.assign(matches[1].first, matches[1].second);

			tmp.assign(matches[2].first, matches[2].second);
			double modMass = strtod(tmp.c_str(), NULL);
			tmp.clear();

			tmp.assign(matches[3].first, matches[3].second);
			isVar = tmp[0];

			addAAmass(AA, modMass, isVar);
			tmp.clear();
		}


		//spectrum_query line
		else if( boost::regex_match(line, matches, spectrum_query_regex)) {
			score = 0; // prep for next iteration

			mds = new matchDataStruct();
			mds->specId.assign(matches[1].first, matches[1].second); // specId

			tmp = "";
			tmp.assign(matches[2].first, matches[2].second); // scan number
			mds->scanNum = atoi(tmp.c_str());
			tmp.clear();

			tmp.assign(matches[3].first, matches[3].second); // charge state
			mds->charge = atoi(tmp.c_str());
			tmp.clear();
			score++;
		}


		//search_hit line
		else if( boost::regex_match(line, matches, search_hit_regex) ) {
			mds->peptide.assign(matches[1].first, matches[1].second);

			tmp.assign(matches[2].first, matches[2].second);
			mass = strtod(tmp.c_str(), NULL);
			mds->mass = mass;
			tmp.clear();
			score++;
		}

		//mod_aminoacid_mass line
		else if( boost::regex_match(line, matches, modification_mass_regex)) {
			tmp.assign(matches[1].first, matches[1].second);
			modPos = atoi(tmp.c_str()) - 1; // zero-bases indexing
			tmp.clear();

			tmp.assign(matches[2].first, matches[2].second);
			mass = strtod(tmp.c_str(), NULL);
			tmp.clear();

			mds->mods[ modPos ] = mass; // record modification position
			score++;
		}

		else if( boost::regex_match(line, matches, nterm_mod_regex)) {
			tmp.assign(matches[1].first, matches[1].second);
			mass = strtod(tmp.c_str(), NULL);
			tmp.clear();

			mds->mods[ -1 ] = mass; // -1 means n-terminus
		}


		// whatever the PSM selection criterion is, we will store it in
		// the iniProb variable of the PSMClass
		else if( boost::regex_match(line, matches, *scoring_regex_ptr) ) {
			tmp.assign(matches[1].first, matches[2].second);

			double d = strtod(tmp.c_str(), NULL);
			if(g_scoringMethod == 2) {
				if(d < TINY_NUM) d = TINY_NUM; // prevents errors when taking long
				mds->iniProb = -1.0 * log(d);
			}
			else mds->iniProb = d;
			tmp.clear();
			score++;
		}


		//end of record, decide if you are going to keep it
		else if( boost::regex_match(line, end_spectrum_query) ) {

			if(score >= 4) {
				// if we are limiting charge states, this line of code will exclude all PSMs
				// that do not match the requested charge state.
				if(g_LIMIT_CHARGE_STATE && (mds->charge != g_CHARGE_STATE)) mds->iniProb = -10000;

				// double check that the amino acid sequence of the peptide is valid
				// some databases have invalid peptide characters
				if( allValidAAs(&mds->peptide) == false ) mds->iniProb = -10000;

				if(mds->iniProb >= g_prob_threshold) { // meets our score threshold

					curPSM = new PSMClass(mds);

					if(curPSM->isKeeper()) PSMvec->push_back(*curPSM);

					// we'll use this spectrum for estimating the model parameters
					if(curPSM->useForModeling()) numPSMs_forModeling++;

//					PepClass *pc = new PepClass(mds);
//					delete(pc); pc = NULL;
				}
			}

			delete(curPSM);
			curPSM = NULL;

			delete(mds);
			mds = NULL;
		}
	}
	in.close();

	if(curPSM != NULL) {
		PSMvec->push_back(*curPSM);
		delete(curPSM); curPSM = NULL;
	}

	numPSM = (signed) PSMvec->size();

	if(numPSM == 0) {
		cerr << "\nERROR!: I was unable to parse out any PSMs!!\nExiting now...\n";
		exit(-1);
	}

	if(numPSMs_forModeling < g_MIN_MODEL_NUM) {
		cerr << "\nERROR!: I was unable to find at least " << g_MIN_MODEL_NUM
			 << " PSMs with a score >= " << g_model_prob << endl
			 << "If you want to analyze the data in '" 
			 << g_srcXMLfile.substr(g_srcXMLfile.find_last_of("/")+1) << "' " 
			 << "you must lower your modeling threshold.\nExiting now...\n\n";
		exit(-1);
	}


	cerr << numPSM << " phospho-peptides read in... (Score >= " << g_prob_threshold << ")\n\n";
}



// Function reads in the spectra from the data stored for each PSMClass object
void PepXMLClass::readInSpectra() {
	map<string, SpecStruct> spectraLib;
	string curSpectrumFilePath;
	map<string, deque<scanIdStruct> >::iterator fileIter;
	deque<scanIdStruct>::iterator curScan;
	deque<PSMClass>::iterator curPSM;

	scanIdStruct *curScanId = NULL;
	SpecStruct *curSpectrum = NULL;
	MSDataFile *msd = NULL;

	// Extract from each PSM, it's parent spectrum file and scan number
	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
		curScanId = new scanIdStruct();

		curScanId->spectrumFile = curPSM->getSpectrumFileName();
		curScanId->scanNumber = curPSM->getScanNumber();
		curScanId->specId = curPSM->getSpecId();

		scanMap[ curScanId->spectrumFile ].push_back(*curScanId);
		delete(curScanId);
		curScanId = NULL;
	}

	if(g_ext == "mgf") {
		parseMGF();
		return;
	}


	// iterate over each file name extracting the relevant spectra from it
	for(fileIter = scanMap.begin(); fileIter != scanMap.end(); fileIter++) {
		cerr << fileIter->first << ": ";

		filesystem::path curFile( fileIter->first.c_str() );
		filesystem::path spectral_dir( g_srcDir.c_str() );
		filesystem::path curFilePath( spectral_dir/curFile );
		curSpectrumFilePath = curFilePath.file_string();

		// msd is the object that holds the current spectrum file in memory
		msd = new MSDataFile( curSpectrumFilePath.c_str() );
		if(!msd->run.spectrumListPtr.get()) {
			cerr << "\nNo spectra found in '" << fileIter->first << "'\n\n";
			exit(-1);
		}

		// get pointers to all the spectra in the current spectrum file
		SpectrumListPtr allSpectra = msd->run.spectrumListPtr;

		// turn our scan numbers into navtieId strings for pwiz
		CVID nativeIdFormat_ = id::getDefaultNativeIDFormat(*msd);

		// now get all the spectra of relevance for this file.
		for(curScan = fileIter->second.begin(); curScan != fileIter->second.end(); curScan++) {

			curSpectrum = new SpecStruct;
			*curSpectrum = getSpectrum(*curScan, allSpectra, nativeIdFormat_);

			spectraLib[ curScan->specId ] = *curSpectrum;
			delete(curSpectrum);
			curSpectrum = NULL;
		}

		delete(msd); msd = NULL; // this should 'close' the spectrum file object
		curSpectrumFilePath.clear();
		cerr << fileIter->second.size() << " spectra read in.\n";
		logF << fileIter->first << ": " << fileIter->second.size() << " spectra read in.\n";
	}

	// assign the data in spectraLib to each of the PSMClass objects in PSMvec
	map<string, SpecStruct>::iterator curSS;
	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
		curSS = spectraLib.find( curPSM->getSpecId() );
		curPSM->recordSpectrum( curSS->second );
	}
	spectraLib.clear(); // clean up

	cerr << endl; // prettier stderr
}




/*
 * Function to read in MGF files. The Proteowizard library is having issues with these
 * so we have to use our own code for it.
 */
void PepXMLClass::parseMGF() {

	map<string, deque<scanIdStruct> >::iterator m;
	deque<scanIdStruct>::iterator curScan;
	map<string, SpecStruct> *spectrumMapPtr = NULL;
	map<string, SpecStruct>::iterator ssIter;
	deque<PSMClass>::iterator curPSM;
	vector<string> *vPtr = NULL;
	SpecStruct *curSpectrum = NULL;
	string curSpectrumFilePath;
	string specId, line;
	ifstream mgf;
	double mz, intensity;
	int ctr;

	for(m = scanMap.begin(); m != scanMap.end(); m++) {

		spectrumMapPtr = new map<string, SpecStruct>;
		ctr = 0;
		cerr << m->first << ": ";

		filesystem::path curFile( m->first.c_str() );
		filesystem::path spectral_dir( g_srcDir.c_str() );
		filesystem::path curFilePath( spectral_dir/curFile );
		curSpectrumFilePath = curFilePath.file_string();

		mgf.open(curSpectrumFilePath.c_str(), ios::in);

		if( !mgf.is_open() ) {
			cerr << "\nERROR: Unable to open '" << curSpectrumFilePath << endl;
			exit(0);
		}


		// iterate line-by-line over current mgf file
		while( !mgf.eof() ) {
			line.clear();
			getline(mgf, line);

			if(line.substr(0,7) == "END IONS") { // end of a spectrum record
				if( !specId.empty() )
					spectrumMapPtr->insert(pair<string, SpecStruct>(specId, *curSpectrum));

				delete(curSpectrum);
				curSpectrum = NULL;
			}

			else if(line.substr(0,10) == "BEGIN IONS") { // beginning of a new spectrum
				curSpectrum = new SpecStruct;
				mz = 0;
				intensity = 0;
				specId.clear();
			}

			else if(line.substr(0,6) == "TITLE=") { // specId line
				specId = line.substr(6);
			}

			else if( isdigit(line[0]) ) { // you've hit a peak entry
				vPtr = new vector<string>(2);
				*vPtr = split_string(line);

				mz = str2dbl( vPtr->at(0) );
				intensity = str2dbl( vPtr->at(1) );

				curSpectrum->mz.push_back(mz);
				curSpectrum->intensity.push_back(intensity);
				delete(vPtr);
			}
		}
		// close current MGF file
		mgf.close();

		// record the last spectrum
		if(curSpectrum != NULL) {
			spectrumMapPtr->insert(pair<string, SpecStruct>(specId, *curSpectrum));
			delete(curSpectrum);
			curSpectrum = NULL;
		}


		// Iterate over the specId values for the current file. 'm' is an iterator
		// that points to the scanMap object parsed from the interact.pep.xml file
		for(curScan = m->second.begin(); curScan != m->second.end(); curScan++) {
			ssIter = spectrumMapPtr->find(curScan->specId);

			// found match
			if(ssIter != spectrumMapPtr->end()) {
				for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
					if(curPSM->getSpecId() == curScan->specId) {
						curPSM->recordSpectrum( ssIter->second );
						ctr++;
					}
				}
			}
		}

		// clean up memory for next iteration
		delete(spectrumMapPtr);
		spectrumMapPtr = NULL;
		cerr << ctr << " spectra read in.\n";
		logF << m->first << ": " << ctr << " spectra read in.\n";
	}
}


/*
 * Function uses the proteo-wizard library to read in the spectrum from the
 * information provided in 'scanInfo'
 */
SpecStruct PepXMLClass::getSpectrum(scanIdStruct curScan, SpectrumListPtr allSpectra, CVID &nativeIdFormat_) {
	SpecStruct ret;
	vector<double> cur_mz;
	vector<double> cur_intensity;
	string scanId = int2string(curScan.scanNumber);
	const bool getBinaryData = true; // you need this variable to get the actual spectrum
	string nativeScanId = id::translateScanNumberToNativeID(nativeIdFormat_, scanId);

	// find the index of the spectrum with this scan number
	size_t foundIndex = allSpectra->find(nativeScanId);

	if(foundIndex != allSpectra->size()) { // found spectra

		SpectrumPtr spectrumPtr = allSpectra->spectrum(foundIndex, getBinaryData);

		// actually get the current spectrum
		auto_ptr<SpectrumInfo> specInfo(new SpectrumInfo());
		specInfo->SpectrumInfo::update(*spectrumPtr, true);

		// fill in MZIntensityPair vector for access to the binary data
		vector<MZIntensityPair> pairs;
		vector<MZIntensityPair>::const_iterator v2_iter;
		spectrumPtr->getMZIntensityPairs(pairs);

		for(v2_iter = pairs.begin(); v2_iter != pairs.end(); v2_iter++) {
			ret.mz.push_back( v2_iter->mz );
			ret.intensity.push_back( v2_iter->intensity );
		}
	}

	return ret;
}




// This function removes PSMs that cannot be scored because they don't have
// enough high-scoring examples for modeling.
void PepXMLClass::prunePSMdeq() {

	deque<int> chargeDeq;
	deque<PSMClass>::iterator curPSM;
	map<int, prunePSMstruct> zMap;
	map<int, prunePSMstruct>::iterator curZ;
	int z;
	string sm = ""; // *S*coring *M*ethod

	cerr << "\nLooking for low quality PSM scores...\n";

	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
		z = curPSM->getCharge();

		curZ = zMap.find(z);
		if(curZ == zMap.end()) {
			zMap[z].totalZ = 0;
			zMap[z].modelZ = 0;
			zMap[z].maxScore = TINY_NUM;
		}

		zMap[z].totalZ++;
		if(curPSM->useForModeling()) zMap[z].modelZ++;
		if( curPSM->getIniProb() > zMap[z].maxScore ) zMap[z].maxScore = curPSM->getIniProb();
	}

	// this just helps report the more clearly to the user the cause of the problem
	if(g_scoringMethod == 0) sm = "Probability";
	if(g_scoringMethod == 1) sm = "Sequest XCorr";
	if(g_scoringMethod == 2) sm = "-log(E-value)";
	if(g_scoringMethod == 3) sm = "Mascot Ion Score";

	for(curZ = zMap.begin(); curZ != zMap.end(); curZ++) {
		// all PSMs for this charge state must be removed
		if(curZ->second.maxScore < g_model_prob) {
			cerr << "Maximum " << sm << " score observed for +" << curZ->first << " "
				 << "charge state PSMs: " << curZ->second.maxScore << endl
				 << "This is below the modeling threshold of: " << g_model_prob << "\n"
			     << "Discarding " << curZ->second.totalZ << " PSM's "
				 << "due to insufficient data for modeling of their charge state.\n"
				 << "If you want to try and 'capture' these cases, adjust "
				 << "your Luciphor parameters\n";

			z = curZ->first;
			for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
				if(curPSM->getCharge() == z) curPSM = PSMvec->erase(curPSM);
			}
		}
	}
}



// This function removes any PSMs having # permutations > g_NUM_PERMS_LIMIT
void PepXMLClass::removeComplexPSMs() {
	deque<PSMClass>::iterator curPSM;
	double numPerms = 0;
	int droppedPSMs = 0;

	// this for-loop syntax is necessary because we are delete from the PSM m
	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); ) {
		numPerms = curPSM->getNumPerms();
		if(numPerms > g_NUM_PERMS_LIMIT) {
			curPSM = PSMvec->erase(curPSM);
			numPSM--;
			droppedPSMs++;
		}
		else curPSM++;
	}
	if(droppedPSMs > 0) {
		cerr << endl << droppedPSMs << " PSMs were removed because of their "
			 << "complexity.\nIf you would like to score these cases, rerun "
			 << "luciphor with the '-k' option increasing the number beyond "
			 << g_NUM_PERMS_LIMIT << endl << endl;
	}
}




// This function gets the modeling parameters needed to score the spectra
void PepXMLClass::acquireModelParameters() {

	deque<PSMClass>::iterator curPSM;
	map<int, modelParamStruct> modelParametersMap;
	map<int, modelParamStruct>::iterator m;
	map<int, double>::iterator cur_mz_err;
	map<int, int> chargeFreqMap;
	map<int, int>::iterator cfm_iter;

	list<double> M_ints_Y, M_ints_B, M_dist_Y, M_dist_B;
	list<double> U_ints, U_dist;
	list<double> *ptr  = NULL;
	list<double>::iterator L;

	set<int> missedChargeStates;
	set<int>::iterator s_iter;

	modelParamStruct *curParams = NULL;


	double meanM_b, meanM_y, varM_b, varM_y, meanU, varU;
	double meanMd_b, meanMd_y, varMd_b, varMd_y, meanUd, varUd;
	double mz_err_tol;
	double x;
	int i;
	int curChargeState = 0;
	int maxChargeState = 0;
	int zN; // used to hold count of number of PSM for a given charge state
	string msg;
	fstream debug_modelData;

	// A pool of reusable threads that will process the data.
	// The default number of threads is 1
	boost::threadpool::pool TP( g_NUM_THREADS );

	// First step is to identify all of the matched peaks from all of the spectra
	// and collect the intensities & distances of the classified peaks.
	// We will only use spectra with a probability >= modelProb
	// to compute the model parameters.
	i = 0;
	cerr << "\nAcquiring model parameters from " << numPSMs_forModeling << " spectra with Score >= "
		 << g_model_prob << endl;

	// We collect model parameters by charge state
	zN = 0;
	maxChargeState = getMaxChargeState();



	// There is a minimum number of PSMs we need to obtain accurate modeling parameters.
	// If we don't get at least 'N' PSMs for a given charge state, we can't model
	// that charge state. The user can force the program to model them.
	for(curChargeState = 2; curChargeState <= maxChargeState; curChargeState++) {

		zN = 0;
		for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
			if( curPSM->useForModeling() ) {
				if(curPSM->getCharge() == curChargeState) zN++;
			}
		}
		chargeFreqMap[ curChargeState ] = zN;
	}

	for(cfm_iter = chargeFreqMap.begin(); cfm_iter != chargeFreqMap.end(); cfm_iter++) {
		cerr << "+" << cfm_iter->first << ":";
		cerr << setw(5);
		cerr << cfm_iter->second << " PSMs for modeling.\n";
	}
	cerr << setw(0);



	for(curChargeState = 2; curChargeState <= maxChargeState; curChargeState++) {
		zN = 0; // keep track of how many "modeling" PSMs you encounter

		M_ints_Y.clear();
		M_ints_B.clear();
		M_dist_Y.clear();
		M_dist_B.clear();
		U_ints.clear();
		U_dist.clear();

		for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

			if(curPSM->getCharge() != curChargeState) continue;
			if( !curPSM->useForModeling() ) continue;
			zN++;
			TP.schedule( boost::bind(&PSMClass::threaded_recordModelingParameters_matched, boost::ref(*curPSM) ));
		}
		TP.wait(); // wait for all the threads to end


		// record the matched distances and intensities
		for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

			if(curPSM->getCharge() != curChargeState) continue;
			if( !curPSM->useForModeling() ) continue;

			// matched peaks
			ptr = curPSM->getParamList('m', 'y', 'i');
			for(L = ptr->begin(); L != ptr->end(); L++) M_ints_Y.push_back(*L);

			ptr = curPSM->getParamList('m', 'b', 'i');
			for(L = ptr->begin(); L != ptr->end(); L++) M_ints_B.push_back(*L);

			ptr = curPSM->getParamList('m', 'y', 'd');
			for(L = ptr->begin(); L != ptr->end(); L++) M_dist_Y.push_back(*L);

			ptr = curPSM->getParamList('m', 'b', 'd');
			for(L = ptr->begin(); L != ptr->end(); L++) M_dist_B.push_back(*L);

			// Unmatched peaks
			ptr = curPSM->getParamList('u', 'u', 'i');
			for(L = ptr->begin(); L != ptr->end(); L++) U_ints.push_back(*L);

			ptr = curPSM->getParamList('u', 'u', 'd');
			for(L = ptr->begin(); L != ptr->end(); L++) U_dist.push_back(*L);


		}

		if(zN >= g_MIN_MODEL_NUM) { // record modeling parameters only if you have data for the current charge state

			if(g_DEBUG_MODE) {

				debug_modelData.open("modelData.debug", ios::out);
				debug_modelData << "dataType\tvalue\n";

				for(L = M_ints_B.begin(); L != M_ints_B.end(); L++) debug_modelData << "M_I_b\t" << *L << endl;
				for(L = M_ints_Y.begin(); L != M_ints_Y.end(); L++) debug_modelData << "M_I_y\t" << *L << endl;
				for(L = M_dist_B.begin(); L != M_dist_B.end(); L++) debug_modelData << "M_D_b\t" << *L << endl;
				for(L = M_dist_Y.begin(); L != M_dist_Y.end(); L++) debug_modelData << "M_D_y\t" << *L << endl;

				for(L = U_ints.begin(); L != U_ints.end(); L++) debug_modelData << "U_I\t" << *L << endl;
				for(L = U_dist.begin(); L != U_dist.end(); L++) debug_modelData << "U_D\t" << *L << endl;

				debug_modelData.close();
			}


			// intensity
			meanM_b = getMean(&M_ints_B);
			meanM_y = getMean(&M_ints_Y);
			meanU = getMean(&U_ints);

			varM_b  = getVar(&M_ints_B);
			varM_y  = getVar(&M_ints_Y);
			varU  = getVar(&U_ints);


			// distance
			const double CID_ADJUST = 16.0/25.0; //hwchoi adjustment factor 4/5 to std. dev

			meanMd_b = getMean(&M_dist_B);
			varMd_b  = getVar(&M_dist_B) * CID_ADJUST;

			meanMd_y = getMean(&M_dist_Y);
			varMd_y  = getVar(&M_dist_Y) * CID_ADJUST;

			meanUd = 0;
			varUd  = getVar(&U_dist);



			cerr << endl
				 << "Z: +" << curChargeState << "    # PSM: " << zN << " with Score >= " << g_model_prob << endl
				 << "------------------------------------------------------------------------\n"
				 << "Z: +" << curChargeState << "    b-ion Intensity Matched (mean, stdev): ("
				 << meanM_b << ", " << sqrt(varM_b) << "); N = " << M_ints_B.size() << endl
				 << "Z: +" << curChargeState << "    y-ion Intensity Matched (mean, stdev): ("
				 << meanM_y << ", " << sqrt(varM_y) << "); N = " << M_ints_Y.size() << endl
				 << "Z: +" << curChargeState <<"     b-ion Distance Matched (mean, stdev): ("
				 << meanMd_b << ", " << sqrt(varMd_b) << "); N = " << M_dist_B.size() << endl
				 << "Z: +" << curChargeState <<"     y-ion Distance Matched (mean, stdev): ("
				 << meanMd_y << ", " << sqrt(varMd_y) << "); N = " << M_dist_Y.size() << endl
				 << "Z: +" << curChargeState << "        Intensity Unmatched (mean, stdev): (" << meanU
				 << ", " << sqrt(varU) << "); N = " << U_ints.size() << endl
				 << "Z: +" << curChargeState << "         Distance Unmatched (mean, stdev): (" << meanUd
				 << ", " << sqrt(varUd) << "); N = " << U_dist.size() << endl << endl;

			logF << endl << "Model parameters (Z = " << curChargeState << ")\n"
				 << "# PSMs:\t" << zN << endl
				 << "b-ion intensity Matched (mean, sd):\t(" << meanM_b << ", " << sqrt(varM_b) << "); N = " << M_ints_B.size() << endl
				 << "y-ion intensity Matched (mean, sd):\t(" << meanM_y << ", " << sqrt(varM_y) << "); N = " << M_ints_Y.size() << endl
				 << "b-ion distance Matched (mean, sd):\t(" << meanMd_b << ", " << sqrt(varMd_b) << "); N = " << M_dist_B.size() << endl
				 << "y-ion distance Matched (mean, sd):\t(" << meanMd_y << ", " << sqrt(varMd_y) << "); N = " << M_dist_Y.size() << endl
				 << "Unmatched peak Intensity (mean, sd):\t(" << meanU << ", " << sqrt(varU) << "); N = " << U_ints.size() << endl
			     << "Unmatched peak Distance (mean, sd):\t(" << meanUd << ", " << sqrt(varUd) << "); N = " << U_dist.size() << endl;

			curParams = new modelParamStruct();
			curParams->matched_mean_b = meanM_b;
			curParams->matched_var_b = varM_b;

			curParams->matched_mean_y = meanM_y;
			curParams->matched_var_y = varM_y;

			curParams->unMatched_mean = meanU;
			curParams->unMatched_var = varU;

			curParams->matched_dist_mean_b = meanMd_b;
			curParams->matched_dist_var_b = varMd_b;

			curParams->matched_dist_mean_y = meanMd_y;
			curParams->matched_dist_var_y = varMd_y;

			// these two are only used for CID data
			curParams->unMatched_dist_mean = meanUd;
			curParams->unMatched_dist_var = varUd;


			// record modeling parameters into each PSM
			for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
				if(curPSM->getCharge() != curChargeState) continue;
				curPSM->setParamStruct(curParams);
			}

			// record modeling parameters into modelParametersMap for this pepXMLClass object
			modelParametersMap[ curChargeState ] = *curParams;
		}
		else {
			missedChargeStates.insert(curChargeState);
		}

		// just some clean up
		M_ints_B.clear();
		M_ints_Y.clear();
		M_dist_B.clear();
		M_dist_Y.clear();
		U_ints.clear();
		U_dist.clear();

		delete(curParams); curParams = NULL;

	} // end model building over given charge state

	if(g_captureChargeStateModel) {
		if(!missedChargeStates.empty()) cerr << endl; // prettier output

		// iterate over the charge states and identify any that were not recorded
		// in the modelParametersMap. For those cases, assign them the modeling parameters
		// of a lower charge state. This is not a great solution, but it's currently the
		// best option if you don't have enough data to model a given charge state.
		for(s_iter = missedChargeStates.begin(); s_iter != missedChargeStates.end(); s_iter++) {

			for(int Z = (*s_iter - 1); Z >= 2; Z--) {
				m = modelParametersMap.find( Z );

				if(m != modelParametersMap.end()) { // you found the nearest modeled charge state
					curParams = new modelParamStruct();
					*curParams = m->second;

					for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
						if(curPSM->getCharge() == *s_iter) curPSM->setParamStruct( curParams );
					}
					delete(curParams); curParams = NULL;

					cerr << "+" << *s_iter << " PSMs will be scored using parameters from "
						 << "+" << Z << " PSMs\n";
					break; // exit the loop for this missed charge state.
				}
			}
		}

	}
	cerr << endl; // prettier output


}



// This function gets the modeling parameters needed to score the spectra.
// Use this for HCD data only.
void PepXMLClass::acquireModelParameters_HCD() {

	deque<PSMClass>::iterator curPSM;
	deque<PSMClass> *decoyDeq = NULL;
	map<int, double>::iterator cur_mz_err;
	list<double> M_ints, M_dist, U_ints, U_dist, M_ints_b, M_ints_y;
	list<double> *ptr  = NULL;
	list<double>::iterator L;
	modelParamStruct *curParams = NULL;
	PSMClass *decoyPSM = NULL;


	double meanM_ints, varM_ints, meanU_ints, varU_ints;
	double meanM_dist, varM_dist;
	double mz_err_tol;
	double x;
	int i, zN;
	string msg;
	fstream debug_modelData;


	// A pool of reusable threads that will process the data.
	// The default number of threads is 1
	boost::threadpool::pool TP( g_NUM_THREADS );

	// First step is to identify all of the matched peaks from all of the spectra
	// and collect the intensities of the classified peaks.
	// We will only use spectra with a probability >= modelProb
	// to compute the model parameters.
	i = 0;
	cerr << "\nAcquiring model parameters from " << numPSMs_forModeling << " spectra with Score >= "
		 << g_model_prob << endl;


	zN = 0; // keep track of how many "modeling" PSMs you encounter

	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

		if( !curPSM->useForModeling() ) continue;
		zN++;
		TP.schedule( boost::bind(&PSMClass::threaded_recordModelingParameters_matched, boost::ref(*curPSM) ));
	}
	TP.wait(); // wait for all the threads to end


	// record the matched distances and intensities
	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

		if( !curPSM->useForModeling() ) continue;


		// matched peak data
		ptr = curPSM->getParamList('m', 'y', 'i');
		for(L = ptr->begin(); L != ptr->end(); L++) {
			M_ints.push_back(*L);
			M_ints_y.push_back(*L);
		}

		ptr = curPSM->getParamList('m', 'b', 'i');
		for(L = ptr->begin(); L != ptr->end(); L++) {
			M_ints.push_back(*L);
			M_ints_b.push_back(*L);
		}

		ptr = curPSM->getParamList('m', 'y', 'd');
		for(L = ptr->begin(); L != ptr->end(); L++) M_dist.push_back(*L);

		ptr = curPSM->getParamList('m', 'b', 'd');
		for(L = ptr->begin(); L != ptr->end(); L++) M_dist.push_back( *L );


		// Unmatched peak data
		ptr = curPSM->getParamList('u', 'u', 'i');
		for(L = ptr->begin(); L != ptr->end(); L++) U_ints.push_back(*L);

		ptr = curPSM->getParamList('u', 'u', 'd');
		for(L = ptr->begin(); L != ptr->end(); L++) U_dist.push_back( *L );
	}




	if(zN >= g_MIN_MODEL_NUM) { // record modeling parameters only if you have data for the current charge state

		if(g_DEBUG_MODE) {
			// Debug function that prints out the values used for modeling.
			debug_modelData.open("modelData.debug", ios::out);
			debug_modelData << "dataType\tvalue\n";

			for(L = M_ints.begin(); L != M_ints.end(); L++) debug_modelData << "M_I\t" << *L << endl;
			for(L = M_dist.begin(); L != M_dist.end(); L++) debug_modelData << "M_D\t" << *L << endl;
			for(L = U_ints.begin(); L != U_ints.end(); L++) debug_modelData << "U_I\t" << *L << endl;
			for(L = U_dist.begin(); L != U_dist.end(); L++) debug_modelData << "U_D\t" << *L << endl;

			debug_modelData.close();
		}


		// intensity
		meanM_ints = getMean(&M_ints);
		meanU_ints = getMean(&U_ints);

		varM_ints = getVar(&M_ints);
		varU_ints = getVar(&U_ints);

		// distance
		meanM_dist = getMode(&M_dist);
		varM_dist  = getVar(&M_dist);

		double meanU_dist = getMode(&U_dist);
		double varU_dist  = getVar(&U_dist);

		cerr << "\n# PSM: " << zN << " with Score >= " << g_model_prob << endl;


		cerr << "Distance Matched  (mode, stdev):   (" <<  meanM_dist << ", " << sqrt(varM_dist) << "); N = " << M_dist.size() << endl
			 << "Intensity Matched (mean, stdev):   (" <<  meanM_ints << ", " << sqrt(varM_ints) << "); N = " << M_ints.size() << endl
			 << "Intensity Unmatched (mean, stdev): (" << meanU_ints << ", " << sqrt(varU_ints) << "); N = " << U_ints.size() << endl;


		//cerr << "Distance Unmatched (mean, stdev):  (" << meanU_dist << ", " << sqrt(varU_dist) << "); N = " << U_dist.size() << endl;


		logF << endl
			 << "# PSM: " << zN << " with Score >= " << g_model_prob << endl
			 << "Intensity Matched (mean, stdev): (" <<  meanM_ints << ", " << sqrt(varM_ints) << "); N = " << M_ints.size() << endl
			 << "Distance Matched  (mean, stdev): (" <<  meanM_dist << ", " << sqrt(varM_dist) << "); N = " << M_dist.size() << endl
			 << "Intensity Unmatched (mean, stdev): (" << meanU_ints << ", " << sqrt(varU_ints) << "); N = " << U_ints.size() << endl;
			 //<< "Distance Unmatched (mean, stdev):  (" << meanU_dist << ", " << sqrt(varU_dist) << "); N = " << U_dist.size() << endl;


		curParams = new modelParamStruct();
		curParams->matched_mean_b = meanM_ints;
		curParams->matched_var_b = varM_ints;
		curParams->matched_mean_y = meanM_ints;
		curParams->matched_var_y = varM_ints;

		curParams->unMatched_mean = meanU_ints;
		curParams->unMatched_var = varU_ints;

		curParams->matched_dist_mean_b = meanM_dist;
		curParams->matched_dist_var_b = varM_dist;
		curParams->matched_dist_mean_y = meanM_dist;
		curParams->matched_dist_var_y = varM_dist;

		// these two are only used for CID data so just set them to zero here
		curParams->unMatched_dist_mean = 0;
		curParams->unMatched_dist_var = 0;


		/**********************************************************************
		 * use nonparametric model to define the matched peak distribution
		 * (Hyungwon's code)
		**********************************************************************/
		estimateNonparamInt_b(&M_ints_b, curParams, (g_MZ_ERR * 0.5));
		estimateNonparamInt_y(&M_ints_y, curParams, (g_MZ_ERR * 0.5));
		estimateNonparamInt_U(&U_ints, curParams, (g_MZ_ERR * 0.5));

		estimateNonparamDist(&M_dist, curParams, (g_MZ_ERR * 0.5) );
		estimateNonparamDist_U(&U_dist, curParams, (g_MZ_ERR * 0.5) );
		/**********************************************************************
		 * End nonparametric code
		**********************************************************************/


		// record modeling parameters into each PSM
		for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) curPSM->setParamStruct(curParams);

	}

	if(g_DEBUG_MODE) printDensity(curParams);

	delete(curParams); curParams = NULL;

	cerr << endl; // prettier output


}





// This function does the final scoring
void PepXMLClass::scoreSpectra() {

	if(g_scoreSelect) {
		cerr << "\nScoring " << g_PSMscoreSet.size() << " phospho spectra...\n";
	}
	else  cerr << "\nScoring " << PSMvec->size() << " phospho spectra...\n";

	deque<PSMClass>::iterator curPSM;

	g_progressCtr = 0;
	// A pool of reusable threads that will process the data.
	// The default number of threads is 1
	boost::threadpool::pool TP( g_NUM_THREADS );


	// now use the model parameters to score the spectra
	cerr << "Scoring...  ";
	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

		if(g_scoreSelect) { // consider only PSMs in user list
			set<string>::iterator s;
			s = g_PSMscoreSet.find(curPSM->getSpecId());
			if(s == g_PSMscoreSet.end()) continue;
		}

		TP.schedule(boost::bind( &PSMClass::threaded_scorePSM, boost::ref(*curPSM) ));
	}
	TP.wait(); // wait for all jobs in threadpool queue to finish.

	cerr << "Done.\n\n";

}



// Function actually writes the luciphor results to disk
void PepXMLClass::writeLuciphorResults() {
	deque<PSMClass>::iterator curPSM;

	// consolidate spectra that match to the same phospho peptide sequence
	if( !g_FULL_MONTY && !g_SITE_LEVEL_SCORING ) consolidatePSMs();

	cerr << "Writing results to '" << g_outputName << "'\n";

	// write results to disk
	ofstream outf;
	outf.open( g_outputName.c_str(), ios::out );
	if( !outf ) {
		cerr << "\nERROR: Unable to create output file '" << g_outputName << "'\n"
			 << "Exiting...\n";
		exit(-1);
	}


	if( g_FULL_MONTY ) {

		outf << "specId\t"
			 << "peptideSeq\t"
			 << "predictedPep_1\t"
			 << "predictedPep_2\t";

		if(g_scoringMethod == 0) outf << "pepProphetProb\t";
		if(g_scoringMethod == 1) outf << "Xcorr\t";
		if(g_scoringMethod == 2) outf << "negLogEvalue\t";
		if(g_scoringMethod == 3) outf << "IonScore\t";


		outf << "numRPS\t" // number of reported phospho sites
			 << "numPPS\t" // number of potential phospho sites

			 << "localFLR\t"
			 << "globalFLR\t"
			 << "delta_score\t"
			 << "Luciphor_Score_1\t"
			 << "Luciphor_Score_2\t"

			 << "isDecoy1\t"
			 << "isDecoy2\t"
			 << "totalNumPeaks\t"
			 << "numMatchedPeaks1\t"
			 << "numMatchedPeaks2\t"
			 << "fractionIonsMatched1\t"
			 << "scoreTime"
			 << endl;

	}
	else if( g_SITE_LEVEL_SCORING ) {
		outf << "specId\t";

		if(g_scoringMethod == 0) outf << "pepProphetProb\t";
		if(g_scoringMethod == 1) outf << "Xcorr\t";
		if(g_scoringMethod == 2) outf << "negLogEvalue\t";
		if(g_scoringMethod == 3) outf << "IonScore\t";

		outf << "numRPS\t"
			 << "numPPS\t"
			 << "predictedPeptide\t"
			 << "deltaScore\t"
			 << "siteScores\n";
	}
	else { // default output
		outf << "repSpecId\t"
			 << "origTPPpep\t"
			 << "predictedPep_1\t"
			 << "predictedPep_2\t";

		if(g_scoringMethod == 0) outf << "pepProphetProb\t";
		if(g_scoringMethod == 1) outf << "Xcorr\t";
		if(g_scoringMethod == 2) outf << "negLogEvalue\t";
		if(g_scoringMethod == 3) outf << "IonScore\t";

		outf << "nss\t"
			 << "numRPS\t" // number of reported phospho sites
			 << "numPPS\t" // number of potential phospho sites
			 << "localFLR\t"
			 << "globalFLR\t"
			 << "delta_score\t"
			 << "Luciphor_score_1\t"
			 << "Luciphor_score_2\t"
			 << "isDecoy1\t"
			 << "isDecoy2\t"
			 << "totalNumPeaks\t"
			 << "numMatchedPeaks1\t"
			 << "numMatchedPeaks2"
			 << endl;
	}


	cerr << "Writing results to disk... ";
	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

		if(g_scoreSelect) { // consider only PSMs in user list
			set<string>::iterator s;
			s = g_PSMscoreSet.find(curPSM->getSpecId());
			if(s == g_PSMscoreSet.end()) continue;
		}

		if(g_writeDTA) curPSM->writeSpectrumToDisk();

		curPSM->write_results( outf );
	}
	cerr << "Done.\n\n";

	outf.close();
}



// Function returns the maximum charge state observed in the PSMVec vector
int PepXMLClass::getMaxChargeState() {
	deque<PSMClass>::iterator curPSM;
	int ret = 0;

	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++)
		if(curPSM->getCharge() > ret) ret = curPSM->getCharge();

	return ret;
}



// Function opens an output stream for writing that will record the search parameters
// and statistical data for this run of Luciphor
void PepXMLClass::openLogFile() {

	//define the log file's name and open the file stream for writing
	string logName = "luciphor";
	size_t found = g_outputName.find_last_of(".");
	if(found != string::npos) logName = g_outputName.substr(0, found);
	logName += ".luciphor-" + getTimeStamp() + ".log";


	// change the name of the log file if the program is running Ascore
	if(g_runAscoreInstead) {
		logName = "ascore";
		found = g_outputName.find_last_of(".");
		if(found != string::npos) logName = g_outputName.substr(0, found);
		logName += ".ascore-" + getTimeStamp() + ".log";
	}

	logF.open(logName.c_str(), ios::out);
	if(!logF) {
		cerr << "\nERROR: Unable to create log file '" << logName << "'\n";
		exit(-1);
	}

	// report basic info into the log file
	 logF << "Luciphor BUILD: " << g_BUILD_TIME << endl
	 	 << "Input file:\t" << fileName << endl
		 << "# PSMs:\t" << numPSM << endl;

	logF << getExecutionParameters();

}



// Function closes the log file
void PepXMLClass::closeLogFile(string timeTxt) {
	logF << "\nTotal runtime (HH:MM::SS.fffff):\t" << timeTxt << endl << endl;
	logF.close();
}



// Function consolidates PSMClass objects that share the same bestScore sequence
void PepXMLClass::consolidatePSMs() {

	deque<PSMClass>::iterator curPSM;
	scoreStruct bestScore, curBestScore;
	map<string, repPSMStruct>::iterator mIter;
	repPSMStruct *x = NULL;

	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
		curBestScore = curPSM->getScoreStruct();

		mIter = repPSM_map.find( curBestScore.seq );
		if(mIter == repPSM_map.end()) { // new entry

			x = new repPSMStruct;
			x->repSpecId = curPSM->getSpecId();
			x->bestScore = curBestScore.score;
			x->allSpecIds.push_back( curPSM->getSpecId() );
			repPSM_map[ curBestScore.seq ] = *x;
			delete(x); x = NULL;
		}
		else { // repeat sequence

			if(curBestScore.score > mIter->second.bestScore) {
				mIter->second.bestScore = curBestScore.score;
				mIter->second.repSpecId = curPSM->getSpecId();
			}

			mIter->second.allSpecIds.push_back( curPSM->getSpecId() );
		}
	}


	// now prune the data in PSMvec
	cerr << "Reducing search space... " << PSMvec->size() << " --> ";
	int i = 0;
	deque<PSMClass> keepers;
	map<string, repPSMStruct>::iterator curRepPSM;
	string repId, curId;


	for(curRepPSM = repPSM_map.begin(); curRepPSM != repPSM_map.end(); curRepPSM++) {
		repId = curRepPSM->second.repSpecId;
		i = (signed) curRepPSM->second.allSpecIds.size();

		for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
			curId = curPSM->getSpecId();
			if( repId.compare(curId) == 0 ) {
				curPSM->setNSS( i );
				keepers.push_back( *curPSM );
			}
		}
	}
	cerr << keepers.size() << endl;

	PSMvec->clear();
	*PSMvec = keepers;
	keepers.clear();

}




// Function processed the data using Gygi's Ascore algorithm
void PepXMLClass::process_with_Ascore() {

	deque<PSMClass>::iterator curPSM, nextPSM;
	g_totalNumPSM = (signed) PSMvec->size();
	g_progressCtr = 0;


	// A pool of reusable threads that will process the data.
	// The default number of threads is 1
	boost::threadpool::pool TP( g_NUM_THREADS );


	cerr << "\nRunning Ascore... ";
	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

		if(g_scoreSelect) { // consider only PSMs in user list
			set<string>::iterator s;
			s = g_PSMscoreSet.find(curPSM->getSpecId());
			if(s == g_PSMscoreSet.end()) continue;
		}

		curPSM->setSpectrumPtr("raw");
		curPSM->generatePermutations();

		TP.schedule(boost::bind( &PSMClass::runAscore, boost::ref(*curPSM) ));
	}
	TP.wait(); // wait for all jobs in threadpool queue to finish

	cerr << "Done.\n\n";

	// write results to disk
	ofstream outf;
	outf.open( g_outputName.c_str(), ios::out );

	if(!outf) {
		cerr << "\nERROR: Unable to create output file '" << g_outputName << "'\n"
			 << "Exiting...\n";
		exit(-1);
	}
	cerr << "Writing Ascore results to '" << g_outputName << "'\n\n";

	outf << "specId\t"
		 << "origPeptideSeq\t"
		 << "AscoreSeq1\t"
		 << "AscoreSeq2\t";

	if(g_scoringMethod == 0) outf << "pepProphetProb\t";
	if(g_scoringMethod == 1) outf << "Xcorr\t";
	if(g_scoringMethod == 2) outf << "negLogEvalue\t";
	if(g_scoringMethod == 3) outf << "IonScore\t";

	outf << "numRPS\t" // number of reported phopsho sites
		 << "numPPS\t" // number of potential phospho site (ie: # STY characters)
		 << "peakDepth\t"
		 << "PeptideScore\t"
		 << "Ascore\t"
		 << "totalNumPeaks\t"
		 << "numMatchedPeaks1\t"
		 << "numMatchedPeaks2\n";

	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

		if(g_scoreSelect) { // consider only PSMs in user list
			set<string>::iterator s;
			s = g_PSMscoreSet.find(curPSM->getSpecId());
			if(s == g_PSMscoreSet.end()) continue;
		}

		if(g_writeDTA) curPSM->writeAscoreSpectrum();
		curPSM->write_ascore_results( outf );
	}

	outf.close();

}



// Function computes the False Localization Rate (FLR) for phospho-site assigments
// for the PSMs of this file
void PepXMLClass::calcFLR() {

	std::deque<PSMClass>::iterator curPSM;
	std::deque<flrStruct> *ptr = new std::deque<flrStruct>;
	double maxDeltaScore = -1.0;
	flrStruct *curFLR = NULL;
	FLRClass *flr = NULL;


	if(g_scoreSelect) {
		int N = 0;
		set<string>::iterator S;
		for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
			S = g_PSMscoreSet.find( curPSM->getSpecId() );
			if( S != g_PSMscoreSet.end() ) N++;
		}

		cerr << "\nScored PSMs = " << N << endl;
		if(N < g_MIN_NUM_PSM_FOR_FLR) {
			cerr << "Not enough PSMs scored to accurately estimate the FLR. "
				 << "(Minimum is " << g_MIN_NUM_PSM_FOR_FLR << ")\n"
				 << "FLR estimation will not be performed.\n\n";


			curFLR = new flrStruct;
			curFLR->globalFLR = -1;
			curFLR->localFLR  = -1;
			curFLR->prob      = -1;
			for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++)
				curPSM->setFLR( curFLR );

			delete(curFLR); curFLR = NULL;

			return;
		}
	}

	cerr << "Estimating False Localization Rate (FLR).\n";

	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
		curFLR = new flrStruct;
		*curFLR = curPSM->getFLRdata();
		ptr->push_back( *curFLR );
		if( curFLR->deltaScore > maxDeltaScore ) maxDeltaScore = curFLR->deltaScore;

		delete(curFLR); curFLR = NULL;
	}

	flr = new FLRClass(ptr, maxDeltaScore);
	flr->evalTickMarks();
	flr->calcBothFDRs();
	flr->setMinorMaps();
	flr->performMinorization();
	flr->assignFDRvalues();

	// get the FLR values out of the FLR object
	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
		curFLR = new flrStruct;
		*curFLR = flr->getFLR( curPSM->getSpecId() );
		curPSM->setFLR( curFLR );
		delete(curFLR);
	}

	delete(flr);
	delete(ptr);
}



// Function computes a new m/z error tolerance for the given data set based upon
// the matched peak distances. The g_MZ_ERR variable gets updated here
void PepXMLClass::adjustMZerr() {

	deque<PSMClass>::iterator curPSM;
	list<double>::iterator L;
	list<double> *ptr = NULL;
	deque<double> distDeque;
	deque<double>::iterator D;
	double N, isotopicAdj;
	double confInt = 0.99;

	double alpha = 1.0 - confInt;
	double alpha2 = alpha * 0.5;
	double x1, x2, new_mz_err;
	int i, j;


	// A pool of reusable threads that will process the data.
	// The default number of threads is 1
	boost::threadpool::pool TP( g_NUM_THREADS );


	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

		if( !curPSM->useForModeling() ) continue;

		TP.schedule( boost::bind(&PSMClass::threaded_recordModelingParameters_matched, boost::ref(*curPSM) ));
	}
	TP.wait(); // wait for all the threads to end

	// record the matched distances retaining only distances
	// that are less than 'max_mz_err'
	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {

		if( !curPSM->useForModeling() ) continue;

		ptr = curPSM->getParamList('m', 'y', 'd');
		for(L = ptr->begin(); L != ptr->end(); L++) {
			if( !isInfinite(*L) && !dbl_isnan(*L) )  distDeque.push_back(*L);
		}

		ptr = curPSM->getParamList('m', 'b', 'd');
		for(L = ptr->begin(); L != ptr->end(); L++) {
			if( !isInfinite(*L) && !dbl_isnan(*L) ) distDeque.push_back(*L);
		}
		ptr = NULL;
		x1 = 0;

		curPSM->clear(); // prep this PSM object for the next scoring iteration
	}

	if( !distDeque.empty() ) {

		for(i = 0; i < (signed)distDeque.size(); i++) {
			cout << distDeque.at(i) << endl;
		}
		exit(0);
	}
}




// Function resets the mz_err window to the value of g_MZ_ERR
void PepXMLClass::resetMZerr() {
	deque<PSMClass>::iterator curPSM;

	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
		curPSM->updateMZerr( (g_MZ_ERR * 0.5) );
	}
}


//// Function just prints out the number of permutations a given PSM will have
//void PepXMLClass::getCounts() {
//	deque<PSMClass>::iterator curPSM;
//
//	cout << "specId\tpeptide\tpepLen\tnumSTY\tnumPhosphoSites\tnumPerm\tnumDecoys\n";
//	for(curPSM = PSMvec->begin(); curPSM != PSMvec->end(); curPSM++) {
//		cout << curPSM->countInfo();
//	}
//	exit(0);
//}

