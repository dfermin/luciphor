/*
 * globals.cpp
 *
 *  Created on: Mar 18, 2011
 *      Author: dfermin
 */

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <map>
#include <deque>
#include <ctime>
#include <fstream>
#include <sstream>
#include <ctype.h>
#include <algorithm>
#include <getopt.h>  // required for getOpt_long() function
#include <cmath> // needed for power function
#include <boost/filesystem/operations.hpp>
#include <boost/thread.hpp> // for multi-threading
#include "globals.hpp"


using namespace std;


string g_scoringMethod = "default";
string g_srcXMLfile;
string g_ext = "mzXML";
string g_srcDir;
string g_outputName;
string g_cmdLineArgs;
string g_BUILD_TIME = "";

bool g_runAscoreInstead = false;
bool g_userDefinedOutput = false;
bool g_FULL_MONTY = false;
bool g_IS_HCD = false;
bool g_writeDTA = false;
bool g_DEBUG_MODE = false;
bool g_NO_NL_PEAKS = false;
bool g_NO_NL_PEAKS_MODEL = false;
bool g_NL_MODEL_ONLY = false;
bool g_randDecoyAA = true; // by default we search decoys
bool g_singleLetter = false;
bool g_captureChargeStateModel = false;
bool g_LIMIT_CHARGE_STATE = false;
bool g_WRITE_TOP_TWO = false;
bool g_useOnlySiteDetermIons = false;
bool g_usePPM = false;

double MIN_MZ = 100.0; // lowest m/z value we'll consider
double g_MZ_ERR = 0.5;
double g_prob_threshold = 0.05;
double g_model_prob = 0.95;
double g_NUM_PERMS_LIMIT = pow(2.0, 14); // maximum number of permutations to consider

set<string> g_PSMscoreSet;
bool g_scoreSelect = false;

int g_NUM_THREADS = 1;
int g_intensityType = 2;
int g_CHARGE_STATE = 0;
int g_MIN_MODEL_NUM = 10;

map<char, char> decoyAA;
map<char, double> AAmass;
map<char, string> modAAmap;


// List long ( double dashed  ) command line options here
// For tutorial see: http://www.ibm.com/developerworks/aix/library/au-unix-getopt.html
static const struct option longOpts[] = {
		{ "Ascore", no_argument, NULL, 'A' },
		{ "single", no_argument, NULL, 0 },
		{ "capture", no_argument, NULL, 'c' },
		{ "noDecoys", no_argument, NULL, 0 },
		{ "hcd", no_argument, NULL, 0 },
		{ "siteDetermIons", no_argument, NULL, 0 },
        { "ppm", no_argument, NULL, 0 },
		{ "debug", no_argument, NULL, 0 },
		{ "Sl", required_argument, NULL, 0 }
};



void print_usage() {
	cerr << "\nUSAGE: luciphor -i -d -w [ -e -f -m -p -o (-M -P) -t -b -T -A -n -c -Z -k]\n"
		 << "   -i <interact.pep.xml>     pepXML file to analyze\n"
		 << "   -d <path>                 the path to the spectral data files\n"
		 << "   -w <# greater than 0>     The MS2 fragment ion tolerance in Daltons (default is " << g_MZ_ERR << " Da.)\n"

		 << "\nOptional parameters\n"
		 << "   -p <0.0-1.0>              min. Peptide probability threshold (default is " << g_prob_threshold << ")\n"
		 << "   -m <0.0-1.0>              min. Peptide probability to use for parameter estimations (default is " << g_model_prob << ")\n"
		 << "   -e <mzML, mzXML>          File extension used for reading spectral data files (default is mzXML)\n"
		 << "   -f                        \"The Full Monty\": Report scores for all spectra not just representative ones\n"
		 << "   -A / --Ascore             Run our version of Ascore algorithm instead (m/z window fixed at 0.6 for this)\n"
		 << "   -c / --capture            Capture PSM's you can't model, use the model parameters for the next closest charge state\n"
		 << "   -Z <1,2,3...>             Only model and score spectra of the given charge state\n"
		 << "   -o <output_file_name>     Write results to this file name (default output name is based upon input pepXML file name)\n\n"
		 << "   -T <1-N>                  Number of threads to use (default is 1)\n"

		 << "   --hcd                     Tells me that the input data is HCD and to adjust my parameters accordingly\n\n"

		 << "   -k <256-N>                Consider PSMs with 'k' number of permutations (default is " << g_NUM_PERMS_LIMIT << ", minimum is 256)\n"
		 << "                             *Note: increasing this parameter beyond " << g_NUM_PERMS_LIMIT << " requires >8GB of RAM and many hours to run\n\n"



		 << "   -n <1,2,3>                Tells me how to deal with neutral loss fragment ions\n"
		 << "                             By default, neutral loss peaks are used for both modeling and scoring\n"
		 << "                             1 = do not use them at all\n"
		 << "                             2 = use neutral loss peaks ONLY for model parameter acquisition\n"
		 << "                             3 = use them only in final scoring of phospho-peptides (ie: not for model parameter acquisition)\n\n"

		 << "   -t <1,2>                  Tells me to write the spectra to disk with each assigned peak annotated\n"
		 << "                             1 = raw spectrum, 2 = scaled intensities (0-100)\n" //3 = log(quantile normalized intensities)\n"
		 << "   -b                        When writing spectra to disk, write *BOTH* luciphor peptide predictions to disk\n\n"

		 << "   -P <score_method>         Use this option if you want to pick spectra based upon\n"
		 << "                             a search engine score instead of PeptideProphet probabilities\n"
		 << "                             *** You MUST use -M _WITH_ -P ***\n"
		 << "                             -P sequest=<Xcorr>   Use the Sequest XCorr values\n"
		 << "                             -P xtandem=<-log(Evalue)>  Use negative log of X-tandem expect score: -log(Evalue)\n"
		 << "                             -P mascot=<ionScore> Use the Mascot ion-score\n\n"

		 << "   -M <score_method>         Use this option if you want to pick spectra for modeling based upon \n"
		 << "                             a search engine score instead of probabilities.\n"
		 << "                             *** You MUST use -P _WITH_ -M ***\n"
		 << "                             -M sequest=<Xcorr>   Use the Sequest XCorr values\n"
		 << "                             -M xtandem=<-log(Evalue)>  Use negative log of X-tandem expect score: -log(Evalue)\n"
		 << "                             -M mascot=<ionScore> Use the Mascot ion-score\n\n"

		 << "   --Sl <file>               Score and report only results for the PSMs in this file, one PSM ID per line\n"
		 << "   --noDecoys                Do *NOT* estimate False Localization Rate (FLR) using decoy phosphorylation sites (Luciphor only option)\n"
		 << "   --single                  Represent modified residues by a single character (usually the AA letter in lower case)\n"
		 << "   --siteDetermIons          Use only site determining ions in final scoring (Luciphor only option)\n\n"

		 << "   --debug                   Debug mode, for developers. Produces extra files to track down errors.\n"

		 << endl;
}



// function parses the command line options given by the user
void parse_command_line_args(int argc, char *argv[]) {

	if(argc < 5) {
		print_usage();
		exit(0);
	}

	string x;
	for(int i = 1; i < argc; i++) {
		x.clear();
		x.assign( argv[i] );
		g_cmdLineArgs += x + " ";
	}
	x.clear();


	string *scoringStr = NULL; // for -P option
	string *modelingStr = NULL; // for -M option
	string *decoyStr = NULL; // for decoy amino acid
	string PSMfile;

	int nn = 0;
	int c;

	int longIndex;
	while( (c = getopt_long(argc, argv, "m:p:i:w:d:e:T:t:P:M:o:n:Z:k:fAcb", longOpts, &longIndex)) != -1 ) {
		switch(c) {

		case 'k':
			g_NUM_PERMS_LIMIT = atof(optarg);
			break;
		case 'Z':
			g_LIMIT_CHARGE_STATE = true;
			g_CHARGE_STATE = atoi(optarg);
			break;
		case 'R':
			g_randDecoyAA = false;
			break;
		case 's':
			g_singleLetter = true;
			break;
		case 'c':
			g_captureChargeStateModel = true;
			break;
		case 'T':
			g_NUM_THREADS = atoi(optarg);
			break;
		case 'f':
			g_FULL_MONTY = true;
			break;
		case 'n':
			nn = atoi(optarg);
			if(nn == 1) g_NO_NL_PEAKS = true;
			if(nn == 2) g_NL_MODEL_ONLY = true;
			if(nn == 3) g_NO_NL_PEAKS_MODEL = true;
			break;
		case 't':
			g_writeDTA = true;
			g_intensityType = atoi(optarg);
			break;
		case 'b':
			g_WRITE_TOP_TWO = true;
			break;
		case 'm':
			g_model_prob = atof(optarg);
			break;
		case 'p':
			g_prob_threshold = atof(optarg);
			break;
		case 'o':
			g_outputName = optarg;
			g_userDefinedOutput = true;
			break;
		case 'A':
			g_runAscoreInstead = true;
			break;
		case 'i':
			g_srcXMLfile = optarg;
			break;
		case 'w':
			g_MZ_ERR = atof(optarg);
			break;
		case 'd':
			g_srcDir = optarg;
			break;
		case 'e':
			g_ext = optarg;
			break;
		case 'M':
			modelingStr = new string( optarg );
			break;
		case 'P':
			scoringStr = new string( optarg );
			break;

		case 0: /* for long options */
			for(int k = 0; k < 8; k++) {
				if( strcmp("Ascore", longOpts[longIndex].name ) == 0 ) {
					g_runAscoreInstead = true;
				}
				if( strcmp("capture", longOpts[longIndex].name) == 0 ) {
					g_captureChargeStateModel = true;
				}
				if( strcmp("hcd", longOpts[longIndex].name) == 0 ) {
					g_IS_HCD = true;
				}
				if( strcmp("noDecoys", longOpts[longIndex].name) == 0 ) {
					g_randDecoyAA = false;
				}
				if( strcmp("single", longOpts[longIndex].name) == 0 ) {
					g_singleLetter = true;
				}
				if( strcmp("ppm", longOpts[longIndex].name) == 0 ) {
					g_usePPM = true;
				}
				if( strcmp("siteDetermIons", longOpts[longIndex].name) == 0 ) {
					g_useOnlySiteDetermIons = true;
				}
				if( strcmp("debug", longOpts[longIndex].name) == 0 ) {
					g_DEBUG_MODE = true;
				}
				if( strcmp("Sl", longOpts[longIndex].name) == 0 ) {
					g_scoreSelect = true;
					PSMfile = optarg;
				}
			}
			break;
		default:
			print_usage();
			exit(0);
		}
	}

	if(g_MZ_ERR <= 0 ) {
		cerr << "\nERROR: You must specify the fragment ion tolerance with '-w'\n\n";
		exit(0);
	}

	if(g_ext.length() < 1) {
		cerr << "\nERROR: I need to know what type of files I'll be parsing "
			 << "out of! use -e \n\n";
		exit(0);
	}

	if(g_srcXMLfile.length() < 1) {
		cerr << "\nERROR: '-i' is a required input parameter.\n\n";
		exit(0);
	}

	if( g_runAscoreInstead && g_captureChargeStateModel ) {
		cerr << "\nERROR: You can't run -A with -c\n"
			 << "The Ascore algorithm does not use modeling.\n\n";
			exit(0);
	}

	if(g_NUM_PERMS_LIMIT < 256) { // minimum permutation threshold

		cerr << "\nERROR: Option '-k' is below minimum value of 256.\n";
		exit(0);
	}

	if( g_runAscoreInstead ) { // no decoy modeling for Ascore
		g_randDecoyAA = false;
	}

	if( g_WRITE_TOP_TWO && !g_writeDTA) {
		cerr << "\nERROR: -b requires the -t option\n\n";
		exit(0);
	}

	if( g_LIMIT_CHARGE_STATE && (g_CHARGE_STATE < 1) ) {
		cerr << "\nERROR: '" << g_CHARGE_STATE << "' < 1, "
			 << "you must provide a provide a positive integer for this variable.\n\n";
		exit(0);
	}

	// you can't capture neighboring charge states if you chose to score only
	// PSMs of a specific charge state
	if(g_LIMIT_CHARGE_STATE) g_captureChargeStateModel = false;


	// this handles alternative scoring function error checking
	if( (scoringStr != NULL) && (modelingStr == NULL) ) {
		cerr << "\nERROR: -M is required with -P option\n\n";
		print_usage();
		exit(0);
	}
	else if( (scoringStr == NULL) && (modelingStr != NULL) ) {
		cerr << "\nERROR: -P is required with -M option\n\n";
		print_usage();
		exit(0);
	}
	else if( (scoringStr != NULL) && (modelingStr != NULL) ) {
		parse_alternative_scoring(scoringStr, modelingStr);
	}


	if(g_NUM_THREADS > (signed) boost::thread::hardware_concurrency() ) {
		cerr << "\nERROR: Requested number of threads (" << g_NUM_THREADS
			 << ") exceeds your available number of threads: "
			 << (signed) boost::thread::hardware_concurrency()
			 << "\nReseting threads to 1\n";
		g_NUM_THREADS = 1;
	}

	if(g_DEBUG_MODE) { // in order to write out debug files, we need to limit program to 1 thread
		g_NUM_THREADS = 1;
	}


	if(g_IS_HCD) {
		MIN_MZ = 100.0; // we can go pretty low in the m/z scale with HCD data
		g_MIN_MODEL_NUM = 50; // since we don't separate based on charge state, we need to increase this
	}

	if( !g_userDefinedOutput ) {
		// set the output file's name based upon input file's name
		boost::filesystem::path curFile( g_srcXMLfile.c_str() );
		string basename = curFile.leaf();
		size_t found = basename.find('.');
		if(found != string::npos) {
			g_outputName = basename.substr(0, found);
		}

		if(g_runAscoreInstead) g_outputName += "_ascore_results.tsv";
		else g_outputName += "_lucipher_results.tsv";
	}




	/***********************
	 * Report to the user what parameters will be used
	 ***********************/
	cerr << "\n==============================================================\n"
	     << "Run parameters:\n"
		 << "Data type: " << g_scoringMethod << endl
		 << "Modeling threshold >= " << g_model_prob << endl
		 << "Scoring threshold  >= " << g_prob_threshold << endl
		 << "Permutation limit: " << g_NUM_PERMS_LIMIT << endl
		 << "Spectral format:   " << g_ext << endl
		 << "Threads: " << g_NUM_THREADS << endl;

	cerr << "Algorithm: ";
	if(g_IS_HCD) cerr << "HCD Mode\n";
	else if(g_runAscoreInstead) cerr << "A-Score (Modeling thresholds will be ignored)\n";
	else cerr << "CID Mode\n";

	cerr << "Decoys: " << g_randDecoyAA << endl;

	cerr << "Fragment ion tolerance: ";
	if(g_runAscoreInstead) cerr << "0.6 Da\n";
	else {
		cerr << g_MZ_ERR;
		if(g_usePPM) cerr << " PPM\n";
		else cerr << " Da\n";
	}

	cerr << "Output format: ";
	if(g_FULL_MONTY) cerr << "Full Monty\n";
	else cerr << "Default\n";

	cerr << "Output to: " << g_outputName << endl;

	cerr << "\nOther options used (if any):\n";

	if(g_useOnlySiteDetermIons) {
		cerr << "\tOny site determining ions will be used for scoring\n";
	}


	if(g_NO_NL_PEAKS) {
		cerr << "\tNeutral loss fragment ions will be ignored\n";
	}

	if(g_NO_NL_PEAKS_MODEL) {
		cerr << "\tNeutral loss fragment ions will *NOT* be used when acquiring model parameters\n";
	}
	if(g_NL_MODEL_ONLY) {
		cerr << "\tNeutral loss fragment ions will *ONLY* be used for acquiring model parameters\n";
	}

	if(g_LIMIT_CHARGE_STATE) {
		cerr << "\tWill only model and score +" << g_CHARGE_STATE << " PSMs\n";
	}

	if(g_captureChargeStateModel && !g_runAscoreInstead) {
		cerr << "\tWill use nearest neighbor parameters for unmodeled charge states.\n";
	}


	if( g_scoreSelect ) {
		parsePSMfile(PSMfile);
		cerr << "\tOnly the " << g_PSMscoreSet.size() << " PSMs from '"
			 << PSMfile << "' will be scored\n";
	}

	if(g_writeDTA) {
		cerr << "\tWriting matched peaks to disk. Peak format: ";
		if(g_intensityType == 1) cerr << "RAW Intensities";
		if(g_intensityType == 2) cerr << "Scaled Intensities (0-100)";
		if(g_intensityType == 3) cerr << "Median Normalized Intensities";
	}

	if(g_WRITE_TOP_TWO) cerr << "\tWriting matched peaks for BOTH predictions\n";


	if(g_DEBUG_MODE) cerr << "\tRunning in debug mode (Limited to 1 thread)...\n";

	cerr << "==============================================================\n";
}


// Function returns a string that reports all of the options that are being
// used for this execution of luciphor
string getExecutionParameters() {
	string ret;

	ret += "Run Date:\t" + getTimeStamp() + "\n" 
		+ "Scoring method:\t" + g_scoringMethod + "\n"
		+ "Modeling threshold:\t" + dbl2string(g_model_prob) + "\n"
		+ "Scoring threshold:\t" + dbl2string(g_prob_threshold) + "\n"
		+ "Permutation limit:\t" + dbl2string(g_NUM_PERMS_LIMIT) + "\n"
		+ "Number of threads:\t" + int2string(g_NUM_THREADS) + "\n";

	ret += "Fragment ion tolerance:\t " + dbl2string(g_MZ_ERR);
	if(g_usePPM) ret +=  " PPM\n";
	else ret += " Da\n";

	ret += "Output file name:\t" + g_outputName + "\n";

	ret += "Command line options used:\t";

	ret += g_cmdLineArgs + "\n";
	ret += "\n";
	return ret;
}



// Function returns true if the passed string contains only the
// twenty valid amino acid characters.
bool allValidAAs(string *srcPtr) {
	bool ret = true;
	int badCount = 0;

	int N = (signed) srcPtr->length();
	for(int i = 0; i < N; i++) {
		char c = toupper( srcPtr->at(i) );

		switch (c) {
			case 'B':
				badCount++;
				break;

			case 'J':
				badCount++;
				break;

			case 'O':
				badCount++;
				break;

			case 'U':
				badCount++;
				break;

			case 'X':
				badCount++;
				break;

			case 'Z':
				badCount++;
				break;

			default:
				break;
		}
	}

	if(badCount > 0) ret = false;
	return ret;
}




// function parses the alternative scoring string data and records the
// scoring info into global variables
void parse_alternative_scoring(string *scoringStr, string *modelingStr) {

	int iS = 0, iM = 0;
	string xCorr_str, deltaCn_str;


	iS = scoringStr->find('=');
	iM = scoringStr->find('=');

	if(iS == (signed)string::npos){
		cerr << "\nERROR -P option is incorrectly formatted\n\n";
		print_usage();
		exit(0);
	}

	if(iM == (signed)string::npos){
		cerr << "\nERROR -M option is incorrectly formatted\n\n";
		print_usage();
		exit(0);
	}

	string s = scoringStr->substr(0, iS);
	string m = modelingStr->substr(0, iM);

	if(s.compare(m) != 0) {
		cerr << "\nERROR -P and -M must be the same scoring method\n\n";
		print_usage();
		exit(0);
	}

	string method = scoringStr->substr(0, iS);

	string cutoffS, cutoffM;
	cutoffS = scoringStr->substr(iS+1);
	cutoffM = modelingStr->substr(iM+1);

	g_scoringMethod = method;

	// X!Tandem expect scores are read in as -log(expect)
	// this means all score comparisons are in the >= category like Sequest/Mascot/Probabilities
	if(g_scoringMethod.compare("xtandem") == 0) {
		g_prob_threshold = str2dbl(cutoffS);
		g_model_prob = str2dbl(cutoffM);
	}

	if(g_scoringMethod.compare("sequest") == 0) {
		g_prob_threshold = str2dbl(cutoffS);
		g_model_prob = str2dbl(cutoffM);
	}

	if(g_scoringMethod.compare("mascot") == 0) {
		g_prob_threshold = str2dbl(cutoffS);
		g_model_prob = str2dbl(cutoffM);
	}

}




// function takes in a new amino acid to add to the AAmass map
// these new amino acid values come from the pepXML file
void addAAmass(string AA, double mass, char isVar) {
	char ch, newCh;
	map<char, string>::iterator curMod;
	string str, massStr;
	string altChars = "98765432";


	if(isVar == 'Y') ch = tolower(AA[0]); // represent variable mods with lowercase letters
	else ch = AA[0];

	if(isVar == 'Y') {
		curMod = modAAmap.find(ch);
		if(curMod == modAAmap.end()) { // new modification to add
			AAmass[ ch ] = mass;
			str = AA;
			massStr = "[" + int2string( (int) round_dbl(mass,0) ) + "]";
			str += massStr;
			modAAmap[ ch ] = str; //example: q = Q[111]


			if(g_DEBUG_MODE) cerr << "1. Added variable mod '" << ch << "': " << modAAmap[ ch ] << endl;
		}
		else {
			// this is a new modification but the letter that represents it
			// has already been used, we need to replace the letter with an
			// alternate character

			// iterate over each potential character until you find one
			// you haven't used yet. Then associate this variable modification
			// with that character
			for(int i = 0; i < (signed) altChars.length(); i++) {
				newCh = altChars[i];
				curMod = modAAmap.find(newCh);

				if(curMod == modAAmap.end()) {
					AAmass[ newCh ] = mass;
					str = AA;
					massStr = "[" + int2string( (int) round_dbl(mass,0) ) + "]";
					str += massStr;
					modAAmap[ newCh ] = str;
				}
			}

			if(g_DEBUG_MODE) cerr << "2. Added variable mod '" << ch << "': " << modAAmap[ ch ] << endl;
		}
	} // end if(isVar == 'Y')
	else {
		AAmass[ ch ] = mass; // fixed modifications replace existing masses
	}
}


// function initializes the AA mass map
void initialize_AA_masses() {

	string str;
	map<char, char>::iterator curDecoy;
	double mass;
	char c;

    AAmass['A'] = 71.03711;
    AAmass['R'] = 156.1011;
    AAmass['N'] = 114.04293;
    AAmass['D'] = 115.02694;
    AAmass['C'] = 103.00919;
    AAmass['E'] = 129.04259;
    AAmass['Q'] = 128.05858;
    AAmass['G'] = 57.02146;
    AAmass['H'] = 137.05891;
    AAmass['I'] = 113.08406;
    AAmass['L'] = 113.08406;
    AAmass['K'] = 128.09496;
    AAmass['M'] = 131.04049;
    AAmass['F'] = 147.06841;
    AAmass['P'] = 97.05276;
    AAmass['S'] = 87.03203;
    AAmass['T'] = 101.04768;
    AAmass['W'] = 186.07931;
    AAmass['Y'] = 163.06333;
    AAmass['V'] = 99.06841;

    // phospho-AA's
    AAmass['s'] = 87.03203 + 79.9663;
    AAmass['t'] = 101.04768 + 79.9663;
    AAmass['y'] = 163.06333 + 79.9663;

/******************************************************************************/

	// decoy AA modifications
	AAmass['2'] = AAmass['A'] + 79.9663;
	AAmass['3'] = AAmass['R'] + 79.9663;
	AAmass['4'] = AAmass['N'] + 79.9663;
	AAmass['5'] = AAmass['D'] + 79.9663;
	AAmass['6'] = AAmass['C'] + 79.9663;
	AAmass['7'] = AAmass['E'] + 79.9663;
	AAmass['8'] = AAmass['Q'] + 79.9663;
	AAmass['9'] = AAmass['G'] + 79.9663;
	AAmass['0'] = AAmass['H'] + 79.9663;
	AAmass['@'] = AAmass['I'] + 79.9663;
	AAmass['#'] = AAmass['L'] + 79.9663;
	AAmass['$'] = AAmass['K'] + 79.9663;
	AAmass['%'] = AAmass['M'] + 79.9663;
	AAmass['&'] = AAmass['F'] + 79.9663;
	AAmass[';'] = AAmass['P'] + 79.9663;
	AAmass['?'] = AAmass['W'] + 79.9663;
	AAmass['~'] = AAmass['V'] + 79.9663;

	// maps between known amino acid char's and our decoy char set
	decoyAA['A'] = '2';
	decoyAA['R'] = '3';
	decoyAA['N'] = '4';
	decoyAA['D'] = '5';
	decoyAA['C'] = '6';
	decoyAA['E'] = '7';
	decoyAA['Q'] = '8';
	decoyAA['G'] = '9';
	decoyAA['H'] = '0';
	decoyAA['I'] = '@';
	decoyAA['L'] = '#';
	decoyAA['K'] = '$';
	decoyAA['M'] = '%';
	decoyAA['F'] = '&';
	decoyAA['P'] = ';';
	decoyAA['W'] = '?';
	decoyAA['V'] = '~';


	// record the decoy AA character strings
    for(curDecoy = decoyAA.begin(); curDecoy != decoyAA.end(); curDecoy++) {
    	str = curDecoy->first;
    	mass = AAmass[ curDecoy->second ];
    	str += "[" + int2string( (int) round_dbl(mass, 0) ) + "]";

    	c = decoyAA[ curDecoy->first ]; // get the single-character representation of this decoy
    	modAAmap[ c ] = str;
    }


/******************************************************************************/
}



// Function substitutes single letter characters for their full mass-qualified string
string repModAAchar(string *srcPtr) {
	string ret;
	int N = srcPtr->length();
	char curChar;
	map<char, string>::iterator curMod;

	if(g_singleLetter) ret = *srcPtr; // just keep the string as is
	else {
		for(int i = 0; i < N; i++) {
			curChar = srcPtr->at(i);
			curMod = modAAmap.find(curChar);

			if(curMod == modAAmap.end()) ret += curChar; // unmodified AA
			else {
				ret += curMod->second;
			}
		}
	}

	return ret;
}



// Function generates a randomly phosphorylated decoy peptide.
// All non-STY AA characters are considered (so long as they are not already
// modified).
string genRandDecoyPeptide(string srcSeq, int numSites) {
	string seq = upperCaseSTY( srcSeq );
	int seqLen = (signed) srcSeq.length(); // want to limit decoy search space
	char c;
	int N = 0;
	string ret; // return variable
	map<char, char>::iterator decoyIter;


	// Store the srcSeq as an array in a deque
	deque<char> *strAry = NULL;
	strAry = new deque<char>();
	for(int i = 0; i < seqLen; i++) strAry->push_back( seq.at(i) );


	// record the indexes of all candidate residues in 'strAry' that
	// can be used for decoy modeling
	deque<int> *d = NULL;
	d = new deque<int>();

	for(int i = 0; i < seqLen; i++) {
		c = seq.at(i);
		if( (c != 'S') && (c != 'T') && (c != 'Y') && (c != 'X') && (!islower(c)) ) {
		//if( (c != 'S') && (c != 'T') && (c != 'Y') && (c != 'X') ) {

			// if you got this far, that means the character at this index
			// is okay to use for making a decoy phospho-peptide
			d->push_back(i);
		}
	}
	
	N = (signed) d->size(); // find out how many candidate residues you have

	if(N < numSites) {
		// There are not enough candidate residues in this peptide to generate
		// a decoy version of it.
		ret = "no_valid_decoys";
	}
	else {
		for(int i = 0; i < 3; i++) random_shuffle(d->begin(), d->end()); // randomly permute the candidate residues

		for(int j = 0; j < numSites; j++) {
			int p = d->at(j);
			//c = toupper( strAry->at(p) ); // get the character at this coordinate
			c = strAry->at(p); // get the character at this coordinate

			decoyIter = decoyAA.find(c);

			if(decoyIter != decoyAA.end()) strAry->at(p) = decoyIter->second;
			else {
				cerr << "\nERROR!: Unable to generate a decoy sequence for '" << srcSeq << "'\n";
				cerr << "Unable to continue. Exiting now.\n\n";
				exit(0);
			}
		}

		ret = ""; // now reassemble the string from strAry
		for(int k = 0; k < seqLen; k++) ret += strAry->at(k);

	}

	delete(d); d = NULL;
	delete(strAry); strAry = NULL;

	return ret;
}



// Function returns true if the given peptide contains a decoy modification
bool isDecoyPep(string *srcPtr) {
	bool ret = false;
	int N = srcPtr->length();
	char curChar;
	size_t found;

	for(int i = 0; i < N; i++) {
		curChar = srcPtr->at(i);

		if( !isalpha(curChar) ) {
			ret = true;
			break;
		}
	}

	return ret;
}


// function returns the passed string with all letters as uppercase characters
string uc(string src) {
	string ret = "";
	size_t N = src.length();

	for(int i = 0; i < N; i++) ret += toupper(src.at(i));

	return ret;
}



// function to round doubles to N number of decimal places
double round_dbl(double r, int places) {
	double off = pow((double)10, places);
	return ( round(r * off) / off );
}



// function splits a string into an array and returns it
vector<string> split_string(string txt) {
	vector<string> ret;
	string buf; // a buffer string
	stringstream ss(txt); // inserts 'txt' into string_stream 'ss'

	//everytime a whitespace char is encounter, a new split occurs
	while(ss >> buf) ret.push_back(buf);

	return ret;
}



// Function converts an int into a string
string int2string(int i) {
	string ret;

	stringstream ss;
	ss << i;
	ret = ss.str();
	return ret;
}



// Function converts a double into a string
string dbl2string(double d) {
	string ret;

	stringstream ss;
	ss << d;
	ret = ss.str();
	return ret;
}



// function converts a string into a double
double str2dbl(string ch) {
	double ret;

	istringstream iss(ch);
	iss >> ret;

	return ret;
}



// Recursive function to compute the Factorial for the given number
double Factorial(double x) {
	double ret=1;
	double i;

	if( x <= 1 ) return 1;

	for (i=1; i<=x; i++) ret *= i;

	return ret;
}



// Function returns the number of combinations of (n choose k) a given number pair
double combinatorial(double n, double k) {
	double ret = 1;

	// n = number of characters in the alphabet 'set'
	// k = desired length of the words to be constructed

	if(n <= 1)
		ret = 1;
	else {
		double diff = n - k;
		ret = Factorial(n) / ( Factorial(k) * Factorial(diff) );
	}

	return ret;
}



// Function converts all S,T, or Y lower case letters to caps
string upperCaseSTY(string src) {
	string ret;
	ret = src;

	for(int i = 0; i < (signed)src.length(); i++) {
		if(src[i] == 's') ret[i] = 'S';
		else if(src[i] == 't') ret[i] = 'T';
		else if(src[i] == 'y') ret[i] = 'Y';
		else ret[i] = src[i];
	}

	return ret;
}


// Function checks to see if a given double is 'nan'
bool dbl_isnan(double var) {
    volatile double d = var;
    return d != d;
}


// Function returns true if the given string contains a phosphorylated AA
// This means the string contains 's/t/y'
bool containsPhosphoSite(string txt) {
	bool ret = false;
	int ctr = 0;
	string ss;
	size_t found;

	// need to remove the fragment indication prefix (ie: y^# or b^#) first
	found = txt.find(":");
	if(found != string::npos) ss = txt.substr(found+1);
	else ss = txt;

	for(int i = 0; i < (signed)ss.length(); i++) {
		if(ss.at(i) == 's') ctr++;
		else if(ss.at(i) == 't') ctr++;
		else if(ss.at(i) == 'y') ctr++;
		//else if(ss.at(i) == g_decoyAA) ctr++;
	}

	if(ctr > 0) ret = true;

	return ret;
}


// Function returns true if the given string contains an S,T, or Y AA
bool containsSTY(string txt) {
	bool ret = false;
	int ctr = 0;
	string ss;
	size_t found;

	// need to remove the fragment indication prefix (ie: y^# or b^#) first
	found = txt.find(":");
	if(found != string::npos) ss = txt.substr(found+1);
	else ss = txt;

	string decoyChars = "234567890@#$%&;?~";

	for(int i = 0; i < (signed)ss.length(); i++) {
		if(toupper(ss.at(i)) == 'S') ctr++;
		else if(toupper(ss.at(i)) == 'T') ctr++;
		else if(toupper(ss.at(i)) == 'Y') ctr++;
		else {
			// this handles the case of a decoy character
			found = decoyChars.find( ss.at(i) );
			if(found != string::npos) ctr++;
		}
	}

	if(ctr > 0) ret = true;

	return ret;
}





// Function returns a timestamp in the form of a string
string getTimeStamp() {

	time_t rawtime;
	struct tm *timeInfo;

	time(&rawtime);
	timeInfo = localtime(&rawtime);

	int year = 1900 + timeInfo->tm_year;
	int month = timeInfo->tm_mon;
	int day = timeInfo->tm_mday;

	// this appends a zero prefix to the date, should it be < 10
	string dayStr;
	if(day < 10) dayStr = "0" + int2string(day);
	else dayStr = int2string(day);

	string monthStrs[12] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};

	string ret = int2string(year) + monthStrs[ month ] + dayStr;

	return ret;
}





// Function generates all binary combinations of length 'L' having 'k' switches
// flipped from 0 to 1
set< vector<int> > generateBinaryGrid(double L, double k) {

	// number of bit strings of length L you'll need to generate
	double numIter = pow((double) 2, (double) L);

	set< vector<int> > grid;
	vector<int> curSet;

	int N = (int) combinatorial(L, k); // figure out how big the grid will be
	int v_sum;

	curSet.reserve( (int) L ); // reserve only enough elements for this job

	while((signed) grid.size() < N) {
		for(int i = 0; i < (int) numIter; i++) { // loop over permutations

			curSet.clear(); // prep vector for next iteration
			v_sum = 0;

			for(int j = 0; j < (int) L; j++) { // loop over length of binary string
				if(i & (int) pow(2, (double) j) ) curSet.push_back(1);
				else curSet.push_back(0);
			}

			// get sum of current perumtation
			for(int i = 0; i< (signed) curSet.size(); i++) v_sum += curSet[i];

			if(v_sum == k) grid.insert(curSet);
		}
	}

	return grid;
}


// Function returns true if the given string contains the given character
bool hasChar(string txt, string ch) {
	size_t found;
	bool ret = false;

	found = txt.find(ch);
	if(found != string::npos) { ret = true; }

	return ret;
}



// Function parses the given file recording the PSM identifiers in it
void parsePSMfile(string srcFile) {

	g_PSMscoreSet.clear();

	ifstream inF;
	inF.open(srcFile.c_str(), ios::in);
	if(!inF.is_open()) {
		cerr << "\nERROR: PSM file '" << srcFile << "' not found\n";
		exit(-1);
	}

	string line;
	while(inF.good()) {
		line = ""; // prep for next iteration
		getline(inF, line);
		if(line.at(0) == '#') continue; // skip comment lines

		g_PSMscoreSet.insert(line);
	}
	inF.close();
}




// Function to reverse a given string
string reverseString(string input) {
	string ret;
	string::reverse_iterator c;
	for(c = input.rbegin(); c < input.rend(); c++) ret += *c;
	return ret;
}



// Function checks to see if a given double is infinite or not
bool isInfinite(double pV) {
	bool ret = false;
	if( fabs(pV) == numeric_limits<double>::infinity() ) ret = true;
	return ret;

}


// Function prints a progress bar to stderr
void printProgress(string txt, int ctr, int N) {
	double p = round_dbl((((double) ctr / (double) N) * 100.0), 0);

	cerr << txt << "  " << p << "% done\r";
}



// Function that returns true if a given file name exists
bool fileExists(const string& filename) {
   fstream file;
   bool ret = false;
   file.open(filename.c_str(), ios_base::out | ios_base::in); // will not create the file
   if(file.is_open()) ret = true;

   file.close();
   return ret;
}

