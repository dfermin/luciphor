/*
 * globas.hpp
 *
 *  Created on: Mar 18, 2011
 *      Author: dfermin
 */

#ifndef GLOBAS_HPP_
#define GLOBAS_HPP_

#define DEBUG

#include <vector>
#include <string>
#include <map>
#include <set>
#include <list>

using namespace std;


/****************
 * Global variables go here
 ****************/
extern string g_srcXMLfile; // name of pepXML file to parse by regex
extern string g_srcDir; // holds name of folder with the spectral files
extern string g_ext; // holds spectrum file format (mzXML or mzML, etc..)
extern string cmdLineArgs; // holds user-given command line arguments
extern string g_outputName;  // holds users' chosen output tag
extern string g_PSMfile; // holds list of PSMs that user wants scored
extern bool g_writeDTA; // if true, user wants to write spectra to disk
extern bool g_userDefinedOutput; // true means the user gave a name for the output file
extern double g_prob_threshold; // min. probability a peptide must have to be parsed
extern double g_model_prob; // min. prob. a peptide must have to be used for modeling
extern double g_MZ_ERR; // // fragment ion mass tolerance
extern double g_MIN_I; // min. intensity a peak needs to be considered for matching
extern double g_MIN_DIST; // min. distance (in log scale) tolerated for unmatched peaks
extern map<char, double> AAmass; // hold masses for amino acids
extern bool g_DEBUG_MODE; // true means the program is running in debug mode
extern bool g_FULL_MONTY; // true means the program should score and report all spectra
extern bool g_IS_HCD; // true means the data being processed is HCD data
extern bool g_NO_NL_PEAKS; // true means that neutral loss fragment ions will not be considered at all
extern bool g_NO_NL_PEAKS_MODEL; // true means that neutral loss fragment ions will not be used for modeling
extern bool g_NL_MODEL_ONLY; // true means that neutral loss peaks will only be considered for building the model
extern bool g_runAscoreInstead; // true means the user wants to run the Ascore algorithm on the data
extern bool g_randDecoyAA;
extern bool g_scoreSelect;
extern bool g_singleLetter;
extern bool g_LIMIT_CHARGE_STATE;
extern bool g_WRITE_TOP_TWO;
extern bool g_captureChargeStateModel; // in case you can't model a charge state, use the next closest charge states' parameters
extern bool g_useOnlySiteDetermIons; // when true, only site determining ions will be used in scoring a PSM (Luciphor only option)
extern bool g_usePPM;
extern bool g_removePrecursorNL;  // true means we remove MH-H2O and MH-H3PO4 peaks from spectrum
extern string g_BUILD_TIME;

extern string g_scoringMethod; // used to manage scoring metric
extern double g_dist_adj;
extern double g_NUM_PERMS_LIMIT;
extern double MIN_MZ;   // lowest m/z value we'll consider
extern int g_NUM_THREADS;
extern int g_intensityType;
extern int g_progressCtr;
extern int g_totalNumPSM;
extern int g_CHARGE_STATE;

extern int g_MIN_MODEL_NUM; // minimum number of PSMs you need for modeling
extern map<char, string> modAAmap;
extern map<char, char> decoyAA;
extern set<string> g_PSMscoreSet; //set of PSMs the user wants scored


// Global constants
const double TINY_NUM = 1e-10; // represents a tiny number
const double TINY_DEN = 1.0;
const double BIG_NUM = 1e10; // represents a big number
const double PPM = 1e-6; // parts per million
const int g_MIN_NUM_PSM_FOR_FLR = 50; // min. number of PSMs needed for FLR calculation

const double H = 1.00728;
const double H2O = 18.00947;
const double NH3 = 17.026549;
const double H3PO4 = 79.96633 + H2O;
const double OH = 17.00219;
const double e = 0.00054858026;


/*
 * Global function declarations
 */
void print_usage();
void parse_command_line_args(int argc, char *argv[]);
void initialize_AA_masses();
void addAAmass(string AA, double mass, char isVar);
void printProgress(string txt, int ctr, int N);
void parse_alternative_scoring(string *scoringStr, string *modelingStr);
void parsePSMfile(string srcFile);
void getSiteDeterminingIons(string pep1, string pep2, set<string> &ions1, set<string> &ions2);


double round_dbl(double r, int places);
double str2dbl(string ch);
double getIonMass(string seq, char ionType);

vector<string> split_string(string line);

string upperCaseSTY(string src);
string uc(string src);
string int2string(int i);
string dbl2string(double d);
string getTimeStamp();
string genRandDecoyPeptide(string srcSeq, int numSites);
string getExecutionParameters();
string repModAAchar(string *srcPtr, double nterm_mass);
string reverseString(string input);

double Factorial(double x); // used to compute factorials (n!)
double combinatorial(double n, double k);  //used to compute combinatorials (n choose k)

bool dbl_isnan(double var);
bool isInfinite(double pV);
bool hasChar(string txt, string ch);
bool containsPhosphoSite(string txt);
bool containsSTY(string txt);
bool isDecoyPep(string *srcPtr);
bool allValidAAs(string *srcPtr);
bool fileExists(const string& filename);

set< vector<int> > generateBinaryGrid(double L, double k);

#endif /* GLOBAS_HPP_ */
