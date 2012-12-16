/*
 * PepXMLClass.hpp
 *
 *  Created on: Mar 18, 2011
 *      Author: dfermin
 */

#ifndef PEPXMLCLASS_HPP_
#define PEPXMLCLASS_HPP_


#include <string>
#include <map>
#include <set>
#include <deque>
#include <fstream>
#include <iostream>
#include "globals.hpp"
#include "PSMClass.hpp"
#include "structs.hpp"
#include "nonParamDist.hpp"

// for reading spectra files use pwiz library
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "pwiz/data/msdata/MSData.hpp"

#include "boost/thread.hpp"

using namespace std;
using namespace pwiz::msdata;
using namespace pwiz::cv;
using namespace boost;


class PepXMLClass;


class PepXMLClass {
private:
	string fileName;
	ofstream logF; // used to write out summary of execution
	int numPSM;
	int numPSMs_forModeling;
	deque< PSMClass > *PSMvec;
	map<string, deque<scanIdStruct> > scanMap;
	map<string, repPSMStruct> repPSM_map;
	map<string, int> repPSMs;
	map<int, double> mzErrMap; // k = charge, v = calculated mz_err tolerance



public:
	void parsePepXMLfile();
	void readInSpectra();
	void quantileNormalize();
	void medianNormalize();
	void scoreSpectra();
	void removeComplexPSMs();
	void acquireModelParameters();
	void acquireModelParameters_HCD();
	void writeLuciphorResults();
	void openLogFile();
	void closeLogFile(string timeTxt);
	SpecStruct getSpectrum(scanIdStruct curScan, SpectrumListPtr allSpectra, CVID &nativeIdFormat_);
	int getMaxChargeState();
	void consolidatePSMs();
	void process_with_Ascore();
	void prunePSMdeq();
	void calcFLR();
	void adjustMZerr();
	void resetMZerr();
};

#endif /* PEPXMLCLASS_HPP_ */
