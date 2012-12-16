//============================================================================
// Name        : luciPHOr2.cpp
// Author      : Damian Fermin, Ph.D
// Version     :
// Copyright   : 2011
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <map>
#include <deque>
#include <boost/date_time/posix_time/posix_time.hpp> // to calculate run time
#include "globals.hpp"
#include "PepXMLClass.hpp"


using namespace std;



int main(int argc, char *argv[]) {

    char *build_date = __DATE__;
    char *build_time = __TIME__;

    g_BUILD_TIME.assign(build_date);
    g_BUILD_TIME.append("-");
    g_BUILD_TIME.append(build_time);

	cerr << "\nLuciphor BUILD: " << g_BUILD_TIME << endl << endl;

	parse_command_line_args(argc, argv);
	initialize_AA_masses();

	boost::posix_time::ptime start_time(boost::posix_time::microsec_clock::local_time());

	srand(time(NULL)); // initialize random number generator with time seed

	PepXMLClass *pepXMLfile = NULL;

	pepXMLfile = new PepXMLClass();

	pepXMLfile->parsePepXMLfile();
	pepXMLfile->removeComplexPSMs();

	pepXMLfile->openLogFile();

	pepXMLfile->readInSpectra();

	if(g_runAscoreInstead) pepXMLfile->process_with_Ascore();
	else {
		if(g_captureChargeStateModel == false) pepXMLfile->prunePSMdeq();

		if(g_IS_HCD) pepXMLfile->acquireModelParameters_HCD();
		else pepXMLfile->acquireModelParameters();

		pepXMLfile->scoreSpectra();

		if(g_randDecoyAA) pepXMLfile->calcFLR();

		pepXMLfile->writeLuciphorResults();
	}

	boost::posix_time::ptime end_time(boost::posix_time::microsec_clock::local_time());
	boost::posix_time::time_duration delta( end_time - start_time );

	pepXMLfile->closeLogFile( boost::posix_time::to_simple_string(delta) );

	delete(pepXMLfile); pepXMLfile = NULL;


	cerr << "\n\nDone!\n"
		 << "Total run time (HH:MM::SS.fffff):  "
		 << boost::posix_time::to_simple_string(delta) << endl << endl;

	return 0;
}
