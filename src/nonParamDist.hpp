/*
 * nonParamDist.hpp
 *
 *  Created on: Nov 14, 2012
 *      Author: dfermin
 */

#ifndef NONPARAMDIST_HPP_
#define NONPARAMDIST_HPP_

#include <deque>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include "structs.hpp";
#include "globals.hpp";

using namespace std;



double getMean2(deque<double> *x);
double getVariance2(deque<double> *ptr);

double normalDensity(double curTickMark, double curScore, double h);
void percentileTrim(deque<double> *ptr, double percent);


/* Intensity */
double getLogNPdensityInt_b(double x, modelParamStruct *paramPtr);
double getLogNPdensityInt_y(double x, modelParamStruct *paramPtr);
double getLogNPdensityInt_U(double x, modelParamStruct *paramPtr);

void estimateNonparamInt_b(list<double> *X, modelParamStruct *paramPtr, double mz_err);
void estimateNonparamInt_y(list<double> *X, modelParamStruct *paramPtr, double mz_err);
void estimateNonparamInt_U(list<double> *X, modelParamStruct *paramPtr, double mz_err);


/* Distance */
double getLogNPdensityDist(double x, modelParamStruct *paramPtr);
double getLogNPdensityDist_U(double x, modelParamStruct *paramPtr);

void estimateNonparamDist(list<double> *X, modelParamStruct *paramPtr, double mz_err);
void estimateNonparamDist_U(list<double> *X, modelParamStruct *paramPtr, double mz_err);


double mz_adjust_dist(double x);
void printDensity(modelParamStruct *paramPtr);

#endif /* NONPARAMDIST_HPP_ */


