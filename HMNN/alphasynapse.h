/* 
 * File:   alphasynapse.h
 * Author: diogopcv
 *
 * Created on 4 de Fevereiro de 2010, 13:24
 */

#ifndef _ALPHASYNAPSE_H
#define	_ALPHASYNAPSE_H

#include <vector>
#include <math.h>
#include "neuron.h"

class neuron;
using namespace std;

class alphasynapse {
public:
	alphasynapse();
	alphasynapse(double tauArg, double gmaxArg, double delayArg, double EsynArg);
	~alphasynapse();
	void setpar(double tauArg, double gmaxArg, double delayArg, double EsynArg);
	void setShortTerm(int argShort, double ap, double ataux);
	void evaluate(double time);
	void addevent(double spk);
	double getE();
	double getGsyn();
	void seth(double hArg);   
	void setK(double k_, int pos);
	void calcK(int pos);
	double getGaux();
protected:    	
	void startKs();
	vector <double> * spikes;       
	double * k;
	double tau;
	double delay;
	double gmax;
	double gsyn;
	double gsynaux;
	double Esyn;
	double h;   
};

#endif	/* _ALPHASYNAPSE_H */

