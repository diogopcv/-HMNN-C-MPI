/* 
 * File:   alphasynapse.cpp
 * Author: diogopcv
 * 
 * Created on 4 de Fevereiro de 2010, 13:24
 */

#include "alphasynapse.h"

alphasynapse::alphasynapse() {
    gmax = 1.0;
    tau = 1.0;
    delay = 5.0;
    Esyn = 50.0;
    spikes = new vector <double>;
    gsyn = 0.0;      
    k = new double[4];
    startKs();  
    seth(0.1);        
}

alphasynapse::alphasynapse(double tauArg, double gmaxArg, double delayArg, double EsynArg) {
    gmax = gmaxArg;
    tau = tauArg;
    delay = delayArg;
    gsyn = 0.0;
    Esyn = EsynArg;
    spikes = new vector <double>;
    gsyn = 0.0;
    k = new double[4];
    startKs();
    seth(0.1);   
}  

alphasynapse::~alphasynapse(){
 	delete[] k;
 	delete spikes;
}

void alphasynapse::setpar(double tauArg, double gmaxArg, double delayArg, double EsynArg) {
    gmax = gmaxArg;
    tau = tauArg;
    delay = delayArg;
    Esyn = EsynArg;
}

void alphasynapse::addevent(double spk) {
spikes->push_back(spk);
}

void alphasynapse::evaluate(double time) {	 
	gsyn += (k[0] + 2*k[1] + 2*k[2] + k[3])/6;     	
	if (!spikes->empty()){
    	double s = time - spikes->at(0) - delay;          
    	if (s>=0){
    		gsyn += gmax;
    		spikes->erase(spikes->begin());
    	}
	}        	      		
}

void alphasynapse::seth(double hArg) {
    h = hArg;
}

double alphasynapse::getGsyn() {
    return gsyn;
}

double alphasynapse::getE() {
    return Esyn;
}

void alphasynapse::setK(double k_, int pos){
	k[pos] = k_;
}

void alphasynapse::startKs(){
	for(int i = 0; i < 4; i++)
		k[i] = 0.0; 	
}

void alphasynapse::calcK(int pos){
	if (pos == 0){
		gsynaux = gsyn;
		k[0] = - h * (gsynaux/tau);			
	}
	if (pos == 1){
		gsynaux = gsyn + k[0]/2;
		k[1] = - h * (gsynaux/tau);			
	}
	if (pos == 2){
		gsynaux = gsyn + k[1]/2;
		k[2] = - h * (gsynaux/tau);		
	}			
	if (pos == 3){
		gsynaux = gsyn + k[2];
		k[3] = - h * (gsynaux/tau);
	}												
}

double alphasynapse::getGaux(){
	return gsynaux;
}