/* 
 * File:   neuron.cpp
 * Author: diogopcv
 * 
 * Created on 4 de Fevereiro de 2010, 13:22
 */

#include "neuron.h"
using std::vector;

neuron::neuron() {
	_ID = 0;
    saxon = new vector <alphasynapse*>;
    sdend = new vector <alphasynapse*>;
    events = new vector <double>;
}    

neuron::~neuron(){
	delete[] fbuf;
	delete[] w;
	delete[] waux;
	for (int i = 0; i < numeq; i++)
		delete[] ks[i];	
	delete sdend;
	delete events; 	
	delete saxon; 	    	
}     

void neuron::setw0(double * w_) {
	for (int i = 0; i < numeq; i++) {
		w[i] = w_[i];
	}
}

void neuron::clearEvents(){
	events->clear();
	return;
}

double neuron::getW(int ind) {
    return w[ind];
}

void neuron::seth(double hArg) {
    h = hArg;
}

double neuron::calcsyncurrent(float time) {
    double isyn = 0.0;
    int i = 0;
    while (i < sdend->size()) {
        isyn += (sdend->at(i))->getGaux()*(waux[0] - (sdend->at(i))->getE());
        i++;
    }
	return isyn;
}    

double neuron::getIsyn(){
    double isyn = 0.0;
    int i = 0;
    while (i < sdend->size()) {
        isyn += (sdend->at(i))->getGsyn()*(w[0] - (sdend->at(i))->getE());
        i++;
    }
	return isyn;	
}

void neuron::evaluate(double inj, float time) {

	// Calculando k1 
	for (int i = 0; i < numeq; i++)
		waux[i] = w[i];
	for(int i = 0; i < sdend->size(); i++)	 
		(sdend->at(i))->calcK(0);			            
	fx(inj, time);       
	for (int i = 0; i < numeq; i++)
		ks[i][0] = h * fbuf[i];
	
	// Calculando k2
	for (int i = 0; i < numeq; i++)
		waux[i] = w[i] + ks[i][0]/2;
	for(int i = 0; i < sdend->size(); i++)	 
		(sdend->at(i))->calcK(1);	              
	fx(inj, time);
	for (int i = 0; i < numeq; i++)
		ks[i][1] = h * fbuf[i];
	
	// Calculando k3
	for (int i = 0; i < numeq; i++)
		waux[i] = w[i] + ks[i][1]/2;
	for(int i = 0; i < sdend->size(); i++)	 
		(sdend->at(i))->calcK(2);  	             
	fx(inj, time);
	for (int i = 0; i < numeq; i++)
		ks[i][2] = h * fbuf[i];
	
	// Calculando k4
	for (int i = 0; i < numeq; i++)
		waux[i] = w[i] + ks[i][2];
	for(int i = 0; i < sdend->size(); i++)	
		(sdend->at(i))->calcK(3); 	        
	fx(inj, time);
	for (int i = 0; i < numeq; i++)
		ks[i][3] = h * fbuf[i];
	
	// Calculando w's novos
	for (int i = 0; i < numeq; i++)
		w[i] += (ks[i][0] + 2*ks[i][1] + 2*ks[i][2] + ks[i][3])/6;                                              
	
	for(int i = 0; i < sdend->size(); i++) 
		(sdend->at(i))->evaluate(time);             
	
	checkPeak(time);
}

void neuron::addsyndend(alphasynapse * syn) {
	sdend->push_back(syn);
}

void neuron::makeconnection(neuron * dend, alphasynapse * syn) {
    saxon->push_back(syn);
    dend->addsyndend(syn);
}

void neuron::sendevent(double time) {
    int i = 0;
    events->push_back(time);
    while (i < saxon->size()) {
        (saxon->at(i))->addevent(time);
        i++;
    }
}	

vector<double> * neuron::getevents(){
    return events;
}

void neuron::startKs(){
	for(int i = 0; i < numeq; i++)
		for(int j = 0; j < 4; j++)
			ks[i][j] = 0.0; 		
}  	

void neuron::setId(int myId){
	_ID = myId;
}

int neuron::getID(){
	return _ID;
}

void neuron::setExcitatory(){
	typesyn = 1;
}

void neuron::setInhibitory(){
	typesyn = 0;
}

int neuron::getTypeSyn(){
	return typesyn;
}