/*
 *  neuron.h
 *  
 *
 *  Created by Diogo Porfirio de Castro Vieira on 02/12/09.
 *  Copyright 2009 Universidade de Sao Paulo. All rights reserved.
 *
 */
#ifndef _NEURON_
#define _NEURON_
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <complex>
#include <vector>
#include "alphasynapse.h"

class alphasynapse;

class neuron {
public:
    neuron();
    virtual ~neuron();     
    void evaluate(double inj, float time);
    void setw0(double * w);
    double getW(int ind);
    void seth(double hArg);
    void makeconnection(neuron * dend, alphasynapse * syn);    
    void addsyndend(alphasynapse * syn);
    void setId(int myId);
    int getID();
    std::vector<double> * getevents();  
    void setExcitatory();
    void setInhibitory();
    int getTypeSyn();    
    double getIsyn();
    void clearEvents();
protected:
    virtual void checkPeak(float time) = 0;    	
    virtual void fx(double inj, float time) = 0;
	void startKs();   
    void sendevent(double time);
    double calcsyncurrent(float time) ;
    std::vector <alphasynapse*> * saxon;
    std::vector <alphasynapse*> * sdend;
    std::vector <double> * events;
    double ** ks;
    double * fbuf;
    double * w, * waux;
    int _ID;
    double h;    
    int numeq;
    int typesyn;
};

#endif
