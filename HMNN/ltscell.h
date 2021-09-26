//
//  ltscell.h
//  cortex
//
//  Created by Diogo Porfirio de Castro Vieira on 27/11/11.
//  Copyright 2011 Universidade de Sao Paulo. All rights reserved.
//

#ifndef cortex_ltscell_h
#define cortex_ltscell_h

#include "neuron.h"

class ltscell: public neuron {
public:
    ltscell();
    ltscell(double aArg, double bArg, double cArg, double dArg);
    void setpar(double aArg, double bArg, double cArg, double dArg);        
    void setvrest(double vrest);
    void setvtresh(double vtresh);
    void setk(double k_);
    void setcap(double cap_);     
    
protected:
    void checkPeak(float time);
    void fx(double inj,float time);
    double a, b, c, d;
    double cap, k, vr, vt;
};
#endif
