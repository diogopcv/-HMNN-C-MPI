/*
 *  izhi.h
 *  
 *
 *  Created by Diogo Porfirio de Castro Vieira on 02/12/09.
 *  Copyright 2009 Universidade de Sao Paulo. All rights reserved.
 *
 */
#ifndef _IZHI_
#define _IZHI_
#include "neuron.h"

class izhiCom : public neuron {
   public:
    izhiCom();
    izhiCom(double aArg, double bArg, double cArg, double dArg);
    void setpar(double aArg, double bArg, double cArg, double dArg);        
    void setvrest(double vrest);
    void setvtresh(double vtresh);
    void setk(double k_);
    void setcap(double cap_);   
    void turnRS();
    void turnIBS();
    void turnCHS();  
   protected:
    void checkPeak(float time);
    void fx(double inj,float time);
    double a, b, c, d;
    double cap, k, vr, vt;
};
    
#endif
