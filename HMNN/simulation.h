#ifndef _SIMULATION_H
#define	_SIMULATION_H

#define IM1 2147483563 
#define IM2 2147483399 
#define AM (1.0/IM1) 
#define IMM1 (IM1-1) 
#define IA1 40014 
#define IA2 40692 
#define IQ1 53668 
#define IQ2 52774 
#define IR1 12211 
#define IR2 3791 
#define NTAB 32
#define NDIV (1+IMM1/NTAB) 
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include "fscell.h"
#include "ltscell.h"  
#include "izhiCom.h"
#include "alphasynapse.h"
//#include <gsl/gsl_rng.h>

using namespace std;

class simulation {
    public:
		simulation();
		~simulation();
		void run();
		void setSeed(long seed);
		void setGin(double gin);
		void setGexc(double gexc);	
		void setH(double h);
		void setTmax(double tmax);
		void setNNeuron(int num);
		void setM(int m);
		void setProbConn(double prob);
		void setProbExc(double probEx);
		void rasterdata(char * pchar);
		void printTime(string name);
		void printAvgFreq(string name);	
		void createNet();
		void printLfp(const char * pchar);
		void printAvgVolt(const char * pchar);
		void setExcClass(short int exc);
		void setInbClass(short int inb);			
//		void setGenRan (gsl_rng * r);
    protected:
    	float randNext();
//		gsl_rng * r;
    	vector <neuron*> * listn;
    	vector <alphasynapse*> * lists;
    	long idum;
    	int nneuron, m;
    	short int excClass, inbClass;
    	double h, tmax;
    	double gin, gexc;
    	double prob, probEx;
    	double * avgVolt, * lfp;
};
    
#endif
