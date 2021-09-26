#include "simulation.h"

simulation::simulation () {
	tmax = 5000.0;
	h = 0.01;  
	setSeed(-184503872);
	nneuron = 1024;
	m = 0;	
	prob = 0.01;
	probEx = 0.9;
	gin = 0.1;
	gexc = 0.1;
	excClass = 0; 
	inbClass = 0;
	listn = new vector <neuron*>;
	lists = new vector <alphasynapse*>;		
}

simulation::~simulation () {
	neuron * ptrn = NULL;
	alphasynapse * ptrs = NULL;
	for(int i = 0; i < listn->size(); i++){
		ptrn = listn->at(i);
		delete ptrn;
		ptrn = NULL;
	}
	listn->clear();
	delete listn;
	for(int i = 0; i < lists->size(); i++){
		ptrs = lists->at(i);
		delete ptrs;
		ptrs = NULL;		
	}
	lists->clear();
	delete lists;
}

void simulation::run () {
	
	int nMax = (int) ((double)tmax/h) + 1;
	double t = 0.0, stim;	
	avgVolt = new double[nMax];
	lfp = new double[nMax];
	bool reset = true;
	
	for (int k = 0; k < nMax; k++){
		avgVolt[k] = 0.0;
		lfp[k] = 0.0;		
	}
	
//	FILE * pfile = fopen("voltage.dat","w");
	
	for (int k = 0; k < nMax; k++) {
        t += h;
        if (t >= 0.0 && t <= 200){
        	stim = 350.0;            
        }else{
        	if (reset){
        		reset = false;
        		for (int i = 0; i < listn->size(); i++) {
        			(listn->at(i))->clearEvents();
        		}
        	}
        	stim = 0.0;
        }
        for (int i = 0; i < listn->size(); i++) {
	        (listn->at(i))->evaluate(stim * randNext(), t);
	        avgVolt[k] += (listn->at(i))->getW(0);
	        lfp[k] += (listn->at(i))->getIsyn();	        
        }
        avgVolt[k] = avgVolt[k]/listn->size();
        lfp[k] = lfp[k]/listn->size();
//        fprintf(pfile, "%lf\t%lf\t%lf\t%lf\n", t, (listn->at(0))->getW(0), (listn->at(1))->getW(0), (listn->at(251))->getW(0));
    }
    
//    fclose(pfile);   	
}

void simulation::printLfp(const char * pchar){
	int nMax = (int) ((double)tmax/h) + 1;
	FILE * pfile = fopen(pchar,"w");
	for (int k = 0; k < nMax; k++){
		fprintf(pfile,"%lf\n",lfp[k]);	
	}	
	fclose(pfile);
}

void simulation::printAvgVolt(const char * pchar){
	int nMax = (int) ((double)tmax/h) + 1;
	FILE * pfile = fopen(pchar,"w");
	for (int k = 0; k < nMax; k++){
		fprintf(pfile,"%lf\n",avgVolt[k]);	
	}	
	fclose(pfile);	
}	

void simulation::setSeed (long seed) {
	this->idum = seed;
}

//void simulation::setGenRan (gsl_rng * r) {
//	this->r = r;
//}

float simulation::randNext() {
	int j; 
	long k; 
	static long idum2=123456789; 
	static long iy=0; 
	static long iv[NTAB]; 
	float temp;
	if (idum <= 0) { 
		if (-(idum) < 1) idum=1; 
		else idum = -(idum); 
		idum2=(idum); 
		for (j=NTAB+7;j>=0;j--) {
			k=(idum)/IQ1; 
			idum=IA1*(idum-k*IQ1)-k*IR1; 
			if (idum < 0) idum += IM1; 
			if (j < NTAB) iv[j] = idum;
		} 
		iy=iv[0];
	} 
	k=(idum)/IQ1; 
	idum=IA1*(idum-k*IQ1)-k*IR1; 
	if (idum < 0) idum += IM1; 
	k=idum2/IQ2; 
	idum2=IA2*(idum2-k*IQ2)-k*IR2; 
	if (idum2 < 0) idum2 += IM2; 
	j=iy/NDIV; iy=iv[j]-idum2; iv[j] = idum; 
	if (iy < 1) iy += IMM1; 
	if ((temp=AM*iy) > RNMX) return RNMX; 
	else return temp;
}  

void simulation::setGin(double gin){
	this->gin = gin;
}

void simulation::setGexc(double gexc){
	this->gexc = gexc;	
}

void simulation::setH(double h){
	this->h = h;	
}

void simulation::setTmax(double tmax){
	this->tmax = tmax;	
}

void simulation::setNNeuron(int num){
	this->nneuron = num;	
}

void simulation::setM(int m){
	this->m = m;	
}

void simulation::setProbConn(double prob){
	this->prob = prob;	
}

void simulation::setProbExc(double probEx){
	this->probEx = probEx;	
}

void simulation::setExcClass(short int exc){
	this->excClass = exc;
}

void simulation::setInbClass(short int inb){
	this->inbClass = inb;
}

/*
RS = 0
IBS = 1
CHS = 2
FS = 0
LTS = 1 
*/

void simulation::createNet () {
    
    long idum = -1089532789; 
    alphasynapse * syn;
    int i, j, count, sizem = nneuron, num, mod1, mod2, n1, n2, nn;
    double sort;
    int ** listsyn = new int *[1000000];
	for (i=0; i<1000000; i++)
		listsyn[i] = new int[3];		
    
    neuron * cell;
    double rnd;
    for (i = 0; i < nneuron; i++) {      
        num = (int) (randNext()*5);
        if (num == 4){
        	if (inbClass == 0)
        		cell = new fscell();
        	else
				cell = new ltscell();
        	cell->setInhibitory();
        }
        else{
        	cell = new izhiCom();
        	if (excClass == 0)
				((izhiCom *)cell)->turnRS();
			else if (excClass == 1)
				((izhiCom *)cell)->turnIBS();
			else
				((izhiCom *)cell)->turnCHS();
       		cell->setExcitatory();
        }
        cell->seth(h);
        listn->push_back(cell);          
    }
    
    count = 0;
	for (i = 0; i < nneuron; i++) {  
		for (j = 0; j < nneuron; j++) {
			sort = randNext();			
			if (j!=i & sort < prob) {
				listsyn[count][0] = i;
				listsyn[count][1] = j;	
				listsyn[count][2] = 1;		
				count++;
			}
		}	
	}
	listsyn[count][0] = -1;	
	
    for (i = 1; i <= m; i++){
    	sizem = sizem/2;
    	count = 0;
		while(listsyn[count][0] != -1){  									
			n1 = 1; n2 = 1; mod1 = 0; mod2 = 0;			
			while(listsyn[count][1] >= mod2 + sizem){
				mod2 += sizem;
				n2++;	
			}												
			while(listsyn[count][0] >= mod1 + sizem){
				mod1 += sizem;
				n1++;	
			}				
			if (n1!=n2 & listn->at(listsyn[count][0])->getTypeSyn() == 0){	
				nn = mod1 + (int) (randNext()*sizem);
				listsyn[count][1] = nn;				
			}
			if (n1!=n2 &  listn->at(listsyn[count][0])->getTypeSyn() == 1){	
				sort = randNext();
				if (sort < probEx) {
					nn = mod1 + (int) (randNext()*sizem);
					listsyn[count][1] = nn;	
					listsyn[count][2] = 1;
				}
				else {
					listsyn[count][2]++;
				}											
			}							  	
	    	count++;
    	}
    }
    
//    FILE * pfile = fopen("mtx.dat", "w");
    
    count = 0;
    int numexc = 0, numinib = 0;
    while(listsyn[count][0]!=-1){
    	syn = new alphasynapse;
    	syn->seth(h);	
    	lists->push_back(syn);
    	if( listn->at(listsyn[count][0])->getTypeSyn()==0){
    		syn->setpar(6.0, gin, 1.0, -80.0);
//    		fprintf(pfile,"%d\t%d\t%d\n",listsyn[count][0], listsyn[count][1], listn->at(listsyn[count][0])->getTypeSyn());
    	}
    	else{
	    	syn->setpar(5.0, gexc, 1.0, 0.0);
//    		fprintf(pfile,"%d\t%d\t%d\n",listsyn[count][0], listsyn[count][1], listn->at(listsyn[count][0])->getTypeSyn());	    	
    	}	    	
	    listn->at(listsyn[count][0])->makeconnection( listn->at(listsyn[count][1]), syn);		    	
    	count++;
    }
    
//    fclose(pfile);
    
	for (i=0; i<1000000; i++)
		delete[] listsyn[i];
}   
   
void simulation::rasterdata(char * pchar) {
    ofstream fout, fplot;
    fout.open(pchar);
    int i, j, k, max_size = 0;
    
    for (i = 0; i < listn->size(); i++)
        if (!listn->at(i)->getevents()->empty())
        	if (listn->at(i)->getevents()->size() > max_size)
        		max_size = listn->at(i)->getevents()->size();    
	
    for (i = 0; i < listn->size(); i++) {
        if (!listn->at(i)->getevents()->empty()) {
        	fout << i << "\t";
            for (j = 0; j < listn->at(i)->getevents()->size(); j++) {
                fout << listn->at(i)->getevents()->at(j) << "\t";
            }
            if (listn->at(i)->getevents()->size() < max_size){
            	k = max_size - listn->at(i)->getevents()->size();
            	for (j = 0;j < k; j++)
            		fout << 0.0 << "\t";	
            }
            fout << "\n";
        }
    }
    
    /*fout.close();
    fplot.open("plot");
    fplot << "set term X11;\nreset;\nset nokey;\nset xrange [0:5000];\nset yrange [0:1024];\nplot \"" << pchar << "\" using 2:1 lt -1 pt 0";
    for (i = 2; i <= max_size; i++){
        fplot << ", \"" << pchar << "\" using " << i << ":1 lt -1 pt 0";
    }
	fplot << "\nset term postscript enhanced color;\nset output \"" << pchar << ".eps\";\nreplot;";    
    fplot.close();
    system("gnuplot plot");*/  
}

void simulation::printTime(string name){
	FILE * pfile;
	pfile = fopen(name.c_str(), "a");
	neuron * nrn;
	double tmp = 0.0, tmp2 = 0.0, aux = 0.0, last = 0.0, std = 0.0;
	for (int i = 0; i < nneuron; i++){
		nrn = listn->at(i);
		if ((nrn->getevents())->empty()) continue;
		aux = (nrn->getevents())->back();
		tmp += aux;
		tmp2 += aux*aux;
		if (aux > last) last = aux;
	}
	tmp = tmp/(listn->size());
	tmp2 = tmp2/(listn->size());
	std = sqrt(tmp2 - tmp*tmp);
	fprintf(pfile,"%lf\t%lf\t%lf\t%lf\t%lf\n",gin, gexc, last, tmp, std);
	fclose(pfile);
}  	

void simulation::printAvgFreq(string name){
	FILE * pfile;
	pfile = fopen(name.c_str(), "a");
	neuron * nrn;
	double tmp = 0.0, tmp2 = 0.0, aux = 0.0, std = 0.0;
	int count = 0;
	for (int i = 0; i < nneuron; i++){
		nrn = listn->at(i);
		if ((nrn->getevents())->size() < 2) continue;
		for (int j = 1; j < (nrn->getevents())->size(); j++){
			aux = 1000/((nrn->getevents())->at(j) - (nrn->getevents())->at(j-1));
			tmp += aux;
			tmp2 += aux*aux;
			count++;		
		}
	}
	tmp = tmp/count;
	tmp2 = tmp2/count;
	std = sqrt(tmp2 - tmp*tmp);
	fprintf(pfile,"%lf\t%lf\t%lf\t%lf\n",gin, gexc, tmp, std);
	fclose(pfile);
}  
