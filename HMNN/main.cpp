/* 
 * File:   main.cpp
 * Author: diogopcv
 *
 * Created on 4 de Fevereiro de 2010, 13:19
 */
 
/*
RS = 0
IBS = 1
CHS = 2
FS = 0
LTS = 1 
*/

#include <iostream>
#include "simulation.h"
#include "mpi.h"
#include <math.h>
#define DIM 20

using namespace std;

void funcaoNo(){
	double dadoIn[2];
	int dadoOut;
	FILE * pfile;
	MPI::Status status;
	simulation * sim = NULL;
	ostringstream timeFile, freqFile;
	while (1) {
		MPI::COMM_WORLD.Recv(dadoIn, 2, MPI::DOUBLE, 0, MPI::ANY_TAG, status);
		if (dadoIn[0] < 0.0) {
			break;
		}
		double tmax = 2200, h = 0.01;
		int nMax =(int) ((double)tmax/h) + 1; 		 
		simulation * sim = NULL;
		timeFile << "CHSFS_M4_last" << dadoIn[0] << "_" << dadoIn[1] << ".dat";
		freqFile << "CHSFS_M4_freqAvg" << dadoIn[0] << "_" << dadoIn[1] << ".dat";		
		
		for(int i = 0; i < 25; i++) {
			long rnd = -184503872 + i*10000;
			sim = new simulation();
			sim->setM(4);
			sim->setExcClass(2); // CHS
			sim->setInbClass(0); // FS		
			sim->setH(h);
			sim->setTmax(tmax);
			sim->setSeed(rnd);
			sim->setGin(dadoIn[0]);
			sim->setGexc(dadoIn[1]);
			sim->createNet();	
			sim->run();
			sim->printTime(timeFile.str());
			sim->printAvgFreq(freqFile.str());
			delete sim;
			sim = NULL;			
		}
			
		dadoOut = 1;
		MPI::COMM_WORLD.Send(&dadoOut, 1, MPI::INT, 0, status.Get_tag());
		timeFile.str("");
		freqFile.str("");
	}
	return;	
}

int main(int argc, char * argv[]) {
  
  int numtasks, rank, numSims = DIM*DIM, tag=1;
  MPI::Status status;
 
  MPI::Init(argc,argv);
  numtasks = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();
  
  if(rank == 0){
    double dadoOut[2];
    int dadoIn;
    
    // Primeira leva ...
    int proximo_processo = 1, enviados = 0, recebidos = 0;
    while (enviados<numSims && proximo_processo<numtasks) {
      dadoOut[0] = ((double) ((enviados)/DIM)+1)*5.0;
      dadoOut[1] = ((double) ((enviados)%DIM)+1)*1.0;	
      MPI::COMM_WORLD.Send(dadoOut, 2, MPI::DOUBLE, proximo_processo, 0);
      enviados++;
      proximo_processo++;
    }
    MPI::COMM_WORLD.Recv(&dadoIn, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status); 
    recebidos++;  
    
    // O que faltou ...
    while (enviados<numSims || recebidos<numSims) {
     if (enviados<numSims) {
		dadoOut[0] = ((double) ((enviados)/DIM)+1)*5.0;
		dadoOut[1] = ((double) ((enviados)%DIM)+1)*1.0;	
        MPI::COMM_WORLD.Send(dadoOut, 2, MPI::DOUBLE, status.Get_source(), 0);
        enviados++;
      }
      if (recebidos<numSims) {
        MPI::COMM_WORLD.Recv(&dadoIn, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, status);
        recebidos++;
      }
    }
    
    // Finalizando ...
    dadoOut[0] = -1.0;
    dadoOut[1] = -1.0;
    for (int processo_atual=1; processo_atual<numtasks; processo_atual++) {
      MPI::COMM_WORLD.Send(dadoOut, 2, MPI::DOUBLE, processo_atual, 0);
    }
  }
  else{
    funcaoNo();	
  }
  
  MPI::Finalize();
  return 0;
}
