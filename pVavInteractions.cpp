#include <iostream>
#include <vector>

using namespace std;

#include "SIR_MPI.h"
#include "pVavInteractions.h"



void pVav_RungeKutta(Particle* currentParticle, vector<double>& speciesVec, double currentTime, double stoppingTime, double deltaT){
	int numSpecies=speciesVec.size();
	(*currentParticle).unwrap_pVavParameters();
    int n=(int)(((stoppingTime-currentTime))/(deltaT));

    vector<double (*)(Particle*,vector<double>&)> interactionPointer=(*currentParticle).interactionFunctions;

    for(int t=0;t<n;t++){
        vector<double> interSpecies(numSpecies,0);
        vector<double> k1(numSpecies,0);
        for(int i=0;i<(int)k1.size();i++){
            k1[i]=deltaT*interactionPointer[i](currentParticle,speciesVec);
			interSpecies[i]=speciesVec[i]+k1[i]/2.;
        }
        vector<double> k2(numSpecies,0);
        for(int i=0;i<(int)k2.size();i++){
            k2[i]=deltaT*interactionPointer[i](currentParticle,interSpecies);
			interSpecies[i]=speciesVec[i]+k2[i]/2.;
        }
        vector<double> k3(numSpecies,0);
        for(int i=0;i<(int)k3.size();i++){
            k3[i]=deltaT*interactionPointer[i](currentParticle,interSpecies);
			interSpecies[i]=speciesVec[i]+k3[i];
        }
        vector<double> k4(numSpecies,0);
        for(int i=0;i<(int)k4.size();i++){
            k4[i]=deltaT*interactionPointer[i](currentParticle,interSpecies);
        }
        for(int i=0;i<numSpecies;i++){
            speciesVec[i]+=(k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.);
        }
    }
}


//Species order is Syk, Vav, Syk-vav, pVav, SHP1, SHP1-Vav

double dSyk(Particle* inParticle, vector<double>& currentSpecies){
	double complexForm(-1.*(*inParticle).k0*currentSpecies[0]*currentSpecies[1]);
	double complexBreak((*inParticle).k1*currentSpecies[2]);
	double phosphor((*inParticle).k2*currentSpecies[2]);
	return (complexForm+complexBreak+phosphor);
}

double dVav(Particle* inParticle, vector<double>& currentSpecies){
	double complexForm(-1.*(*inParticle).k0*currentSpecies[0]*currentSpecies[1]);
	double complexBreak((*inParticle).k1*currentSpecies[2]);
	double ShpBreak((*inParticle).k5*currentSpecies[5]);
	return (complexForm+complexBreak+ShpBreak);
}
	
double dSV(Particle* inParticle, vector<double>& currentSpecies){
	double complexForm((*inParticle).k0*currentSpecies[0]*currentSpecies[1]);
	double complexBreak(-1.*(*inParticle).k1*currentSpecies[2]);
	double phosphor(-1.*(*inParticle).k2*currentSpecies[2]);
	return (complexForm+complexBreak+phosphor);
}
	
double dpVav(Particle* inParticle, vector<double>& currentSpecies){
	double phosphor((*inParticle).k2*currentSpecies[2]);
	double complexForm(-1.*(*inParticle).k3*currentSpecies[3]*currentSpecies[4]);
	double complexBreak((*inParticle).k4*currentSpecies[5]);
	return (phosphor+complexForm+complexBreak);
}
	
double dSHP1(Particle* inParticle, vector<double>& currentSpecies){
	double complexForm(-1.*(*inParticle).k3*currentSpecies[4]*currentSpecies[3]);
	double complexBreak((*inParticle).k4*currentSpecies[5]);
	double phosphor((*inParticle).k5*currentSpecies[5]);
	return (phosphor+complexBreak+complexForm);
}
	
double dSHP1Vav(Particle* inParticle, vector<double>& currentSpecies){
	double complexForm((*inParticle).k3*currentSpecies[3]*currentSpecies[4]);
	double complexBreak(-1.*(*inParticle).k4*currentSpecies[5]);
	double phosphor(-1.*(*inParticle).k5*currentSpecies[5]);
	return (complexForm+complexBreak+phosphor);
}
	
	
	