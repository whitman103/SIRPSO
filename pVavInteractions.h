#ifndef PVAV_H
#define PVAV_H

#include <iostream>
#include <vector>

#include "SIR_MPI.h"



	
void pVav_RungeKutta(Particle* inParticle, vector<double>& speciesVec, double currentTime, double stoppingTime, double deltaT);

double dSyk(Particle* inParticle, vector<double>& currentSpecies);
double dVav(Particle* inParticle, vector<double>& currentSpecies);
double dSV(Particle* inParticle, vector<double>& currentSpecies);
double dpVav(Particle* inParticle, vector<double>& currentSpecies);
double dSHP1(Particle* inParticle, vector<double>& currentSpecies);
double dSHP1Vav(Particle* inParticle, vector<double>& currentSpecies);



#endif