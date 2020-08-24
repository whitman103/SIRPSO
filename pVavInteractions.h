#ifndef PVAV_H
#define PVAV_H

#include <iostream>
#include <vector>
#include <string>
#include <tuple>

#include "SIR_MPI.h"



	
void pVav_RungeKutta(Particle* inParticle, vector<double>& speciesVec, double currentTime, double stoppingTime, double deltaT);

double dSyk(Particle* inParticle, vector<double>& currentSpecies);
double dVav(Particle* inParticle, vector<double>& currentSpecies);
double dSV(Particle* inParticle, vector<double>& currentSpecies);
double dpVav(Particle* inParticle, vector<double>& currentSpecies);
double dSHP1(Particle* inParticle, vector<double>& currentSpecies);
double dSHP1Vav(Particle* inParticle, vector<double>& currentSpecies);

double fitnessFunction_pVav(vector<vector<vector<int> > >& trueDistribution, vector<vector<vector<int> > >& testDistribution, string fitnessStyle);

vector<vector<vector<int> > > pVav_Gillespie(Particle* inParticle, Gillespie* inReactionObject, vector<double>& reportTimes, vector<int>& specNum, int numOfRuns);

void pVav_ParticleGillespieConvert();



#endif