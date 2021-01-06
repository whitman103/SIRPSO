#include <Rcpp.h>
using namespace Rcpp;

/* Requirements: Rtools and Rcpp. Rcpp can be downloaded from CRAN on R. Rtools is offered on CRAN and is easily found.

*/

/* Functions:

dSyk (species, parameters): return dSyk/dt term of ODE
dVav (species, parameters): return dVav/dt term of ODE
dSV (species, paramaters): return dSykVav/dt term of ODE
dpVav (species, paramaeters): return dpVav/dt term of ODE
dSHP1 (species, parameters): return dSHP/dt term of ODE
dSHP1Vav (species, parmaeters): return dSHP1vav term of ODE
checkInputs (species, paramaters): checks to make sure inputs are of right size and will be accepted by program
rungeKutta_evolve (species, parameters, startingTime, stoppingTime, deltaT): advances system from startingTime to stoppingTime using a time increment of deltaT. The vector "species" is modified in place.

*/

/* Species are:
1 Syk 
2 Vav 
3 Syk-Vav 
4 pVav
5 Shp1
6 Shp1-Vav

Parameters are:
1 k0
2 k1
3 k2
4 k3
5 k4
6 k6

Example usage: 
library(Rcpp)
sourceCpp(r_RungeKutta.cpp)
inputSpecies=c(600,200,0,0,60,0)
inputParameters=c(0.003,0.01,0.1,0.0013904,0.2,0.07)
checkInputs(inputSpecies,inputParameters)
startingTime=0
stoppingTime=3
deltaT=0.01
rungeKutta_evolve(inputSpecies,inputParameters,startingTime,stoppingTime,deltaT)

*/

/* The species input vector is modified in place. */

// [[Rcpp::export]]
double dSyk(NumericVector species, NumericVector parameters){
	double complexForm(-1.*parameters[0]*species[0]*species[1]);
	double complexBreak(parameters[1]*species[2]);
	double phosphor(parameters[2]*species[2]);
	return (complexForm+complexBreak+phosphor);
}

// [[Rcpp::export]]
double dVav(NumericVector species, NumericVector parameters){
	double complexForm(-1.*parameters[0]*species[0]*species[1]);
	double complexBreak(parameters[1]*species[2]);
	double shpBreak(parameters[5]*species[5]);
	return (complexForm+complexBreak+shpBreak);
}

// [[Rcpp::export]]
double dSV(NumericVector species, NumericVector parameters){
	double complexForm(parameters[0]*species[0]*species[1]);
	double complexBreak(-1.*parameters[1]*species[2]);
	double phosphor(-1.*parameters[2]*species[2]);
	return (complexForm+complexBreak+phosphor);
}

// [[Rcpp::export]]
double dpVav(NumericVector species, NumericVector parameters){
	double phosphor(parameters[2]*species[2]);
	double complexForm(-1.*parameters[3]*species[3]*species[4]);
	double complexBreak(parameters[4]*species[5]);
	return (phosphor+complexForm+complexBreak);
}

// [[Rcpp::export]]
double dSHP1(NumericVector species, NumericVector parameters){
	double complexForm(-1.*parameters[3]*species[3]*species[4]);
	double complexBreak(parameters[4]*species[5]);
	double phosphor(parameters[5]*species[5]);
	return (complexForm+complexBreak+phosphor);
}

// [[Rcpp::export]]
double dShp1Vav(NumericVector species, NumericVector parameters){
	double complexForm(parameters[3]*species[3]*species[4]);
	double complexBreak(-1.*parameters[4]*species[5]);
	double phosphor(-1.*parameters[5]*species[5]);
	return (complexForm+complexBreak+phosphor);
}

// [[Rcpp::export]]
int checkInputs(NumericVector species, NumericVector parameters){
	if (species.size()!=6){
		Rcout<<"Incorrect species vector size for Syk-Shp problem"<<"\n";
		Rcout<<"Needs an array of size 6."<<"\n";
		return 0;
	}
	if(parameters.size()!=6){
		Rcout<<"Incorrect parameter vector size for Syk-Shp problem."<<"\n";
		Rcout<<"Needs an array of size 6."<<"\n";
		return 0;
	}
	return 1;
}

// [[Rcpp::export]]
int rungeKutta_evolve(NumericVector species, NumericVector parameters, double currentTime, double stoppingTime, double deltaT){
	int numSpecies=species.size();
	int numOfTimeEvolutions((stoppingTime-currentTime)/deltaT);

	for(int tIndex=0;tIndex<numOfTimeEvolutions;tIndex++){
		NumericVector interSpecies(species.size());
		NumericVector k1(species.size());
		for(int i=0;i<(int)k1.size();i++){
			switch(i){
				case 0:
				k1[i]=deltaT*dSyk(species,parameters);
				break;
				case 1:
				k1[i]=deltaT*dVav(species,parameters);
				break;
				case 2:
				k1[i]=deltaT*dSV(species,parameters);
				break;
				case 3:
				k1[i]=deltaT*dpVav(species,parameters);
				break;
				case 4:
				k1[i]=deltaT*dSHP1(species,parameters);
				break;
				case 5:
				k1[i]=deltaT*dShp1Vav(species,parameters);
				break;
				default:
				break;
			}
			interSpecies[i]=species[i]+k1[i]/2.;
		}
		NumericVector k2(species.size());
		for(int i=0;i<(int)k2.size();i++){
			switch(i){
				case 0:
				k2[i]=deltaT*dSyk(species,parameters);
				break;
				case 1:
				k2[i]=deltaT*dVav(species,parameters);
				break;
				case 2:
				k2[i]=deltaT*dSV(species,parameters);
				break;
				case 3:
				k2[i]=deltaT*dpVav(species,parameters);
				break;
				case 4:
				k2[i]=deltaT*dSHP1(species,parameters);
				break;
				case 5:
				k2[i]=deltaT*dShp1Vav(species,parameters);
				break;
				default:
				break;
			}
			interSpecies[i]=species[i]+k2[i]/2.;
		}
		NumericVector k3(species.size());
		for(int i=0;i<(int)k3.size();i++){
			switch(i){
				case 0:
				k3[i]=deltaT*dSyk(species,parameters);
				break;
				case 1:
				k3[i]=deltaT*dVav(species,parameters);
				break;
				case 2:
				k3[i]=deltaT*dSV(species,parameters);
				break;
				case 3:
				k3[i]=deltaT*dpVav(species,parameters);
				break;
				case 4:
				k3[i]=deltaT*dSHP1(species,parameters);
				break;
				case 5:
				k3[i]=deltaT*dShp1Vav(species,parameters);
				break;
				default:
				break;
			}
			interSpecies[i]=species[i]+k3[i];
		}
		NumericVector k4(species.size());
		for(int i=0;i<(int)k4.size();i++){
			switch(i){
				case 0:
				k4[i]=deltaT*dSyk(species,parameters);
				break;
				case 1:
				k4[i]=deltaT*dVav(species,parameters);
				break;
				case 2:
				k4[i]=deltaT*dSV(species,parameters);
				break;
				case 3:
				k4[i]=deltaT*dpVav(species,parameters);
				break;
				case 4:
				k4[i]=deltaT*dSHP1(species,parameters);
				break;
				case 5:
				k4[i]=deltaT*dShp1Vav(species,parameters);
				break;
				default:
				break;
			}
		}

		for(int i=0;i<numSpecies;i++){
			species[i]+=(k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.);
		}
	}
	return 0;
}
