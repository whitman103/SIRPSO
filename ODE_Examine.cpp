#include <iostream>
#include <fstream>
#include <vector>
#include <string>

	using namespace std;
#include <boost/random/mersenne_twister.hpp>
	boost::mt19937 generator;
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lognormal_distribution.hpp>

#include "SIR_MPI.h"
#include "pVavInteractions.h"

vector<tuple<double,double> > pVavSetBounds(vector<double>& inSpeciesCounts,vector<double>& inTimes);
int loadPvavInputs(vector<double>& speciesVector, vector<double>& initParameters, vector<double>&stoppingTimes, vector<tuple<double,double> >& bounds, string inFile);

int main(){


	//T, I, V, R
	const int numOfSpecies(6);
	const int momentVectorSize(2*numOfSpecies+(numOfSpecies-1)*numOfSpecies/2);
	
	vector<double> speciesVector={600,200,0,0,60,0};

	//k0, k1, k2, k3, k4, k5
	const int numOfParameters(6);
	vector<double> initParameters={0.008,0.1,1.0,0.013904,0.05,0.07};
	vector<tuple<double,double> > twoBounds;
	vector<double> stoppingTimes={0,100};
	int errorFlag=loadPvavInputs(speciesVector,initParameters,stoppingTimes,twoBounds, "0");
	twoBounds=pVavSetBounds(speciesVector,stoppingTimes);

	vector<double (*)(Particle*, vector<double>&)> interactionFuncts={dSyk,dVav,dSV,dpVav,dSHP1,dSHP1Vav};

	boost::normal_distribution<> standardNormal(0,1);
    boost::mt19937 generator2;
    generator2.seed(generator());
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> > normalGenerator(generator2,standardNormal);
	boost::mt19937 exGenerator;
	exGenerator.seed(generator2());

	//Generates cholesky matrix to produce lognormal distributions
	vector<vector<double> > inValues;
	vector<double> means;
	loadCovariance(means, inValues,"outCov");

	vector<boost::mt19937> exNoiseEngines(3);
	vector<boost::variate_generator<boost::mt19937, boost::lognormal_distribution<double> > > extrinsicNoiseGenerators;
	for(int i=0;i<(int)exNoiseEngines.size();i++){
		exNoiseEngines[i].seed(time(NULL));
		double scaledMean(means[i]);
		double scaledSD(inValues[i][i]);
		boost::lognormal_distribution<> currentTest(scaledMean,scaledSD);
		boost::variate_generator<boost::mt19937,boost::lognormal_distribution<double> > createdEngine(exNoiseEngines[i],currentTest);
		extrinsicNoiseGenerators.push_back(createdEngine);
	}

	Particle trueParticle=Particle(numOfParameters,initParameters,twoBounds,interactionFuncts,1.);

	rungeKuttaUpdate(&trueParticle,speciesVector,0,10,0.002);

	return 0;
}

vector<tuple<double,double> > pVavSetBounds(vector<double>& inSpeciesCounts,vector<double>& inTimes){

	double overallRateLimit(600); //10 per second
	vector<tuple<double,double> > outBounds(6);
	double maxTime(inTimes[inTimes.size()-1]);
	double k0Min(0), k0Max(0);
	k0Max=overallRateLimit/(inSpeciesCounts[0]*inSpeciesCounts[1]);
	k0Min=(1./maxTime)/(inSpeciesCounts[0]*inSpeciesCounts[1]);
	outBounds[0]=make_tuple(k0Min,k0Max);
	
	double k1Min(0), k1Max(0);
	double limitSykVav((min(inSpeciesCounts[0],inSpeciesCounts[1])));
	k1Max=overallRateLimit/limitSykVav;
	k1Min=(1./maxTime)/(limitSykVav);
	outBounds[1]=make_tuple(k1Min,k1Max);

	outBounds[2]=make_tuple(k1Min,k1Max);

	double k3Min(0), k3Max(0);
	k3Max=overallRateLimit/(inSpeciesCounts[1]*inSpeciesCounts[4]);
	k3Min=(1./maxTime)/(inSpeciesCounts[1]*inSpeciesCounts[4]);
	outBounds[3]=make_tuple(k3Min,k3Max);

	double k4Min(0), k4Max(0);
	double limitShpVav(min(inSpeciesCounts[1],inSpeciesCounts[4]));
	k4Max=overallRateLimit/limitShpVav;
	k4Min=(1./maxTime)/limitShpVav;
	outBounds[4]=make_tuple(k4Min,k4Max);

	outBounds[5]=make_tuple(k4Min,k4Max);

	return outBounds;
}



int loadPvavInputs(vector<double>& speciesVector, vector<double>& initParameters, vector<double>&stoppingTimes, vector<tuple<double,double> >& bounds, string inFile){
	int errorOut(0);
	ifstream inData(inFile);
	int indexLoop(0);
	double doubleHold(0);
	double doubleHold2(0);
	int intHold(0);
	inData>>indexLoop;
	if(indexLoop!=speciesVector.size()){
		errorOut=1;
		return errorOut;
	}
	for(int i=0;i<indexLoop;i++){
		inData>>intHold;
		speciesVector[i]=intHold;
	}
	inData>>indexLoop;
	if(initParameters.size()!=indexLoop){
		errorOut=1;
		return errorOut;
	}
	for(int i=0;i<indexLoop;i++){
		inData>>doubleHold;
		initParameters[i]=doubleHold;
	}
	inData>>indexLoop;
	stoppingTimes.resize(indexLoop);
	for(int i=0;i<indexLoop;i++){
		inData>>doubleHold;
		stoppingTimes[i]=doubleHold;
	}
	bounds.resize(initParameters.size());
	for(int i=0;i<(int)bounds.size();i++){
		bounds[i]=make_tuple(initParameters[i]/10.,initParameters[i]*10.);
	}
	return errorOut;

}
