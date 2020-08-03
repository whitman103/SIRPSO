#include <iostream>
#include <vector>
	using namespace std;
#include <boost/random/mersenne_twister.hpp>
	boost::mt19937 generator;
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "GillespieFunctions.h"
#include "SIR_MPI.h"

int main(){
	
	generator.seed(time(NULL));
    
    boost::normal_distribution<> standardNormal(0,1);
    boost::mt19937 generator2;
    generator2.seed(generator());
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> > normalGenerator(generator2,standardNormal);
	boost::mt19937 exGenerator;
	exGenerator.seed(generator2());

    string customString("NoDynamicsLater_");
    

    bool exNoise(true);
	//T, I, V, R
	const int numOfSpecies(4);
	vector<double> speciesVector={0.95,.05,1.,0.};
	int scalingFactor(500);
	transform(speciesVector.begin(),speciesVector.end(),speciesVector.begin(),bind(std::multiplies<double>(),std::placeholders::_1,scalingFactor));
	vector<int> intSpecies(speciesVector.begin(),speciesVector.end());
	vector<double> resetSpecies=speciesVector;
	vector<int> intSpeciesReset=intSpecies;
	//beta, delta, c, p, gamma
	const int numOfParameters(5);
	string inParams="paramTests.txt";
	tuple<double,double,double,double,double> initParameters=readParameterData(inParams);
	vector<double> initBounds={.5,.1,.1,0.01,0.05};
	//Number of PSO iterations
	const int numOfIterations(100);
	//Number of Gillespie samples to use for distributions
	const int numOfSamples(500);

	//Number of Particle sets to run
	const int numOfRuns(25);

	//Generates cholesky matrix to produce lognormal distributions
	vector<vector<double> > inValues;
	vector<double> means;
	loadCovariance(means, inValues,"outCov");
	
	vector<vector<double> > outDecomp=generateCholesky(inValues);


	vector<double (*)(Particle*,vector<double>&)> interactionFuncts(numOfSpecies);
	interactionFuncts[0]=firstInteraction;
	interactionFuncts[1]=secondInteraction;
	interactionFuncts[2]=thirdInteraction;
	interactionFuncts[3]=fourthInteraction;


	double timeIncrement(0.002);

	map<string,int> styleMap;
	styleMap["RungeKutta"]=0;
	styleMap["Gillespie"]=1;

	int solutionStyle(1);
	
	vector<double> stoppingTimes;
	for(int i=0;i<10;i++){
		stoppingTimes.push_back(i*10.);
	}

	SolStruct solutionStructure;
	solutionStructure.timeIncrement=timeIncrement;
	solutionStructure.stoppingTimes=stoppingTimes;

	Particle trueParticle=Particle(numOfParameters,initBounds,interactionFuncts,initParameters,scalingFactor);

	
	vector<vector<double> > trueArray(stoppingTimes.size()-1,vector<double>(numOfSpecies,0));
	vector<vector<double> > trueVar=trueArray;
	
	
	intSpecies=intSpeciesReset;
	Gillespie ReactionObject1("SIRCoeffs");
	ReactionObject1.initializeData("SIRConsts",ReactionObject1.reactConsts,intSpeciesReset);
	
	generateDistributions(&trueParticle, &ReactionObject1, stoppingTimes, intSpecies, numOfSamples, &exGenerator, "test.txt");
	
	return 0;
}