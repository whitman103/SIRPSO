#include <chrono>
#include <ctime>
#include <map>
#include <math.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <functional>
#include <map>
#include <mpi.h>
	using namespace std;
#include <boost/random/mersenne_twister.hpp>
	boost::mt19937 generator;
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/lognormal_distribution.hpp>

	
/* Quick functions for portability. */
	
inline double randPull(){
	return (double)generator()/(double)generator.max();
}

inline int randSite(int size){
	return generator()%size;
}


#include "GillespieFunctions.h"
#include "SIR_MPI.h"
#include "pVavInteractions.h"

int loadPvavInputs(vector<double>& speciesVector, vector<double>& initParameters, vector<double>& stoppingTimes, string inFile);
int loadPvavInputs(vector<double>& speciesVector, vector<double>& initParameters, vector<double>&stoppingTimes, vector<tuple<double,double> >& bounds, string inFile);
tuple<vector<double>,vector<double> > loadExperimentalData(string inFile);
vector<tuple<double,double> > pVavSetBounds(vector<double>& inSpeciesCounts,vector<double>& inTimes);

int main(int argc, char** argv){

    generator.seed(time(NULL));
    /*string localPath=std::filesystem::current_path();
	string outputFolder=localPath+"\\DataFolder";
    std::filesystem::create_directory(outputFolder);*/
	string outputFolder="DataFolder_ODEMeans_"+string(argv[3]);
	mkdir(outputFolder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    outputFolder+="//";

    string customString("extrinsicNoise_");
    

    bool exNoise(true);
	//T, I, V, R
	const int numOfSpecies(6);
	
	vector<double> speciesVector={600,200,0,0,60,0};

	//k0, k1, k2, k3, k4, k5
	const int numOfParameters(6);
	vector<double> initParameters={0.008,0.1,1.0,0.013904,0.05,0.07};
	vector<double> stoppingTimes={0,100};
	vector<tuple<double,double> > twoBounds;
	int errorFlag=loadPvavInputs(speciesVector,initParameters,stoppingTimes,twoBounds, argv[2]);
	twoBounds=pVavSetBounds(speciesVector,stoppingTimes);

	if(errorFlag==1){
		cout<<"loading data didn't work"<<endl;
		return 0;
	}
	
	int scalingFactor(1);
	transform(speciesVector.begin(),speciesVector.end(),speciesVector.begin(),bind(std::multiplies<double>(),std::placeholders::_1,scalingFactor));
	vector<int> intSpecies(speciesVector.begin(),speciesVector.end());
	vector<double> resetSpecies=speciesVector;
	vector<int> intSpeciesReset=intSpecies;
	int numOfParticles(stoi(argv[1]));
	//Number of PSO iterations
	const int numOfIterations(200);
	//Number of Gillespie samples to use for distributions
	const int numOfSamples(200);

	//Number of Particle sets to run
	const int numOfRuns(200);

	//Generates cholesky matrix to produce lognormal distributions
	vector<vector<double> > inValues;
	vector<double> means;
	loadCovariance(means, inValues,"outCov");
	//vector<vector<double> > outDecomp=generateCholesky(inValues);


	vector<double (*)(Particle*, vector<double>&)> interactionFuncts={dSyk,dVav,dSV,dpVav,dSHP1,dSHP1Vav};

	boost::normal_distribution<> standardNormal(0,1);
    boost::mt19937 generator2;
    generator2.seed(generator());
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> > normalGenerator(generator2,standardNormal);
	boost::mt19937 exGenerator;
	exGenerator.seed(generator2());

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

	ofstream fuckOff("fuckoff.txt");
	for(int i=0;i<100;i++){
		fuckOff<<extrinsicNoiseGenerators[0]()<<" ";
	}
	fuckOff.close();
	
    
	double timeIncrement(0.002);

	map<string,int> styleMap;
	styleMap["RungeKutta"]=0;
	styleMap["Gillespie"]=1;

	int solutionStyle(0);

	SolStruct solutionStructure;
	solutionStructure.timeIncrement=timeIncrement;
	solutionStructure.stoppingTimes=stoppingTimes;

	Particle trueParticle=Particle(numOfParameters,initParameters,twoBounds,interactionFuncts,scalingFactor);

	vector<vector<double> > trueArray(stoppingTimes.size()-1,vector<double>(numOfSpecies,0));
    vector<vector<vector<int> > > trueDistributions(stoppingTimes.size()-1,vector<vector<int> > (numOfSpecies,vector<int> (numOfRuns,0)));

	switch(solutionStyle){
		case 0:
		{
			if(!exNoise){
				resetSpecies=speciesVector;
				trueArray=generateData(&trueParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
			}
			else{
				
				for(int sample=0;sample<numOfSamples;sample++){
					speciesVector=resetSpecies;
					speciesVector[0]=extrinsicNoiseGenerators[0]();
					speciesVector[1]=extrinsicNoiseGenerators[1]();
					speciesVector[4]=extrinsicNoiseGenerators[2]();
					vector<vector<double> > noiseData=generateData(&trueParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
					for(int i=0;i<(int)trueArray.size();i++){
						for(int j=0;j<(int)trueArray[i].size();j++){
							trueArray[i][j]+=noiseData[i][j]/(double)numOfSamples;
						}
					}
				}
			}
		}
		break;
		case 1:
		{
			if(!exNoise){
				intSpecies=intSpeciesReset;
				Gillespie ReactionObject1("pVavCoeffs");
				ReactionObject1.initializeData("pVavConsts",ReactionObject1.reactConsts,intSpeciesReset);
				trueDistributions=pVav_Gillespie(&trueParticle, &ReactionObject1, stoppingTimes, intSpecies, numOfSamples);
			}
			else{
				intSpecies=intSpeciesReset;
				Gillespie ReactionObject1("pVavCoeffs");
				ReactionObject1.initializeData("pVavConsts",ReactionObject1.reactConsts,intSpeciesReset);
				trueDistributions=pVav_Gillespie(&trueParticle, &ReactionObject1, stoppingTimes, intSpecies, numOfSamples);//Needs to be updated for extrinsic noise

			}
		}
		break;
		default:
		return 0;
	}
	
	double inDelta(10);
	vector<FuzzyTree> FuzzyStructure(numOfParticles,FuzzyTree(inDelta));

	ofstream outRunge;
	if(!exNoise){
		outRunge.open(outputFolder+customString+"outRunge_noNoise.txt");
	} else{
		outRunge.open(outputFolder+customString+"outRunge_testNoise.txt");
	}
	for(int i=0;i<(int)trueParticle.currentSolution.size();i++){
		outRunge<<trueParticle.currentSolution[i]<<" ";
	}
	outRunge<<endl;
	
	
	
	int sizeOfParameterVector(initParameters.size());
	double fitnessCollection[numOfParticles];
	double parameterPassVector[sizeOfParameterVector];
	double parameterMatrixHold[sizeOfParameterVector*numOfParticles];
	
	int nTasks(-1);
	int taskID(-1);

	Gillespie threadReaction("pVavCoeffs");
	threadReaction.initializeData("pVavConsts",threadReaction.reactConsts,intSpeciesReset);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskID);
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	
	ofstream outFile;
	if(taskID==0){
		outFile.open(outputFolder+customString+"bestParticles"+".txt");
	}
	
	generator.seed(taskID);
	for(int i=0;i<(int)exNoiseEngines.size();i++){
		exNoiseEngines[i].seed(taskID*nTasks*i);
	}

	for(int run=0;run<numOfRuns;run++){
		double globalBestFitness(1e26);

		FuzzyTree fuzzyStruct(FuzzyStructure[taskID]);
		Particle threadParticle=Particle(numOfParameters,initParameters,twoBounds,interactionFuncts,scalingFactor);
		for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
			auto[lowerBound,upperBound]=twoBounds[i];
			threadParticle.currentSolution[i]=pow(10,randPull()*log10(upperBound/lowerBound))*lowerBound;
		}
		
        vector<vector<vector<int> > > testDistributions;
		switch(solutionStyle){
			case 0:
			{
				//Initialize
				vector<vector<double> > testArray(stoppingTimes.size()-1,vector<double>(numOfSpecies,0));
				vector<vector<double> > resetTestArray=testArray;
				speciesVector=resetSpecies;
				if(!exNoise){
					testArray=generateData(&threadParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
				}
				else{
					for(int sample=0;sample<numOfSamples;sample++){
						speciesVector=resetSpecies;
						speciesVector[0]=extrinsicNoiseGenerators[0]();
						speciesVector[1]=extrinsicNoiseGenerators[1]();
						speciesVector[4]=extrinsicNoiseGenerators[2]();
						vector<vector<double> > noiseData=generateData(&threadParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
						for(int i=0;i<(int)testArray.size();i++){
							for(int j=0;j<(int)testArray[i].size();j++){
								testArray[i][j]+=noiseData[i][j]/(double)numOfSamples;
							}
						}
					}
				}

				threadParticle.currentFitness=fitnessFunction(trueArray,testArray);
				threadParticle.bestFitness=threadParticle.currentFitness;
				threadParticle.bestSolution=threadParticle.currentSolution;
				for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
					parameterPassVector[i]=threadParticle.bestSolution[i];
				}
				
				
				MPI_Gather(&threadParticle.currentFitness, 1, MPI_DOUBLE, fitnessCollection, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Gather(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, parameterMatrixHold, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				
				if(taskID==0){
					checkForNewGlobalBest(fitnessCollection, parameterMatrixHold, parameterPassVector, numOfParticles, globalBestFitness, numOfParameters);
				}
				
				MPI_Bcast(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);
				
				threadParticle.twoBoundPerformUpdate(&generator,parameterPassVector,&fuzzyStruct);
				

				//Iterate
				for(int iteration=0;iteration<numOfIterations;iteration++){
					speciesVector=resetSpecies;
					testArray=resetTestArray;
					if(!exNoise){
						testArray=generateData(&threadParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
					}
					else{
						for(int sample=0;sample<numOfSamples;sample++){
							speciesVector=resetSpecies;
							speciesVector[0]=extrinsicNoiseGenerators[0]();
							speciesVector[1]=extrinsicNoiseGenerators[1]();
							speciesVector[4]=extrinsicNoiseGenerators[2]();
							vector<vector<double> > noiseData=generateData(&threadParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
							for(int i=0;i<(int)testArray.size();i++){
								for(int j=0;j<(int)testArray[i].size();j++){
									testArray[i][j]+=noiseData[i][j]/(double)numOfSamples;
								}
							}
						}
					}
					
					threadParticle.currentFitness=fitnessFunction(trueArray,testArray);
					if(threadParticle.currentFitness<threadParticle.bestFitness){
						threadParticle.bestSolution=threadParticle.currentSolution;
						threadParticle.bestFitness=threadParticle.currentFitness;
					}
					for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
						parameterPassVector[i]=threadParticle.currentSolution[i];
					}

					MPI_Gather(&threadParticle.currentFitness, 1, MPI_DOUBLE, fitnessCollection, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Gather(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, parameterMatrixHold, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					if(taskID==0){
						checkForNewGlobalBest(fitnessCollection, parameterMatrixHold, parameterPassVector, numOfParticles, globalBestFitness,numOfParameters);
					}
					MPI_Bcast(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);

					threadParticle.twoBoundPerformUpdate(&generator,parameterPassVector,&fuzzyStruct);
					
				}
				
			}
			break;
			/*case 1:
			{
				vector<vector<double> > testArray;
				
                vector<vector<vector<int> > > testDistributions;
			
				intSpecies=intSpeciesReset;
				if(!exNoise){
					intSpecies=intSpeciesReset;
					testDistributions=pVav_Gillespie(&threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples);
				}
				else{
					intSpecies=intSpeciesReset;
					testDistributions=pVav_Gillespie(&threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples);//needs to be updated for extrinsic noise
				}

				

				threadParticle.currentFitness=fitnessFunction_pVav(trueDistributions,testDistributions,"pVav_Means_MultipleTrans");
				threadParticle.bestFitness=threadParticle.currentFitness;
				threadParticle.bestSolution=threadParticle.currentSolution;
				for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
					parameterPassVector[i]=threadParticle.bestSolution[i];
				}
				
				MPI_Gather(&threadParticle.currentFitness, 1, MPI_DOUBLE, fitnessCollection, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Gather(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, parameterMatrixHold, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				
				if(taskID==0){
					checkForNewGlobalBest(fitnessCollection, parameterMatrixHold, parameterPassVector, numOfParticles, globalBestFitness, numOfParameters);
				}
				
				MPI_Bcast(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);
				
				threadParticle.performUpdate(&generator,parameterPassVector,&fuzzyStruct);

				//iterate
				for(int iteration=0;iteration<numOfIterations;iteration++){
					
					if(!exNoise){
						intSpecies=intSpeciesReset;
						testDistributions=pVav_Gillespie(&threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples);
					}
					else{
						intSpecies=intSpeciesReset;
						testDistributions=pVav_Gillespie(&threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples);//needs to be updated for extrinsic noise
					}

					threadParticle.currentFitness=fitnessFunction_pVav(trueDistributions,testDistributions,"pVav_Means_MultipleTrans");
					if(threadParticle.currentFitness<threadParticle.bestFitness){
						threadParticle.bestSolution=threadParticle.currentSolution;
						threadParticle.bestFitness=threadParticle.currentFitness;
					}
					for(int i=0;i<(int)threadParticle.bestSolution.size();i++){
						parameterPassVector[i]=threadParticle.bestSolution[i];
					}

					MPI_Gather(&threadParticle.currentFitness, 1, MPI_DOUBLE, fitnessCollection, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Gather(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, parameterMatrixHold, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					if(taskID==0){
						checkForNewGlobalBest(fitnessCollection, parameterMatrixHold, parameterPassVector, numOfParticles, globalBestFitness,numOfParameters);
					}
					MPI_Bcast(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					
					threadParticle.performUpdate(&generator,parameterPassVector,&fuzzyStruct);


					
				}
			}*/
			default:
			break;
		}
		if(taskID==0){
			outFile<<globalBestFitness<<" ";
			for(int i=0;i<numOfParameters;i++){
				outFile<<parameterPassVector[i]<<" ";
			}
			outFile<<endl;
		}

	}
	outFile.close();
    
	MPI_Finalize();

    return 0;
}

int loadPvavInputs(vector<double>& speciesVector, vector<double>& initParameters, vector<double>& stoppingTimes, string inFile){
	int errorOut(0);
	ifstream inData(inFile);
	int indexLoop(0);
	double doubleHold(0);
	int intHold(0);
	ofstream errorFile("FuckedUp.txt");
	inData>>indexLoop;
	if(indexLoop!=speciesVector.size()){
		errorOut=1;
		errorFile<<"First"<<endl;
		return errorOut;
	}
	for(int i=0;i<indexLoop;i++){
		inData>>intHold;
		speciesVector[i]=intHold;
	}
	inData>>indexLoop;
	if(initParameters.size()!=indexLoop){
		errorOut=1;
		errorFile<<"Second"<<endl;
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
	return errorOut;
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

tuple<vector<double>,vector<double> > loadExperimentalData(string inFile){
	int numOfTimes(5);
	tuple<vector<double>,vector<double> > outMeansAndSDs;
	vector<double> means(numOfTimes);
	vector<double> sds(numOfTimes);
	ifstream readIn(inFile);
	double doubleHold(0);
	for(int i=0;i<5;i++){
		readIn>>doubleHold;
		means[i]=doubleHold;
		readIn>>doubleHold;
		sds[i]=doubleHold;
	}
	get<0>(outMeansAndSDs)=means;
	get<1>(outMeansAndSDs)=sds;
	return outMeansAndSDs;
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








