#include <windows.h>
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
#include <filesystem>
#include <string>
#include <map>
#include <mpi.h>
	using namespace std;
#include <boost/random/mersenne_twister.hpp>
	boost::mt19937 generator;
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

	
/* Quick functions for portability. */
	
inline double randPull(){
	return (double)generator()/(double)generator.max();
}

inline int randSite(int size){
	return generator()%size;
}


#include "SIR.h"



int main(int argc, char* argv[]){

    generator.seed(time(NULL));
    string localPath=std::filesystem::current_path();
	string outputFolder=localPath+"\\DataFolder";
    std::filesystem::create_directory(outputFolder);
    outputFolder+="\\";
    boost::normal_distribution<> standardNormal(0,1);
    boost::mt19937 generator2;
    generator2.seed(generator());
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> > normalGenerator(generator,standardNormal);

    string customString("NoDynamicsLater_");
    

    int noise(0);
	//T, I, V, R
	const int numOfSpecies(4);
	vector<double> speciesVector={0.95,.05,1.,0.};
	vector<double> resetSpecies=speciesVector;
	//beta, delta, c, p, gamma
	const int numOfParameters(5);
	tuple<double,double,double,double,double> initParameters=make_tuple(0.05,0.05,0.0005,0.02,0.05);
	vector<double> initBounds={1,.5,.5,0.5,0.5};
	const int numOfParticles(20);
	//Number of PSO iterations
	const int numOfIterations(100);

	//Number of Particle sets to run
	const int numOfRuns(40);

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

	vector<double> stoppingTimes={0,10,20,40,60,80};
	double timeIncrement(0.002);

	map<string,int> styleMap;
	styleMap["RungeKutta"]=0;
	styleMap["Gillespie"]=1;
	int solutionStyle(0);
	SolStruct solutionStructure;
	solutionStructure.timeIncrement=timeIncrement;
	solutionStructure.stoppingTimes=stoppingTimes;

	Particle trueParticle=Particle(numOfParameters,initBounds,interactionFuncts,initParameters);

	
	vector<vector<double> > trueArray(stoppingTimes.size()-1,vector<double>(numOfSpecies,0));
	
	switch(solutionStyle){
		case 0:
		{
			if(noise==0){
				trueArray=generateData(&trueParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
			}
			else{
				const int numSamples(2500);
				for(int sample=0;sample<numSamples;sample++){
					vector<double> baseNormal(numOfSpecies,0);
					for(int i=0;i<(int)baseNormal.size();i++){
						baseNormal[i]=normalGenerator();
					}
					baseNormal=transformInit(baseNormal, inValues, speciesVector);
					baseNormal[0]=min(0.99,baseNormal[0]);
					baseNormal[1]=1-baseNormal[0];//Keeps total number of cells equal to 1
					baseNormal[2]=speciesVector[2]*(1+.2*(2*randPull()-1));
					baseNormal[3]=0;
					vector<vector<double> > noiseData=generateData(&trueParticle,baseNormal,&solutionStructure,styleMap["RungeKutta"]);
					for(int i=0;i<(int)trueArray.size();i++){
						for(int j=0;j<(int)trueArray[i].size();j++){
							trueArray[i][j]+=noiseData[i][j]/(double)numSamples;
						}
					}
				}
			}
		}
		break;
		default:
		return 0;
	}
		
	



	double inDelta(10);
	vector<FuzzyTree> FuzzyStructure(numOfParticles,FuzzyTree(inDelta));

	ofstream outRunge;
	switch(noise){
		case 0:
		outRunge.open(outputFolder+customString+"outRunge_noNoise.txt");
		break;
		case 1:
		outRunge.open(outputFolder+customString+"outRunge_testNoise.txt");
		break;
		default:
		return 0;
	}
	for(int i=0;i<(int)trueParticle.currentSolution.size();i++){
		outRunge<<trueParticle.currentSolution[i]<<" ";
	}
	outRunge<<endl;
	
	double globalBestFitness(1e26);
	double globalBestSolution([numOfParameters]);
	int sizeOfParameterVector(numOfParameters);
	double fitnessCollection([numOfParticles]);
	double parameterPassVector[(sizeOfParameterVector)];
	double parameterMatrixHold[(sizeOfParameterVector)*numOfParticles];
	
	int nTasks(-1);
	int taskID(-1);
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskID);
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	
	
	

	for(int run=0;run<numOfRuns;run++){

		FuzzyTree fuzzyStruct(FuzzyStructure[taskID]);
		Particle threadParticle=Particle(numOfParameters,initBounds,interactionFuncts,initParameters);

		for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
			threadParticle.currentSolution[i]=randPull()*initBounds[i];
		}
		

		switch(solutionStyle){
			case 0:
			{
				//Initialize

				speciesVector=resetSpecies;
				vector<vector<double> > testArray=generateData(&threadParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);

				threadParticle.currentFitness=fitnessFunction(trueArray,testArray);
				threadParticle.bestFitness=threadParticle.currentFitness;
				for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
					parameterPassVector[i]=threadParticle.bestSolution[i];
				}
				
				
				MPI_Gather(&threadParticle.currentFitness, 1, MPI_DOUBLE, fitnessCollection, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Gather(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, parameterMatrixHold, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				
				if(taskID==0){
					checkForNewGlobalBest(fitnessCollection, parameterMatrixHold, parameterPassVector, numOfParticles, globalBestFitness);
				}
				
				MPI_Bcast(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);
				
				threadParticle.performUpdate(&generator,globalBestSolutions,&fuzzyStruct);
				

				//Iterate
				for(int iteration=0;iteration<numOfIterations;iteration++){

					speciesVector=resetSpecies;
					vector<vector<double> > testArray=generateData(&threadParticle,speciesVector,&fuzzyStruct,styleMap["RungeKutta"]);

					threadParticle.currentFitness=fitnessFunction(trueArray,testArray);
					if(threadParticle.currentFitness<threadParticle.bestFitness){
						threadParticle.bestSolution=threadParticle.currentSolution;
						threadParticle.bestFitness=threadParticle.currentFitness;
					}
					for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
						parameterPassVector[i]=threadParticle.bestSolution[i];
					}

					MPI_Gather(&threadParticle.currentFitness, 1, MPI_DOUBLE, fitnessCollection, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Gather(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, parameterMatrixHold, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					if(taskID==0){
						checkForNewGlobalBest(fitnesscollection, parameterMatrixHold, parameterPassVector, numOfParticles, globalBestFitness);
					}
					MPI_Bcast(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);

					for(int particle=0;particle<(int)particleList.size();particle++){
						particleList[particle].performUpdate(&generator,globalBestSolutions,&FuzzyStructure[particle]);
					}
				}
				
			}
			break;
		}

		for(int i=0;i<(int)globalBestSolutions.size();i++){
			outRunge<<globalBestSolutions[i]<<" ";
		}
		outRunge<<endl;
	}

	outRunge.close();
    

    return 0;
}






