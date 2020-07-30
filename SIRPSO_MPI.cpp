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

	
/* Quick functions for portability. */
	
inline double randPull(){
	return (double)generator()/(double)generator.max();
}

inline int randSite(int size){
	return generator()%size;
}


#include "GillespieFunctions.h"
#include "SIR_MPI.h"



int main(int argc, char** argv){

    generator.seed(time(NULL));
    /*string localPath=std::filesystem::current_path();
	string outputFolder=localPath+"\\DataFolder";
    std::filesystem::create_directory(outputFolder);*/
	string outputFolder="DataFolder_ODENoise"+string(argv[3]);
	mkdir(outputFolder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    outputFolder+="//";
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
	tuple<double,double,double,double,double> initParameters=make_tuple(0.1,0.05,0.05,0.0025,0.025);
	vector<double> initBounds={.5,.1,.1,0.01,0.05};
	int numOfParticles(stoi(argv[1]));
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

	string argString(argv[2]);
	vector<double> stoppingTimes=readVectorFile(argString+".txt");
	double timeIncrement(0.002);

	map<string,int> styleMap;
	styleMap["RungeKutta"]=0;
	styleMap["Gillespie"]=1;

	int solutionStyle(1);

	SolStruct solutionStructure;
	solutionStructure.timeIncrement=timeIncrement;
	solutionStructure.stoppingTimes=stoppingTimes;

	Particle trueParticle=Particle(numOfParameters,initBounds,interactionFuncts,initParameters,scalingFactor);

	
	vector<vector<double> > trueArray(stoppingTimes.size()-1,vector<double>(numOfSpecies,0));
	vector<vector<double> > trueVar=trueArray;
	
	switch(solutionStyle){
		case 0:
		{
			if(!exNoise){
				resetSpecies=speciesVector;
				trueArray=generateData(&trueParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
			}
			else{
				const int numSamples(2500);
				for(int sample=0;sample<numSamples;sample++){
					vector<double> baseNormal(numOfSpecies,0);
					for(int i=0;i<(int)baseNormal.size();i++){
						baseNormal[i]=normalGenerator();
					}
					baseNormal=transformInit(baseNormal, inValues, speciesVector, &generator);
					transform(baseNormal.begin(),baseNormal.end(),baseNormal.begin(),bind(std::multiplies<double>(),std::placeholders::_1,(double)scalingFactor));
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
		case 1:
		{
			if(!exNoise){
				intSpecies=intSpeciesReset;
				Gillespie ReactionObject1("SIRCoeffs");
				ReactionObject1.initializeData("SIRConsts",ReactionObject1.reactConsts,intSpeciesReset);
				tie(trueArray,trueVar)=generateGillespieData(&trueParticle, &ReactionObject1, stoppingTimes, intSpecies, numOfSamples);
			}
			else{
				intSpecies=intSpeciesReset;
				Gillespie ReactionObject1("SIRCoeffs");
				ReactionObject1.initializeData("SIRConsts",ReactionObject1.reactConsts,intSpeciesReset);
				tie(trueArray,trueVar)=generateGillespieData(&trueParticle, &ReactionObject1, stoppingTimes, intSpecies, numOfSamples, &exGenerator);
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
	
	
	
	int sizeOfParameterVector(numOfParameters);
	double fitnessCollection[numOfParticles];
	double parameterPassVector[sizeOfParameterVector];
	double parameterMatrixHold[sizeOfParameterVector*numOfParticles];
	
	int nTasks(-1);
	int taskID(-1);

	Gillespie threadReaction("SIRCoeffs");
	threadReaction.initializeData("SIRConsts",threadReaction.reactConsts,intSpeciesReset);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskID);
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	
	ofstream outFile;
	ofstream monitorFile;
	if(taskID==0){
		outFile.open(outputFolder+customString+"bestParticles"+".txt");
	}
	

	for(int run=0;run<numOfRuns;run++){
		double globalBestFitness(1e26);

		FuzzyTree fuzzyStruct(FuzzyStructure[taskID]);
		Particle threadParticle=Particle(numOfParameters,initBounds,interactionFuncts,initParameters,scalingFactor);
		generator.seed(taskID);
		for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
			threadParticle.currentSolution[i]=randPull()*initBounds[i];
		}
		

		switch(solutionStyle){
			case 0:
			{
				//Initialize
				vector<vector<double> > testArray(stoppingTimes.size()-1,vector<double>(numOfSpecies,0));
				speciesVector=resetSpecies;
				if(!exNoise){
					testArray=generateData(&threadParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
				}
				else{
					const int numSamples(2500);
					for(int sample=0;sample<numSamples;sample++){
						vector<double> baseNormal(numOfSpecies,0);
						for(int i=0;i<(int)baseNormal.size();i++){
							baseNormal[i]=normalGenerator();
						}
						baseNormal=transformInit(baseNormal, inValues, speciesVector,&generator);
						vector<vector<double> > noiseData=generateData(&threadParticle,baseNormal,&solutionStructure,styleMap["RungeKutta"]);
						for(int i=0;i<(int)testArray.size();i++){
							for(int j=0;j<(int)testArray[i].size();j++){
								testArray[i][j]+=noiseData[i][j]/(double)numSamples;
							}
						}
					}
				}

				threadParticle.currentFitness=fitnessFunction(trueArray,testArray);
				cout<<"first fitness is: "<<threadParticle.currentFitness<<endl;
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
				

				//Iterate
				for(int iteration=0;iteration<numOfIterations;iteration++){
					speciesVector=resetSpecies;
					if(!exNoise){
						testArray=generateData(&threadParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
					}
					else{
						const int numSamples(2500);
						for(int sample=0;sample<numSamples;sample++){
							vector<double> baseNormal(numOfSpecies,0);
							for(int i=0;i<(int)baseNormal.size();i++){
								baseNormal[i]=normalGenerator();
							}
							baseNormal=transformInit(baseNormal, inValues, speciesVector, &generator);
							vector<vector<double> > noiseData=generateData(&threadParticle,baseNormal,&solutionStructure,styleMap["RungeKutta"]);
							for(int i=0;i<(int)testArray.size();i++){
								for(int j=0;j<(int)testArray[i].size();j++){
									testArray[i][j]+=noiseData[i][j]/(double)numSamples;
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
						for(int i=0;i<numOfParticles;i++){
							monitorFile<<fitnessCollection[i]<<" ";
						}
						monitorFile<<endl;
					}
					MPI_Bcast(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);

					
					threadParticle.performUpdate(&generator,parameterPassVector,&fuzzyStruct);
					
				}
				
			}
			break;
			case 1:
			{
				vector<vector<double> > testArray;
				vector<vector<double> > testVar;
			
				intSpecies=intSpeciesReset;
				if(!exNoise){
					intSpecies=intSpeciesReset;
					tie(testArray,testVar)=generateGillespieData(&threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples);
				}
				else{
					intSpecies=intSpeciesReset;
					tie(testArray,testVar)=generateGillespieData(&threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples, &exGenerator);
				}

				

				threadParticle.currentFitness=fitnessFunction(trueArray,testArray,testVar,trueVar);
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
					intSpecies=intSpeciesReset;
					if(!exNoise){
						intSpecies=intSpeciesReset;
						tie(testArray,testVar)=generateGillespieData(&threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples);
					}
					else{
						intSpecies=intSpeciesReset;
						tie(testArray,testVar)=generateGillespieData(&threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples, &exGenerator);
					}

					threadParticle.currentFitness=fitnessFunction(trueArray,testArray,trueVar,testVar);
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
						checkForNewGlobalBest(fitnessCollection, parameterMatrixHold, parameterPassVector, numOfParticles, globalBestFitness,numOfParameters);
					}
					MPI_Bcast(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					
					threadParticle.performUpdate(&generator,parameterPassVector,&fuzzyStruct);


					
				}
			}
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
	monitorFile.close();
    
	MPI_Finalize();
	

    return 0;
}






