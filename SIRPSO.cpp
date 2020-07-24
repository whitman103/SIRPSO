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
#include <functional>
#include <string>
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


#include "SIR_MPI.h"
#include "GillespieFunctions.h"

void sendVectorToArray(vector<double> inVec, double* inArray);
void sendArrayToVector(vector<double>& inVec, double* inArray);

int main(){

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
    

    for(int noise=0;noise<2;noise++){
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
        const int numOfParticles(5);
        //Number of PSO iterations
        const int numOfIterations(100);

        const int numOfSamples(500);

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
        int solutionStyle(1);
        SolStruct solutionStructure;
        solutionStructure.timeIncrement=timeIncrement;
        solutionStructure.stoppingTimes=stoppingTimes;

        Particle trueParticle=Particle(numOfParameters,initBounds,interactionFuncts,initParameters,scalingFactor);

        
        vector<vector<double> > trueArray(stoppingTimes.size()-1,vector<double>(numOfSpecies,0));
        vector<vector<double> > trueVar=trueArray;
        
        switch(solutionStyle){
            case 0:
            {/*
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
                }*/
            }
            break;
            case 1:
            {
                Gillespie ReactionObject1("SIRCoeffs");
                intSpecies=intSpeciesReset;
                ReactionObject1.initializeData("SIRConsts",ReactionObject1.reactConsts,intSpeciesReset);
                vector<double> resetConsts=ReactionObject1.reactConsts;
                tie(trueArray,trueVar)=generateGillespieData(&trueParticle, &ReactionObject1, stoppingTimes, intSpecies, numOfSamples);
            }
            break;
            default:
            return 0;
        }
            
        


        double inDelta(10);
        vector<FuzzyTree> FuzzyStructure(numOfParticles,FuzzyTree(inDelta));

        Gillespie threadReaction("SIRCoeffs");
        intSpecies=intSpeciesReset;
        threadReaction.initializeData("SIRConsts",threadReaction.reactConsts,intSpecies);

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

        for(int i=0;i<(int)trueArray.size();i++){
            for(int j=0;j<(int)trueArray[i].size();j++){
                outRunge<<trueArray[i][j]<<" ";
            }
            outRunge<<endl;
        }
        outRunge<<"Above is True data"<<endl;

        ofstream monitorFile("fitnessMonitor.txt");

        vector<vector<double> > outArray;

        for(int run=0;run<numOfRuns;run++){


            vector<Particle> particleList(numOfParticles,Particle(numOfParameters,initBounds,interactionFuncts,initParameters,scalingFactor));
            for(int i=0;i<(int)particleList.size();i++){
                vector<double> interParams(numOfParameters);
                for(int param=0;param<(int)interParams.size();param++){
                    //particleList[i].currentSolution[param]=randPull()*initBounds[param];
                }
            }
            
            
            vector<double> globalBestSolutions(numOfParameters,0);
            double parameterArray[numOfParameters];

            switch(solutionStyle){
                case 0:
                {
                    int bestParticle(0);
                    //Initialize
                    vector<double> fitnessHolder(numOfParticles,0);

                    for(int particle=0;particle<(int)particleList.size();particle++){

                        speciesVector=resetSpecies;
                        Particle* interParticle=&particleList[particle];
                        vector<vector<double> > testArray=generateData(interParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);

                        (*interParticle).currentFitness=fitnessFunction(trueArray,testArray);
                        (*interParticle).bestFitness=(*interParticle).currentFitness;
                        fitnessHolder[particle]=(*interParticle).currentFitness;
                    }

                    bestParticle=min_element(fitnessHolder.begin(),fitnessHolder.end())-fitnessHolder.begin();
                    globalBestSolutions=particleList[bestParticle].currentSolution;
                    sendVectorToArray(globalBestSolutions,parameterArray);
                    
                    for(int particle=0;particle<(int)particleList.size();particle++){
                        particleList[particle].performUpdate(&generator,parameterArray,&FuzzyStructure[particle]);
                    }
                    sendArrayToVector(globalBestSolutions,parameterArray);

                    //Iterate
                    for(int iteration=0;iteration<numOfIterations;iteration++){
                        for(int particle=0;particle<(int)particleList.size();particle++){

                            speciesVector=resetSpecies;
                            Particle* interParticle=&particleList[particle];
                            vector<vector<double> > testArray=generateData(interParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);

                            (*interParticle).currentFitness=fitnessFunction(trueArray,testArray);
                            if(((*interParticle).currentFitness<(*interParticle).bestFitness)){
                                (*interParticle).bestSolution=(*interParticle).currentSolution;
                                (*interParticle).bestFitness=(*interParticle).currentFitness;
                            }
                            fitnessHolder[particle]=(*interParticle).currentFitness;
                        }

                        bestParticle=min_element(fitnessHolder.begin(),fitnessHolder.end())-fitnessHolder.begin();
                        sendVectorToArray(globalBestSolutions,parameterArray);
                    
                        for(int particle=0;particle<(int)particleList.size();particle++){
                            particleList[particle].performUpdate(&generator,parameterArray,&FuzzyStructure[particle]);
                        }
                        sendArrayToVector(globalBestSolutions,parameterArray);
                    }
                }
                break;
                case 1:
                {
                    vector<vector<double> > testArray;
                    vector<vector<double> > testVar;
                    vector<double> fitnessHolder(numOfParticles,0);
                    

                    int bestParticle(0);

                    for(int particle=0;particle<(int)particleList.size();particle++){
                        Particle* threadParticle=&particleList[particle];

                        intSpecies=intSpeciesReset;
                        tie(testArray,testVar)=generateGillespieData(threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples);

                        

                        (*threadParticle).currentFitness=fitnessFunction(trueArray,testArray);
                        (*threadParticle).bestFitness=(*threadParticle).currentFitness;
                        (*threadParticle).bestSolution=(*threadParticle).currentSolution;
                        fitnessHolder[particle]=(*threadParticle).currentFitness;
                    }
                    
                    bestParticle=min_element(fitnessHolder.begin(),fitnessHolder.end())-fitnessHolder.begin();
                    
                    globalBestSolutions=particleList[bestParticle].currentSolution;
                    sendVectorToArray(globalBestSolutions,parameterArray);
                    
                    for(int particle=0;particle<(int)particleList.size();particle++){
                        particleList[particle].performUpdate(&generator,parameterArray,&FuzzyStructure[particle]);
                    }
                    sendArrayToVector(globalBestSolutions,parameterArray);
                    

                    for(int iteration=0;iteration<numOfIterations;iteration++){
                        for(int particle=0;particle<(int)particleList.size();particle++){
                            Particle* threadParticle=&particleList[particle];

                            intSpecies=intSpeciesReset;
                            tie(testArray,testVar)=generateGillespieData(threadParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples);

                            

                            (*threadParticle).currentFitness=fitnessFunction(trueArray,testArray);
                            fitnessHolder[particle]=(*threadParticle).currentFitness;
                            if(((*threadParticle).currentFitness<(*threadParticle).bestFitness)){
                                (*threadParticle).bestSolution=(*threadParticle).currentSolution;
                                (*threadParticle).bestFitness=(*threadParticle).currentFitness;
                            }
                            fitnessHolder[particle]=(*threadParticle).currentFitness;
                        }

                        bestParticle=min_element(fitnessHolder.begin(),fitnessHolder.end())-fitnessHolder.begin();
                        sendVectorToArray(globalBestSolutions,parameterArray);
                        
                    
                        for(int particle=0;particle<(int)particleList.size();particle++){
                            particleList[particle].performUpdate(&generator,parameterArray,&FuzzyStructure[particle]);
                        }
                        sendArrayToVector(globalBestSolutions,parameterArray);

                        
                        
                        for(int i=0;i<numOfParticles;i++){
                            monitorFile<<fitnessHolder[i]<<" ";
                        }
                        monitorFile<<endl;
                        
                    }
                }
                break;
            }

            for(int i=0;i<(int)globalBestSolutions.size();i++){
                outRunge<<globalBestSolutions[i]<<" ";
            }
            outRunge<<endl;
            Particle outParticle=Particle(numOfParameters,initBounds,interactionFuncts,initParameters,scalingFactor);
            outParticle.currentSolution=globalBestSolutions;
            intSpecies=intSpeciesReset;
            vector<vector<double> > testArray;
            vector<vector<double> > testVar;
            tie(testArray,testVar)=generateGillespieData(&outParticle, &threadReaction, stoppingTimes, intSpecies, numOfSamples);

            for(int i=0;i<(int)testArray.size();i++){
                for(int j=0;j<(int)testArray[i].size();j++){
                    outRunge<<testArray[i][j]<<" ";
                }
                outRunge<<endl;
            }
            outRunge<<endl;
            
        }

        outRunge.close();
    }
    

    return 0;
}

void sendVectorToArray(vector<double> inVec, double* inArray){
    for(int i=0;i<(int)inVec.size();i++){
        inArray[i]=inVec[i];
    }
}

void sendArrayToVector(vector<double>& inVec, double* inArray){
    for(int i=0;i<(int)inVec.size();i++){
        inVec[i]=inArray[i];
    }
}





