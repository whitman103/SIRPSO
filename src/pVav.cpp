#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <tuple>
#include <string>
#include <filesystem>
#include <functional>
using namespace std;
#include <boost/random/mersenne_twister.hpp>
boost::mt19937 generator;


#include "GillespieFunctions.h"
#include "SIR_MPI.h"
#include "pVavInteractions.h"


int main() {


    generator.seed(time(NULL));

    string localPath = std::filesystem::current_path();
    string outputFolder = localPath + "\\DataFolder";
    std::filesystem::create_directory(outputFolder);
    outputFolder += "\\";


    const int numOfParameters(6);
    const int numOfSpecies(6);
    vector<double> paramValues = {0.003, 0.01, 0.1, 0.0013904, 0.2, 0.07};
    vector<double> speciesVec = {600, 200, 0, 0, 60, 0};
    vector<int> intSpeciesReset(speciesVec.begin(), speciesVec.end());
    vector<double (*)(Particle*, vector<double>&)> interactionFunctions = {dSyk, dVav, dSV, dpVav, dSHP1, dSHP1Vav};
    int scalingSize(1);


    Particle pVavParticle = Particle(numOfParameters, paramValues, interactionFunctions, scalingSize);

    ofstream outData(outputFolder + "testData.txt");
    double timeIncrement(0.001);
    for (double t = 0; t < 2; t += timeIncrement) {
        pVav_RungeKutta(&pVavParticle, speciesVec, 0, timeIncrement, 0.001);
        outData << (t + 1) << " ";
        for (int i = 0; i < (int)speciesVec.size(); i++) {
            outData << speciesVec[i] << " ";
        }
        outData << endl;
    }
    outData.close();

    /*

    Gillespie ReactionObject("pVavCoeffs");

    vector<int> specNum=intSpeciesReset;
    ReactionObject.initializeData("pVavConsts",ReactionObject.reactConsts,intSpeciesReset);

    const int numOfRuns(2000);
    vector<double> reportTimes={0,300};

    vector<vector<vector<int> > > outData(reportTimes.size()-1,vector<vector<int> > (specNum.size(), vector<int> (numOfRuns,0)));

    const int numOfDists(10);
    vector<int> pVavValues(numOfDists,0);
    for(int i=0;i<(int)pVavValues.size();i++){
    	pVavValues[i]=(i+1)*15;
    }



    for(int pVavIndex=0;pVavIndex<(int)pVavValues.size();pVavIndex++){
    	intSpeciesReset[0]=pVavValues[pVavIndex];

    	for(int run=0;run<numOfRuns;run++){
    		double currentTime(0);
    		int reportIndex(0);
    		specNum=intSpeciesReset;// This is where you introduce noise in the initial conditions

    		do{
    			tuple<int,double> hold=ReactionObject.PerformTimeStep2(specNum);
    			currentTime+=get<1>(hold);

    			if(get<1>(hold)<0){
    				for(int index=reportIndex;index<(int)reportTimes.size();index++){
    					for(int i=0;i<(int)specNum.size();i++){
    						outData[index][i][run]=specNum[i];
    					}
    				}
    				reportIndex=(int)reportTimes.size();
    			}
    			else{
    				ReactionObject.specChange(specNum,get<0>(hold),ReactionObject.changeCoeffs);
    				if(currentTime>reportTimes[reportIndex+1]){
    					for(int i=0;i<(int)specNum.size();i++){
    						outData[reportIndex][i][run]=specNum[i];
    					}
    					reportIndex++;
    				}
    			}
    		}while(reportIndex<(int)reportTimes.size()-1);
    	}

    	ofstream outFile(outputFolder+"pVavData_SYK_TurnOn_"+to_string(pVavValues[pVavIndex])+".txt");
    	for(int run=0;run<numOfRuns;run++){
    		for(int time=0;time<(int)reportTimes.size()-1;time++){
    			for(int species=0;species<(int)specNum.size();species++){
    				outFile<<outData[time][species][run]<<" ";
    			}
    		}
    		outFile<<endl;
    	}
    	outFile.close();
    }

    */
    return 0;
}

