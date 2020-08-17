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


int main(){
	
	
	generator.seed(time(NULL));
	
	string localPath=std::filesystem::current_path();
	string outputFolder=localPath+"\\DataFolder";
    std::filesystem::create_directory(outputFolder);
	outputFolder+="\\";

	
	const int numOfParameters(6);
	const int numOfSpecies(6);
	vector<double> paramValues={0.008,0.1,1.0,0.013904,0.05,0.07};
	vector<double> speciesVec={576,200,0,0,384,0};
	vector<double (*)(Particle*, vector<double>&)> interactionFunctions={dSyk,dVav,dSV,dpVav,dSHP1,dSHP1Vav};
	int scalingSize(1);
	
	
	Particle pVavParticle=Particle(numOfParameters,paramValues,interactionFunctions,scalingSize);
	pVav_RungeKutta(&pVavParticle,speciesVec,0,3*3600,0.02);
	
	for(int i=0;i<(int)speciesVec.size();i++){
		cout<<speciesVec[i]<<" ";
	}
	cout<<endl;
	
	
	
	Gillespie ReactionObject("pVavCoeffs");
	vector<int> intSpeciesReset={576,200,0,0,384,0};
	vector<int> specNum=intSpeciesReset;
	ReactionObject.initializeData("pVavConsts",ReactionObject.reactConsts,intSpeciesReset);
	
	const int numOfRuns(500);
	vector<double> reportTimes={0,3*3600};
	
	vector<vector<vector<int> > > outData(reportTimes.size()-1,vector<vector<int> > (specNum.size(), vector<int> (numOfRuns,0)));
	
	
	
	for(int run=0;run<numOfRuns;run++){
		cout<<run<<endl;
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
	
	ofstream outFile(outputFolder+"pVavData.txt");
	for(int run=0;run<numOfRuns;run++){
		for(int time=0;time<(int)reportTimes.size()-1;time++){
			for(int species=0;species<(int)specNum.size();species++){
				outFile<<outData[time][species][run]<<" ";
			}
		}
		outFile<<endl;
	}
	outFile.close();
	
	
	return 0;
}

