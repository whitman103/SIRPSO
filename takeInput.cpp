#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>

#include "SIR_MPI.h"
#include "pVavInteractions.h"
#include "GillespieFunctions.h"

using namespace std;

int loadPvavInputs(vector<double>& speciesVector, vector<double>& initParameters, vector<double>&stoppingTimes, vector<tuple<double,double> >& bounds, string inFile);
void loadFoundParameters(vector<double>& inParameters, string parameterFile);

int main(int argc, char** argv){
	vector<double> speciesVector(6);
	const int numOfParameters(6);
	vector<double> initParameters={0.008,0.1,1.0,0.013904,0.05,0.07};
	vector<double> stoppingTimes={0,100};
	vector<tuple<double,double> > twoBounds;
	int errorFlag=loadPvavInputs(speciesVector,initParameters,stoppingTimes,twoBounds, "InputFolder//outputOne.txt");
	loadFoundParameters(initParameters,"D:\\Downloads\\11_17_2020\\DataFolder_ODEMeans_0\\outFoundParameters.txt");
	stoppingTimes={0,8,32,64,128,256};

	double timeIncrement(0.002);

	SolStruct solutionStructure;
	solutionStructure.timeIncrement=timeIncrement;
	solutionStructure.stoppingTimes=stoppingTimes;

	double scalingFactor(1);

	vector<double (*)(Particle*, vector<double>&)> interactionFuncts={dSyk,dVav,dSV,dpVav,dSHP1,dSHP1Vav};
	
	Particle trueParticle=Particle(numOfParameters,initParameters,twoBounds,interactionFuncts,scalingFactor);

	vector<vector<double> > outData=generateData(&trueParticle,speciesVector,&solutionStructure,0);

	for(int i=0;i<(int)outData.size();i++){
		for(int j=0;j<(int)outData[i].size();j++){
			cout<<outData[i][j]<<" ";
		}
		cout<<endl;
	}

	

	return 0;
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
		cout<<speciesVector.size()<<endl;
		cout<<"Here"<<endl;
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
		cout<<get<0>(bounds[i])<<" "<<get<1>(bounds[i])<<endl;
	}
	return errorOut;

}

void loadFoundParameters(vector<double>& inParameters, string parameterFile){
	ifstream inFile(parameterFile);
	for(int i=0;i<(int)inParameters.size();i++){
		double doubleHold(0);
		inFile>>doubleHold;
		inParameters[i]=doubleHold;
	}
}