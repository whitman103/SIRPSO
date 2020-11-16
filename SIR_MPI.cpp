#include <map>
#include <iostream>
#include <tuple>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include <boost/random/mersenne_twister.hpp>
using namespace std;


#include "SIR_MPI.h"

//Values taken directly from Nobile2017 paper
FuzzyTree::FuzzyTree(double inDelta){
	inertiaMap["Low"]=0.3;
	inertiaMap["Medium"]=0.5;
	inertiaMap["High"]=1.0;
	socialMap["Low"]=1.0;
	socialMap["Medium"]=2.0;
	socialMap["High"]=3.0;
	cognitiveMap["Low"]=0.1;
	cognitiveMap["Medium"]=1.5;
	cognitiveMap["High"]=3.0;
	LMap["Low"]=0.0;
	LMap["Medium"]=0.001;
	LMap["High"]=0.01;
	UMap["Low"]=0.1;
	UMap["Medium"]=0.15;
	UMap["High"]=0.2;
	linguisticMap[0]="Low";
	linguisticMap[1]="Medium";
	linguisticMap[2]="High";
	
	deltaMax=inDelta;
	delta1=0.2*deltaMax;
	delta2=0.4*deltaMax;
	delta3=0.6*deltaMax;
	
	inertia=inertiaMap["Medium"];
	social=socialMap["Medium"];
	cognitive=cognitiveMap["Medium"];
	L=LMap["Medium"];
	U=UMap["Medium"];
};


FuzzyTree::~FuzzyTree(){
};

void FuzzyTree::setParameters(){
	
	calculateDeltaMembershipValues();
	calculatePhiMembershipValues();
	setInertia();
	setSocial();
	setCognitive();
	setL();
	setU();
}

void FuzzyTree::calculateDeltaMembershipValues(){
	if(delta<delta1){
		deltaMembershipValues={1,0,0};
	} 
	else{
		if(delta>delta1&&delta<delta2){
			deltaMembershipValues={(delta2-delta)/(delta2-delta1),(delta-delta1)/(delta2-delta1),0};
		}
		else{
			if(delta>delta2&&delta<delta3){
				deltaMembershipValues={0,(delta3-delta)/(delta3-delta2),(delta-delta2)/(delta3-delta2)};
			}
			else{
				if(delta>delta3){
					deltaMembershipValues={0,0,1};
				}
			}
		}
	}
}

void FuzzyTree::calculatePhiMembershipValues(){
	if(phi<0){
		phiMembershipValues={-phi,1-fabs(phi),0};
	}
	else{
		phiMembershipValues={0,1-fabs(phi),phi};
	}
}

void FuzzyTree::setInertia(){
	inertia=0;
	inertia+=inertiaMap["Low"]*(phiMembershipValues[2]+deltaMembershipValues[0]);
	inertia+=inertiaMap["Medium"]*(phiMembershipValues[1]+deltaMembershipValues[1]);
	inertia+=inertiaMap["High"]*(phiMembershipValues[0]+deltaMembershipValues[2]);
	inertia/=2.;
}

void FuzzyTree::setSocial(){
	social=0;
	social+=socialMap["Low"]*(phiMembershipValues[0]+deltaMembershipValues[1]);
	social+=socialMap["Medium"]*(phiMembershipValues[1]+deltaMembershipValues[0]);
	social+=socialMap["High"]*(phiMembershipValues[2]+deltaMembershipValues[2]);
	social/=2.;
}

void FuzzyTree::setCognitive(){
	cognitive=0;
	cognitive+=cognitiveMap["Low"]*deltaMembershipValues[2];
	cognitive+=cognitiveMap["Medium"]*(phiMembershipValues[2]+phiMembershipValues[1]+deltaMembershipValues[0]+deltaMembershipValues[1]);
	cognitive+=cognitiveMap["High"]*(phiMembershipValues[0]);
	cognitive/=2.;
}

void FuzzyTree::setL(){
	L=0;
	L+=LMap["Low"]*(phiMembershipValues[1]+phiMembershipValues[0]+deltaMembershipValues[2]);
	L+=LMap["Medium"]*(deltaMembershipValues[0]+deltaMembershipValues[1]);
	L+=LMap["High"]*(phiMembershipValues[2]);
	L/=2.;
}

void FuzzyTree::setU(){
	U=0;
	U+=UMap["Low"]*(deltaMembershipValues[0]);
	U+=UMap["Medium"]*(phiMembershipValues[0]+phiMembershipValues[1]+deltaMembershipValues[1]);
	U+=UMap["High"]*(phiMembershipValues[2]+deltaMembershipValues[2]);
	U/=2.;
}

void FuzzyTree::calculatePhi(double lastFitness, double currentFitness){
	phi=0;
	phi=delta/deltaMax;
	phi*=(min(currentFitness,phiNormalization)-min(lastFitness,phiNormalization))/phiNormalization;
}

Particle::Particle(int numOfParameters, vector<double> initBounds, vector<double (*)(Particle*,vector<double>&)> initFunctions, tuple<double,double,double,double,double> initParameters, int scalingSize){
    currentSolution.resize(numOfParameters);
    bestSolution.resize(numOfParameters);
    currentVelocity.resize(numOfParameters);
    interactionFunctions=initFunctions;
	bounds=initBounds;
    bestFitness=0;
    currentFitness=0;
	beta=get<0>(initParameters);
	currentSolution[0]=beta;
	delta=get<1>(initParameters);
	currentSolution[1]=delta;
	c=get<2>(initParameters);
	currentSolution[2]=c;
	p=get<3>(initParameters);
	currentSolution[3]=p;
	gamma=get<4>(initParameters);
	currentSolution[4]=gamma;
	scalingFactor=scalingSize;
}

Particle::Particle(int numOfParameters, vector<double> Parameters, vector<double (*)(Particle*,vector<double>&)> initFunctions, int scalingSize){
	currentSolution.resize(numOfParameters);
    bestSolution.resize(numOfParameters);
    currentVelocity.resize(numOfParameters);
    interactionFunctions=initFunctions;
    bestFitness=0;
    currentFitness=0;
	currentSolution=Parameters;
	scalingFactor=scalingSize;
}

Particle::Particle(int numOfParameters, vector<double> Parameters, vector<double> initBounds, vector<double (*)(Particle*,vector<double>&)> initFunctions, int scalingSize){
	currentSolution.resize(numOfParameters);
    bestSolution.resize(numOfParameters);
    currentVelocity.resize(numOfParameters);
    interactionFunctions=initFunctions;
    bestFitness=0;
    currentFitness=0;
	currentSolution=Parameters;
	scalingFactor=scalingSize;
	bounds=initBounds;
}

Particle::Particle(int numOfParameters, vector<double> Parameters, vector<tuple<double,double> > initBounds, vector<double(*)(Particle*,vector<double>&)> initFunctions, int scalingSize){
	currentSolution.resize(numOfParameters);
    bestSolution.resize(numOfParameters);
    currentVelocity.resize(numOfParameters);
    interactionFunctions=initFunctions;
    bestFitness=0;
    currentFitness=0;
	currentSolution=Parameters;
	scalingFactor=scalingSize;
	twoBounds=initBounds;
}

Particle::~Particle(){
}

void Particle::dumpParticleDetails(ofstream* outStream){
	for(int i=0;i<(int)currentSolution.size();i++){
		(*outStream)<<currentSolution[i]<<" ";
	}
	(*outStream)<<endl;
}

double sgn(double in){
    return (in<0) ? -1 : ((in>0)? 1 : 0);
}


double Particle::performUpdate(boost::mt19937* inRand, double* globalBest, FuzzyTree* fuzzyStruct){
    for(int i=0;i<(int)currentSolution.size();i++){
        double proposedUpdate(0);
        double rand1((double)(*inRand)()/(double)(*inRand).max());
        double rand2((double)(*inRand)()/(double)(*inRand).max());
        proposedUpdate+=(*fuzzyStruct).inertia*currentVelocity[i];
        proposedUpdate+=(*fuzzyStruct).social*rand1*(globalBest[i]-currentSolution[i]);
        proposedUpdate+=(*fuzzyStruct).cognitive*rand2*(bestSolution[i]-currentSolution[i]);
        if(fabs(proposedUpdate)>(*fuzzyStruct).U*bounds[i]){
            proposedUpdate=(*fuzzyStruct).U*bounds[i]*sgn(proposedUpdate);
        }

        if((currentSolution[i]+proposedUpdate)<0){
            currentSolution[i]=1./100.*bounds[i];
            currentVelocity[i]=-1.*proposedUpdate*1./10.;
        } else{
            if((currentSolution[i]+proposedUpdate)>bounds[i]){
                currentSolution[i]=bounds[i];
                currentVelocity[i]=-1.*proposedUpdate*1./10.;
            } else{
                currentSolution[i]+=proposedUpdate;
                currentVelocity[i]=proposedUpdate;
            }
        }
    }
	return 0;
}

double Particle::twoBoundPerformUpdate(boost::mt19937* inRand, double* globalBest, FuzzyTree* fuzzyStruct){
	for(int i=0;i<(int)currentSolution.size();i++){
		auto[lowerBound,upperBound]=twoBounds[i];
        double proposedUpdate(0);
        double rand1((double)(*inRand)()/(double)(*inRand).max());
        double rand2((double)(*inRand)()/(double)(*inRand).max());
        proposedUpdate+=(*fuzzyStruct).inertia*currentVelocity[i];
        proposedUpdate+=(*fuzzyStruct).social*rand1*(globalBest[i]-currentSolution[i]);
        proposedUpdate+=(*fuzzyStruct).cognitive*rand2*(bestSolution[i]-currentSolution[i]);
        if(fabs(proposedUpdate)>(*fuzzyStruct).U*upperBound){
            proposedUpdate=(*fuzzyStruct).U*upperBound*sgn(proposedUpdate);
        }

        if((currentSolution[i]+proposedUpdate)<lowerBound){
            currentSolution[i]=lowerBound;
            currentVelocity[i]=-1.*proposedUpdate*1./10.;
        } else{
            if((currentSolution[i]+proposedUpdate)>upperBound){
                currentSolution[i]=upperBound;
                currentVelocity[i]=-1.*proposedUpdate*1./10.;
            } else{
                currentSolution[i]+=proposedUpdate;
                currentVelocity[i]=proposedUpdate;
            }
        }
    }
	return 0;
}

void Particle::unwrapParameters(){
    beta=currentSolution[0];
    delta=currentSolution[1];
    c=currentSolution[2];
    p=currentSolution[3];
	gamma=currentSolution[4];
}

void Particle::unwrap_pVavParameters(){
	k0=currentSolution[0];
	k1=currentSolution[1];
	k2=currentSolution[2];
	k3=currentSolution[3];
	k4=currentSolution[4];
	k5=currentSolution[5];
}
	
	
//Species are currently T, I, V, R

double firstInteraction(Particle* currentParticle, vector<double>& species){
    return -1.*(*currentParticle).beta*species[0]*species[2]+(*currentParticle).gamma*species[3];
}

double secondInteraction(Particle* currentParticle, vector<double>& species){
    return (*currentParticle).beta*species[0]*species[2]-(*currentParticle).c*species[1];
}

double thirdInteraction(Particle* currentParticle, vector<double>& species){
    return (*currentParticle).p*species[1]-(*currentParticle).delta*species[2];
}

double fourthInteraction(Particle* currentParticle, vector<double>& species){
	return (*currentParticle).c*species[1]-(*currentParticle).gamma*species[3];
}

vector<double> Particle::convertFromParticleToGillespie(){
	unwrapParameters();
	vector<double> outConsts={beta/(double)scalingFactor,gamma,c,p,delta};
	return outConsts;
}

vector<double> Particle::pVavConvertParticleGillespie(){
	unwrap_pVavParameters();
	vector<double> outConsts={k0,k1,k2,k3,k4,k5};
	return outConsts;
}

void Particle::divideBeta(){
	beta=beta/(double)scalingFactor;
}

void Particle::multiplyBeta(){
	beta=beta*(double)scalingFactor;
}

void rungeKuttaUpdate(Particle* currentParticle, vector<double>& speciesVec, double currentTime, double stoppingTime, double deltaT){
    int numSpecies=speciesVec.size();
    (*currentParticle).unwrap_pVavParameters();
    int n=(int)(((stoppingTime-currentTime))/(deltaT));
    vector<double (*)(Particle*,vector<double>&)> interactionPointer=(*currentParticle).interactionFunctions;

    for(int t=0;t<n;t++){
        vector<double> interSpecies(numSpecies,0);
        vector<double> k1(numSpecies,0);
        for(int i=0;i<(int)k1.size();i++){
            k1[i]=deltaT*interactionPointer[i](currentParticle,speciesVec);
			interSpecies[i]=speciesVec[i]+k1[i]/2.;
        }
        vector<double> k2(numSpecies,0);
        for(int i=0;i<(int)k2.size();i++){
            k2[i]=deltaT*interactionPointer[i](currentParticle,interSpecies);
			interSpecies[i]=speciesVec[i]+k2[i]/2.;
        }
        vector<double> k3(numSpecies,0);
        for(int i=0;i<(int)k3.size();i++){
            k3[i]=deltaT*interactionPointer[i](currentParticle,interSpecies);
			interSpecies[i]=speciesVec[i]+k3[i];
        }
        vector<double> k4(numSpecies,0);
        for(int i=0;i<(int)k4.size();i++){
            k4[i]=deltaT*interactionPointer[i](currentParticle,interSpecies);
        }
        for(int i=0;i<numSpecies;i++){
            speciesVec[i]+=(k1[i]/6.+k2[i]/3.+k3[i]/3.+k4[i]/6.);
        }
    }
}

double fitnessFunction(vector<vector<double> >& trueMean, vector<vector<double> >& testMean){
    double interHold(0);
    for(int i=0;i<(int)trueMean.size();i++){
        for(int j=0;j<(int)trueMean[i].size();j++){
            interHold+=pow(trueMean[i][j]-testMean[i][j],2);
        }
    }
    return interHold;
}

double fitnessFunction(vector<vector<double> >& trueMean, vector<vector<double> >& testMean, vector<vector<double> >& trueVar, vector<vector<double> >& testVar){
    double interHold(0);
    for(int i=0;i<(int)trueMean.size();i++){
        for(int j=0;j<(int)trueMean[i].size();j++){
            interHold+=pow(trueMean[i][j]-testMean[i][j],2)+pow(trueVar[i][j]-testVar[i][j],2)/(double)trueMean[i].size();
        }
    }
    return interHold;
}

vector<vector<double> > generateCholesky(vector<vector<double> >& inMatrix){
    vector<vector<double> > choleskyOut(inMatrix.size(),vector<double>(inMatrix[0].size(),0));
    choleskyOut[0][0]=sqrt(inMatrix[0][0]);
    for(int i=0;i<(int)inMatrix.size();i++){
        for(int j=0;j<=i;j++){
            if(i==j&&i!=0){
                double squareHold(0);
                for(int k=0;k<j;k++){
                    squareHold+=pow(choleskyOut[i][k],2);
                }
                choleskyOut[i][i]=sqrt(inMatrix[i][i]-squareHold);
            }
            else{
                double outerProd(0);
                for(int k=0;k<j;k++){
                    outerProd+=choleskyOut[i][k]*choleskyOut[j][k];
                }
            choleskyOut[i][j]=(inMatrix[i][j]-outerProd)/choleskyOut[j][j];
            }
        }
    }
    return choleskyOut;
}


void loadCovariance(vector<double>& outMeans, vector<vector<double> >& inMatrix, string dataPath){
    int inSize(0);
    ifstream inData(dataPath+".txt");
	if(!inData.good()){
		cout<<"failed to load data for covariance matrix"<<endl;
	}
    inData>>inSize;
    outMeans.resize(inSize);
    for(int i=0;i<inSize;i++){
        inData>>outMeans[i];
    }
    vector<vector<double> > interMatrix(inSize,vector<double> (inSize,0));
    for(int i=0;i<(int)interMatrix.size();i++){
        for(int j=0;j<(int)interMatrix[i].size();j++){
            inData>>interMatrix[i][j];
        }
    }
    inData.close();
    inMatrix=interMatrix;
}

vector<double> transformInit(vector<double> inRand, vector<vector<double> >& inCov, vector<double>& inMean, boost::mt19937* generatorIn){
    vector<double> outRand(inRand.size(),0);


    for(int i=0;i<(int)inRand.size();i++){
        for(int j=0;j<(int)inRand.size();j++){
            outRand[i]+=inRand[j]*inCov[i][j];
        }
        outRand[i]+=inMean[i];
    }
	//Explicit changes to the setup 
	outRand[0]=min(0.99, outRand[0]);
	outRand[1]=1-outRand[0];
	outRand[2]=inMean[2]*(1.+2.*(2.*(double)(*generatorIn)()/(double)(*generatorIn).max()-1.));
	outRand[3]=0;

    return outRand;
}


vector<vector<double> > generateData(Particle* inParticle, vector<double>& inSpecies, SolStruct* solChars, int styleFlag){
	switch(styleFlag){
		case 0:
		{
			vector<vector<double> > outData((*solChars).stoppingTimes.size()-1,vector<double>(inSpecies.size(),0));
			double timeIncrement=(*solChars).timeIncrement;
			for(int reportIndex=0;reportIndex<(int)(*solChars).stoppingTimes.size()-1;reportIndex++){
				rungeKuttaUpdate(inParticle,inSpecies,(*solChars).stoppingTimes[reportIndex],(*solChars).stoppingTimes[reportIndex+1],timeIncrement);
				for(int i=0;i<(int)inSpecies.size();i++){
					outData[reportIndex][i]=inSpecies[i];
				}
			}
			return outData;
		}
		default:
		vector<vector<double> > voidCast;
		return voidCast;
	}
}

void checkForNewGlobalBest(double* fitnessCollection, double* parameterMatrixHold, double* parameterPassVector, int numOfParticles, double& globalFitness, int numOfParameters){
	int bestParticle(0);
	double bestFitness(1e13);
	for(int i=0;i<numOfParticles;i++){
		if(fitnessCollection[i]<bestFitness){
			bestFitness=fitnessCollection[i];
			bestParticle=i;
		}
	}
	for(int i=0;i<numOfParameters;i++){
		parameterPassVector[i]=parameterMatrixHold[bestParticle*numOfParameters+i];
	}
	globalFitness=bestFitness;
}
	
vector<vector<vector<double> > > performGillespieSimulation(Particle* inParticle, Gillespie* inReactionObject, vector<double>& reportTimes, vector<int>& specNum, int numOfRuns){
	vector<vector<vector<double> > > outData(reportTimes.size()-1,vector<vector<double> > (specNum.size(), vector<double> (numOfRuns,0)));
	(*inReactionObject).reactConsts=(*inParticle).convertFromParticleToGillespie();
	for(int run=0;run<numOfRuns;run++){
		specNum=(*inReactionObject).resetSpecies;
		double currentTime(0);
		int reportIndex(0);
		do{
			tuple<int,double> hold=(*inReactionObject).PerformTimeStep2(specNum);
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
				(*inReactionObject).specChange(specNum,get<0>(hold),(*inReactionObject).changeCoeffs);
				if(currentTime>reportTimes[reportIndex+1]){
					for(int i=0;i<(int)specNum.size();i++){
						outData[reportIndex][i][run]=specNum[i];
					}
					reportIndex++;
				}
			}
		}while(reportIndex<(int)reportTimes.size()-1);
	}
	
	return outData;
}
	
vector<vector<vector<double> > > performGillespieSimulation(Particle* inParticle, Gillespie* inReactionObject, vector<double>& reportTimes, vector<int>& specNum, int numOfRuns, boost::mt19937* inGenerator){
	vector<vector<vector<double> > > outData(reportTimes.size()-1,vector<vector<double> > (specNum.size(), vector<double> (numOfRuns,0)));
	(*inReactionObject).reactConsts=(*inParticle).convertFromParticleToGillespie();
	for(int run=0;run<numOfRuns;run++){
		specNum=(*inReactionObject).resetSpecies;
		int totalCells(specNum[0]+specNum[1]);
		vector<double> doubleSpecies(specNum.begin(),specNum.end());
		doubleSpecies[1]=doubleSpecies[1]*(1.+.2*(2*(double)(*inGenerator)()/(double)(*inGenerator).max()-1.))+1;
		doubleSpecies[0]=totalCells-doubleSpecies[1];
		doubleSpecies[2]=doubleSpecies[2]*(1.+.3*(2*(double)(*inGenerator)()/(double)(*inGenerator).max()-1.));
		vector<int> backToInt(doubleSpecies.begin(),doubleSpecies.end());
		specNum=backToInt;
		double currentTime(0);
		int reportIndex(0);
		do{
			tuple<int,double> hold=(*inReactionObject).PerformTimeStep2(specNum);
			currentTime+=get<1>(hold);
			
			
			if(get<1>(hold)<0){
				for(int index=reportIndex;index<(int)reportTimes.size()-1;index++){
					for(int i=0;i<(int)specNum.size();i++){
						outData[index][i][run]=specNum[i];
					}
				}
				reportIndex=(int)reportTimes.size();
			}
			else{
				(*inReactionObject).specChange(specNum,get<0>(hold),(*inReactionObject).changeCoeffs);
				if(currentTime>reportTimes[reportIndex+1]){
					for(int i=0;i<(int)specNum.size();i++){
						outData[reportIndex][i][run]=specNum[i];
					}
					reportIndex++;
				}
			}
		}while(reportIndex<(int)reportTimes.size()-1);
	}
	
	return outData;
}

tuple<vector<vector<double> >, vector<vector<double> > > calculateMeansAndVar(vector<vector<vector<double> > >& inDist){
	vector<vector<double> > outMeans(inDist.size(), vector<double> (inDist[0].size(),0));
	vector<vector<double> > outVar=outMeans;
	double interHold(0);
	for(int i=0;i<(int)inDist.size();i++){
		for(int j=0;j<(int)inDist[i].size();j++){
			interHold=0;
			for(int k=0;k<(int)inDist[i][j].size();k++){
				interHold+=inDist[i][j][k]/(double)inDist[i][j].size();
			}
			outMeans[i][j]=interHold;
			interHold=0;
			for(int k=0;k<(int)inDist[i][j].size();k++){
				interHold+=pow(inDist[i][j][k]-outMeans[i][j],2);
			}
			outVar[i][j]=sqrt(interHold)/(double)inDist[i][j].size();
		}
	}
	
	
	return make_tuple(outMeans, outVar);
}

tuple<vector<vector<double> >, vector<vector<double> > > calculateMeansAndVar(vector<vector<vector<int> > >&inDist){
	vector<vector<double> > outMeans(inDist.size(), vector<double> (inDist[0].size(),0));
	vector<vector<double> > outVar=outMeans;
	double interHold(0);
	for(int i=0;i<(int)inDist.size();i++){
		for(int j=0;j<(int)inDist[i].size();j++){
			interHold=0;
			for(int k=0;k<(int)inDist[i][j].size();k++){
				interHold+=inDist[i][j][k]/(double)inDist[i][j].size();
			}
			outMeans[i][j]=interHold;
			interHold=0;
			for(int k=0;k<(int)inDist[i][j].size();k++){
				interHold+=pow(inDist[i][j][k]-outMeans[i][j],2);
			}
			outVar[i][j]=sqrt(interHold)/(double)inDist[i][j].size();
		}
	}
	
	
	return make_tuple(outMeans, outVar);

}

tuple<vector<vector<double> >, vector<vector<double> > > generateGillespieData(Particle* inParticle, Gillespie* inReactionObject, vector<double>& reportTimes, vector<int>& specNum, int numOfRuns){
	vector<vector<vector<double> > > inData=performGillespieSimulation(inParticle,inReactionObject,reportTimes,specNum,numOfRuns);
	
	return calculateMeansAndVar(inData);
}

tuple<vector<vector<double> >, vector<vector<double> > > generateGillespieData(Particle* inParticle, Gillespie* inReactionObject, vector<double>& reportTimes, vector<int>& specNum, int numOfRuns, boost::mt19937* inGenerator){
	vector<vector<vector<double> > > inData=performGillespieSimulation(inParticle,inReactionObject,reportTimes,specNum,numOfRuns,inGenerator);
	
	return calculateMeansAndVar(inData);
}

vector<vector<vector<double> > > generateDistributions(Particle* inParticle, Gillespie* inReactionObject, vector<double>& reportTimes, vector<int>& specNum, int numOfRuns, boost::mt19937* inGenerator, string outFile){
	vector<vector<vector<double> > > inData=performGillespieSimulation(inParticle, inReactionObject, reportTimes, specNum, numOfRuns, inGenerator);
	ofstream outStream(outFile);
	for(int i=1;i<(int)reportTimes.size();i++){
		outStream<<reportTimes[i]<<",";
		for(int j=0;j<specNum.size()-1;j++){
			outStream<<reportTimes[i]<<",";
		}
	}
	outStream<<endl;
	
	for(int i=0;i<numOfRuns;i++){
		for(int j=0;j<(int)reportTimes.size()-1;j++){
			for(int k=0;k<(int)specNum.size();k++){
				outStream<<inData[j][k][i]<<",";
			}
		}
		outStream<<endl;
	}
	outStream.close();
	return inData;
}

vector<double> readVectorFile(string inString){
	ifstream inFile(inString);
	int inHold;
	inFile>>inHold;
	int size(inHold);
	vector<double> inVec(size,0);
	for(int i=0;i<size;i++){
		inFile>>inVec[i];
		cout<<inVec[i]<<" ";
	}
	inFile.close();
	return inVec;
}

tuple<double,double,double,double,double> readParameterData(string inFile){
	tuple<double,double,double,double,double> parameterSet;
	ifstream inData(inFile);
	inData>>get<0>(parameterSet);
	inData>>get<1>(parameterSet);
	inData>>get<2>(parameterSet);
	inData>>get<3>(parameterSet);
	inData>>get<4>(parameterSet);
	
	return parameterSet;
}

