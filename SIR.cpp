#include <map>
#include <iostream>
#include <tuple>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
using namespace std;


#include "SIR.h"

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

Particle::Particle(int numOfParameters, vector<double> initBounds, vector<double (*)(Particle*,vector<double>&)> initFunctions, tuple<double,double,double,double,double> initParameters){
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

}

Particle::~Particle(){
}


void Particle::dumpParticleDetails(string rootFolder, string id){
}

double sgn(double in){
    return (in<0) ? -1 : ((in>0)? 1 : 0);
}



double Particle::performUpdate(boost::mt19937* inRand, vector<double>& globalBest, FuzzyTree* fuzzyStruct){
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

void Particle::unwrapParameters(){
    beta=currentSolution[0];
    delta=currentSolution[1];
    c=currentSolution[2];
    p=currentSolution[3];
	gamma=currentSolution[4];
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

void rungeKuttaUpdate(Particle* currentParticle, vector<double>& speciesVec, double currentTime, double stoppingTime, double deltaT){
    int numSpecies=speciesVec.size();
    (*currentParticle).unwrapParameters();
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

double fitnessFunction(vector<vector<double> >& trueIn, vector<vector<double> >& testIn){
    double interHold(0);
    for(int i=0;i<(int)trueIn.size();i++){
        for(int j=0;j<(int)trueIn[i].size();j++){
            interHold+=pow(trueIn[i][j]-testIn[i][j],2);
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

vector<double> transformInit(vector<double> inRand, vector<vector<double> >& inCov, vector<double>& inMean){
    vector<double> outRand(inRand.size(),0);


    for(int i=0;i<(int)inRand.size();i++){
        for(int j=0;j<(int)inRand.size();j++){
            outRand[i]+=inRand[j]*inCov[i][j];
        }
        outRand[i]+=inMean[i];
    }

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