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
#include <gsl/gsl_linalg.h>
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
double calculateMean(vector<double>& inVector);
vector<vector<double> > invertMatrix(vector<vector<double> >& inData);
double findMaxElement(vector<vector<double> >& inMatrix);
vector<double> generateMomentsVector(vector<vector<double> >& inData, int numSpecies);

int main(int argc, char** argv){

    generator.seed(time(NULL));
    /*string localPath=std::filesystem::current_path();
	string outputFolder=localPath+"\\DataFolder";
    std::filesystem::create_directory(outputFolder);*/
	string outputFolder="DataFolder_GMMIdentity_"+string(argv[3]);
	mkdir(outputFolder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    outputFolder+="//";

    string customString("extrinsicNoise_");
    

    bool exNoise(true);
	bool resetParameter(true);
	bool variancesIncluded(true);
	//T, I, V, R
	const int numOfSpecies(6);
	const int momentVectorSize(2*numOfSpecies+(numOfSpecies-1)*numOfSpecies/2);
	
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
	const int numOfIterations(50);
	//Number of Gillespie samples to use for distributions
	const int numOfSamples(2500);

	//Number of Particle sets to run
	const int numOfRuns(100);

	//Generates cholesky matrix to produce lognormal distributions
	vector<vector<double> > inValues;
	vector<double> means;
	loadCovariance(means, inValues,"outCov");
	//vector<vector<double> > outDecomp=generateCholesky(inValues);
	vector<double> parameterResetValues(initParameters.size(),0);
	if(resetParameter){
		ifstream readIn("InputFolder//start"+string(argv[3])+".txt");
		for(int i=0;i<(int)initParameters.size();i++){
			readIn>>parameterResetValues[i];
		}
		readIn.close();
	}

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

	int mahalanDimension(2*numOfSpecies+(numOfSpecies/2)*(numOfSpecies-1));
	vector<vector<double> > mahalanMetric(mahalanDimension,vector<double>(mahalanDimension,0));
	double sendMetric[mahalanDimension*mahalanDimension];
	for(int i=0;i<(int)mahalanMetric.size();i++){
		mahalanMetric[i][i]=1.;
	}
	for(int i=0;i<pow(mahalanDimension,2);i++){
		sendMetric[i]=0;
	}
    
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
	vector<vector<vector<double> > > trueFullDistribution(numOfSamples);

	switch(solutionStyle){
		case 0:
		{
			for(int sample=0;sample<numOfSamples;sample++){
				speciesVector=resetSpecies;
				speciesVector[0]=extrinsicNoiseGenerators[0]();
				speciesVector[1]=extrinsicNoiseGenerators[1]();
				speciesVector[4]=extrinsicNoiseGenerators[2]();
				trueFullDistribution[sample]=generateData(&trueParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
			}
		}
		break;
		default:
		return 0;
	}
	vector<double> trueMahalan;
	if(variancesIncluded){
		vector<vector<double> > firstDataSet(swapSampleIndices(trueFullDistribution));
		trueMahalan=generateMahalanVector(firstDataSet);
	}

	

	

	
	
	double inDelta(10);
	vector<FuzzyTree> FuzzyStructure(numOfParticles,FuzzyTree(inDelta));

	ofstream outRunge;
	
	outRunge.open(outputFolder+customString+"outRunge_testNoise.txt");
	for(int i=0;i<(int)trueParticle.currentSolution.size();i++){
		outRunge<<trueParticle.currentSolution[i]<<" ";
	}
	outRunge<<endl;
	
	
	
	int sizeOfParameterVector(initParameters.size());
	double fitnessCollection[numOfParticles];
	double parameterPassVector[sizeOfParameterVector];
	double parameterMatrixHold[sizeOfParameterVector*numOfParticles];
	double metricArray[mahalanDimension*mahalanDimension];
	int bestParticleIndex(-1);
	
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
		if(!resetParameter){
			for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
				auto[lowerBound,upperBound]=twoBounds[i];
				threadParticle.currentSolution[i]=pow(10,randPull()*log10(upperBound/lowerBound))*lowerBound;
			}
		}
		else{
			for(int i=0;i<(int)threadParticle.currentSolution.size();i++){
				double lowerBound(parameterResetValues[i]/2.), upperBound(parameterResetValues[i]*2.);
				threadParticle.currentSolution[i]=pow(10,randPull()*log10(upperBound/lowerBound))*lowerBound;
			}
		}
		
		
        vector<vector<vector<int> > > testDistributions;
		switch(solutionStyle){
			case 0:
			{
				//Initialize
				vector<vector<double> > testArray(stoppingTimes.size()-1,vector<double>(numOfSpecies,0));
				vector<vector<vector<double> > > testFullDistributions(numOfSamples);
				vector<vector<double> > resetTestArray=testArray;
				speciesVector=resetSpecies;
				for(int sample=0;sample<numOfSamples;sample++){
					speciesVector=resetSpecies;
					speciesVector[0]=extrinsicNoiseGenerators[0]();
					speciesVector[1]=extrinsicNoiseGenerators[1]();
					speciesVector[4]=extrinsicNoiseGenerators[2]();
					testFullDistributions[sample]=generateData(&threadParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
				}

				vector<vector<double> > firstDataSet(swapSampleIndices(testFullDistributions));
				threadParticle.currentFitness=mahalanFitness(trueMahalan,firstDataSet,mahalanMetric);
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
					for(int sample=0;sample<numOfSamples;sample++){
						speciesVector=resetSpecies;
						speciesVector[0]=extrinsicNoiseGenerators[0]();
						speciesVector[1]=extrinsicNoiseGenerators[1]();
						speciesVector[4]=extrinsicNoiseGenerators[2]();
						testFullDistributions[sample]=generateData(&threadParticle,speciesVector,&solutionStructure,styleMap["RungeKutta"]);
					}
					
					vector<vector<double> > firstDataSet(swapSampleIndices(testFullDistributions));
					threadParticle.currentFitness=fabs(mahalanFitness(trueMahalan,firstDataSet,mahalanMetric));
					
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
						bestParticleIndex=checkForNewGlobalBest(fitnessCollection, parameterMatrixHold, parameterPassVector, numOfParticles, globalBestFitness,numOfParameters);
					}
					MPI_Bcast(parameterPassVector, sizeOfParameterVector, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					threadParticle.twoBoundPerformUpdate(&generator,parameterPassVector,&fuzzyStruct);
					if(iteration==numOfIterations-1&&run==(int)numOfRuns/2){
						if(taskID==bestParticleIndex){
							vector<vector<double> > wMat(momentVectorSize,vector<double>(momentVectorSize,0));
							vector<vector<double> > firstDataSet(swapSampleIndices(testFullDistributions));
							for(int k=0;k<numOfSamples;k++){
								vector<double> interVec(momentVectorSize,0);
								int fillIndex(0);
								for(int i=0;i<numOfSpecies;i++){
									interVec[fillIndex]=firstDataSet[i][k]-trueMahalan[fillIndex];
									fillIndex++;
								}
								for(int i=0;i<numOfSpecies;i++){
									for(int j=i;j<numOfSpecies;j++){
										interVec[fillIndex]=firstDataSet[i][k]*firstDataSet[j][k]-trueMahalan[fillIndex];
										fillIndex++;
									}
								}
								for(int i=0;i<(int)wMat.size();i++){
									for(int j=0;j<(int)wMat[i].size();j++){
										wMat[i][j]+=interVec[i]*interVec[j]/(double)numOfSamples;
									}
								}
							}
							const double normalizationConst(1./numOfSamples);
							for_each(wMat.begin(),wMat.end(),[&normalizationConst](vector<double>& v){transform(v.begin(),v.end(),v.begin(),bind(multiplies<double>(),std::placeholders::_1,normalizationConst));});
							wMat=invertMatrix(wMat);
							double maxValue(1/findMaxElement(wMat));
							for_each(wMat.begin(),wMat.end(),[&maxValue](vector<double>& v){transform(v.begin(),v.end(),v.begin(),bind(multiplies<double>(),std::placeholders::_1,maxValue));});
							mahalanMetric=wMat;
							for(int i=0;i<mahalanDimension;i++){
								for(int j=0;j<mahalanDimension;j++){
									metricArray[i*mahalanDimension+j]=mahalanMetric[i][j];
								}
							}
							MPI_Send(metricArray,mahalanDimension*mahalanDimension,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
						}
						if(taskID==0){
							MPI_Recv(metricArray,mahalanDimension*mahalanDimension, MPI_DOUBLE, bestParticleIndex,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
						}
						MPI_Bcast(metricArray,mahalanDimension*mahalanDimension,MPI_DOUBLE,0,MPI_COMM_WORLD);
						for(int i=0;i<mahalanDimension;i++){
							for(int j=0;j<mahalanDimension;j++){
								mahalanMetric[i][j]=metricArray[i*mahalanDimension+j];
							}
						}
					}
				}
			}
			break;
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

vector<double> generateMomentsVector(vector<vector<double> >& inData, int numSpecies){
	int fillIndex(0);
	vector<double> trueMoments(2*numSpecies+(numSpecies-1)*numSpecies/2);
	for(int i=0;i<numSpecies;i++){
		trueMoments[fillIndex]=calculateMoment(inData[i],1);
		fillIndex++;
	}
	for(int i=0;i<numSpecies;i++){
		for(int j=i;j<numSpecies;j++){
			double sumHold(0);
			for(int k=0;k<(int)inData[i].size();k++){
				sumHold+=inData[i][k]*inData[j][k];
			}
			trueMoments[fillIndex]=sumHold/(double)inData[i].size();
			fillIndex++;
		}
	}
	return trueMoments;
}

double calculateMean(vector<double>& inVector){
	double outValue(0);
	return accumulate(inVector.begin(),inVector.end(),outValue)/(double)inVector.size();
}

double findMaxElement(vector<vector<double> >& inMatrix){
	double outValue(0);
	for(int i=0;i<(int)inMatrix.size();i++){
		auto iter=max_element(inMatrix[i].begin(),inMatrix[i].end());
		if(outValue<*iter){
			outValue=*iter;
		}
	}
	return outValue;
}
// Data input looks like timeIndex, speciesIndex, runIndex

vector<vector<double> > invertMatrix(vector<vector<double> >& inMatrix){
	const int N(inMatrix.size());

	vector<vector<double> > outMatrix=inMatrix;
	double inData[(int)pow(N,2)];
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			inData[i*inMatrix.size()+j]=inMatrix[i][j];
		}
	}
	double inverseMat[(int)pow(N,2)];
	int s;
	gsl_matrix_view m= gsl_matrix_view_array(inData,N,N);
	gsl_matrix_view inv=gsl_matrix_view_array(inverseMat,N,N);
	gsl_permutation *p = gsl_permutation_alloc(N);

	gsl_linalg_LU_decomp(&m.matrix,p,&s);
	gsl_linalg_LU_invert(&m.matrix,p,&inv.matrix);
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			outMatrix[i][j]=gsl_matrix_get(&inv.matrix,i,j);
		}
	}
	return outMatrix;
}








