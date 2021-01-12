#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <functional>
#include <chrono>
#include <gsl/gsl_linalg.h>


#include <boost/random/mersenne_twister.hpp>
	boost::mt19937 generator;
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "Development.h"



using namespace std;



double calculateMean(vector<double>& inVector);
double calculateMoment(vector<double>& inVector, double order);
vector<vector<double> > invertMatrix(vector<vector<double> >& inData);
double findMaxElement(vector<vector<double> >& inMatrix);
vector<double> generateMomentsVector(vector<vector<double> >& inData, int numSpecies);

int main(){

	generator.seed(time(NULL));

	vector<double> trueMeans={2.,10,50,200};
	vector<double> trueSds={0.01,2,25,100};
	vector<double> testMeans=trueMeans;
	vector<double> testSds=trueSds;

	const int N(500);
	int numSpecies(trueMeans.size());

	vector<boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > > trueGenerators;
	vector<boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > > testGenerators;

	vector<boost::mt19937> noiseEngines(numSpecies*2);

	for(int i=0;i<numSpecies;i++){
		noiseEngines[i].seed(i+time(NULL));
		boost::normal_distribution<double> currentDist(trueMeans[i],trueSds[i]);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > createdEngine(noiseEngines[i],currentDist);
		trueGenerators.push_back(createdEngine);
	}

	for(int i=0;i<numSpecies;i++){
		noiseEngines[i+numSpecies].seed(i+time(NULL));
		boost::normal_distribution<double> currentDist(testMeans[i],testSds[i]);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > createdEngine(noiseEngines[i+numSpecies],currentDist);
		testGenerators.push_back(createdEngine);
	}

	vector<vector<double> > trueData(numSpecies,vector<double> (N,0));
	vector<vector<double> > testData(numSpecies,vector<double> (N,0));

	for(int i=0;i<numSpecies;i++){
		for(int j=0;j<N;j++){
			trueData[i][j]=trueGenerators[i]();
			testData[i][j]=testGenerators[i]();
		}
	}

	double trueMean(calculateMoment(trueData[0],1));
	double trueSecond(calculateMoment(trueData[0],2));


	const int momentVectorSize(2*numSpecies+(numSpecies-1)*numSpecies/2);
	vector<double> trueMoments=generateMomentsVector(trueData,numSpecies);

	vector<vector<double> > wMat(momentVectorSize,vector<double>(momentVectorSize,0));
	for(int k=0;k<N;k++){
		vector<double> interVec(momentVectorSize,0);
		int fillIndex(0);
		for(int i=0;i<numSpecies;i++){
			interVec[fillIndex]=testData[i][k]-trueMoments[fillIndex];
			fillIndex++;
		}
		for(int i=0;i<numSpecies;i++){
			for(int j=i;j<numSpecies;j++){
				interVec[fillIndex]=testData[i][k]*testData[j][k]-trueMoments[fillIndex];
				fillIndex++;
			}
		}
		for(int i=0;i<(int)wMat.size();i++){
			for(int j=0;j<(int)wMat[i].size();j++){
				wMat[i][j]+=interVec[i]*interVec[j]/(double)N;
			}
		}
	}
	wMat=invertMatrix(wMat);
	const double normalizationConst(1./N);
	for_each(wMat.begin(),wMat.end(),[&normalizationConst](vector<double>& v){transform(v.begin(),v.end(),v.begin(),bind(multiplies<double>(),std::placeholders::_1,normalizationConst));});
	wMat=invertMatrix(wMat);
	double maxValue(1/findMaxElement(wMat));
	for_each(wMat.begin(),wMat.end(),[&maxValue](vector<double>& v){transform(v.begin(),v.end(),v.begin(),bind(multiplies<double>(),std::placeholders::_1,maxValue));});
	for(int i=0;i<(int)wMat.size();i++){
		for(int j=0;j<(int)wMat[i].size();j++){
			cout<<wMat[i][j]<<" ";
		}
		cout<<endl;
	}


	return 0;
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

double calculateMoment(vector<double>& inVector, double order){
	double outValue(0);
	return accumulate(inVector.begin(),inVector.end(),outValue,[&order](double x,double y){return x+pow(y,order);})/(double)inVector.size();
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

