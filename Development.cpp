#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <chrono>
#include <gsl/gsl_linalg.h>


#include <boost/random/mersenne_twister.hpp>
	boost::mt19937 generator;
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "Development.h"



using namespace std;



double calculateMean(vector<double>& inVector);
vector<vector<double> > invertMatrix(vector<vector<double> >& inData);

int main(){

	generator.seed(time(NULL));

	vector<double> trueMeans={2,10,50,200};
	vector<double> trueSds(trueMeans.size(),0);
	vector<double> testMeans={};
	vector<double> testSds(testMeans.size(),0);

	const int N(300);
	int numSpecies(trueMeans.size());

	vector<boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > > trueGenerators;
	vector<boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > > testGenerators;

	vector<boost::mt19937> noiseEngines(numSpecies*2);

	for(int i=0;i<numSpecies;i++){
		noiseEngines[i].seed(i);
		boost::normal_distribution<double> currentDist(trueMeans[i],trueSds[i]);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > createdEngine(noiseEngines[i],currentDist);
		trueGenerators.push_back(createdEngine);
	}

	for(int i=0;i<numSpecies;i++){
		noiseEngines[i+numSpecies].seed(i);
		boost::normal_distribution<double> currentDist(testMeans[i],testSds[i]);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<double> > createdEngine(noiseEngines[i+numSpecies],currentDist);
		testGenerators.push_back(createdEngine);
	}




	vector<vector<double> > inData={{1,0.6,0},{0,1.5,1},{0,1,1}};
	vector<vector<double> > outMatrix=invertMatrix(inData);


	return 0;
}

double calculateMean(vector<double>& inVector){
	return accumulate(inVector.begin(),inVector.end(),0)/(double)inVector.size();
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

