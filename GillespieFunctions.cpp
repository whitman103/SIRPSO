#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <tuple>
#include <chrono>
#include <random>
#include "GillespieFunctions.h"

inline double tau(double a0, double rand){
	return 1./a0*log(1./rand);
};

Gillespie::Gillespie(string Filename){
	initializeReactions(Filename,changeCoeffs,propCoeffs);
	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed1);
	
};

Gillespie::~Gillespie(){
	changeCoeffs.clear();
	propCoeffs.clear();
}

tuple<int,double> Gillespie::PerformTimeStep(vector<int>& specNum){
	tuple<int,double> output={0,0};
	double r1=(double) generator()/(double)generator.max();
	double r2=(double) generator()/(double)generator.max();
	double a0=sumupslow(specNum,reactConsts,propCoeffs,propVector);
	get<0>(output)=r1choice(a0,r1,propVector);
	get<1>(output)=tau(a0,r2);
	specChange(specNum,get<0>(output),changeCoeffs);
	return output;
}

tuple<int,double> Gillespie::PerformTimeStep2(vector<int>& specNum){
	tuple<int,double> output={0,0};
	double r1=(double) generator()/(double)generator.max();
	double r2=(double) generator()/(double)generator.max();
	double a0=sumupslow(specNum,reactConsts,propCoeffs,propVector);
	get<0>(output)=r1choice(a0,r1,propVector);
	get<1>(output)=tau(a0,r2);
	return output;
}
	
	
	
	
int Gillespie::r1choice(double a0, double rand,vector<double> avector){
	double thresh=rand*a0;
	int indicator(0);
	for(int i=0;i<(int)avector.size();++i){
		if(partial_sum(avector,i-1)<(thresh)&&(thresh)<=partial_sum(avector,i)){
			indicator = i;
			break;
		}
    }
	return indicator;
}
void Gillespie::initializeData (string consts, vector < double >&reactConsts,vector < int >&specNum)
{
  ifstream myfile;
  myfile.open (consts + ".txt");
  if(!myfile.good()){
	  cout<<"Failure to load "+consts<<endl;
  }
  int N;
  myfile >> N;
  int inthold;
  for (int i = 0; i < N; i++)
    {
      myfile >> inthold;
      specNum.push_back (inthold);
    }
  int M;
  myfile >> M;
  double doublehold;
  for (int i = 0; i < M; i++)
    {
      myfile >> doublehold;
      reactConsts.push_back (doublehold);
    }
  myfile.close ();
}

double Gillespie::partial_sum(vector<double> v, int index){
    double sum = 0;
    for (int i = 0; i <= index; ++i){
        sum += v.at(i);
	}
    return sum;
}

double Gillespie::sumupslow(vector<int> specNum, vector<double> reactConsts,vector<vector<int> > propCoeffs,vector<double>& propVector){
	propVector.resize(propCoeffs.size());
	double totalsumhold(0);
		for(int i=0;i<(int) propCoeffs.size();i++){
			double totalCombos(1);
			for(int j=0;j<(int) propCoeffs[i].size();j++){
				if(propCoeffs[i][j]==1){
				totalCombos*=specNum[j];
				}
				else{
					if(propCoeffs[i][j]==2){
						totalCombos*=(double)specNum[j]*(specNum[j])/2.;
					}
				}
			}
			propVector[i]=totalCombos*reactConsts[i];
			totalsumhold+=propVector[i];
		}
	return totalsumhold;
}
	
	
void Gillespie::specChange(vector<int>& specNum, int reactChoice, vector<vector<int> > changeCoeffs){
	for(int i=0;i<(int)changeCoeffs[reactChoice].size();i++){
		specNum[i]-=changeCoeffs[reactChoice][i];
	}
}

void
initializeData (string consts, vector < double >&reactConsts,
		vector < int >&specNum)
{
  ifstream myfile;
  myfile.open (consts + ".txt");
  if(!myfile.good()){
	  cout<<"Failure to load "+consts<<endl;
  }
  int N;
  myfile >> N;
  int inthold;
  for (int i = 0; i < N; i++)
    {
      myfile >> inthold;
      specNum.push_back (inthold);
    }
  int M;
  myfile >> M;
  double doublehold;
  for (int i = 0; i < M; i++)
    {
      myfile >> doublehold;
      reactConsts.push_back (doublehold);
    }
  myfile.close ();
}

void
Gillespie::initializeReactions(string coeffs, vector < vector < int > >&changeCoeffs, vector < vector < int > >&propCoeffs){
  ifstream myfile;
  myfile.open (coeffs + ".txt");
  int inthold;
  myfile >> inthold;
  int N, M;
  N = inthold;
  myfile >> inthold;
  M = inthold;
  int inthold2;
  propCoeffs.resize (M);
  changeCoeffs.resize (M);
  for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			propCoeffs[i].push_back(0);
			changeCoeffs[i].push_back(0);
		}
    }
	for(int i=0;i<M;i++){
		myfile>>inthold;
		int wait=inthold;
		for (int j=0;j<wait;j++){
			myfile>>inthold;
			myfile>>inthold2;
			propCoeffs[i][inthold]=inthold2;
		}
	}
	for(int i=0;i<M;i++){
		myfile>>inthold;
		int wait=inthold;
		for(int j=0;j<wait;j++){
			myfile>>inthold;
			myfile>>inthold2;
			changeCoeffs[i][inthold]=inthold2;
		}
    }
}