#ifndef GILLESPIE_H
#define GILLESPIE_H

#include <vector>
#include <string>
#include <random>
using namespace std;

class Gillespie{
	
	
	public:
	Gillespie(string);
	~Gillespie();
	tuple<int,double> PerformTimeStep(vector<int>&);
	tuple<int,double> PerformTimeStep2(vector<int>&);
	int muChoice(double,double,vector<double>);
	double partial_sum(vector<double>, int);
	double sumupslow(vector<int>,vector<double>,vector< vector<int> >,vector<double>&);
	void specChange(vector<int>&,int,vector<vector<int> >);
	void initializeReactions(string,vector<vector<int> > &,vector<vector<int> >&);
	int r1choice(double, double,vector<double>);
	mt19937 generator;
	void initializeData (string consts, vector < double >&reactConsts,
	vector<int>&specNum);
	vector<double> reactConsts;
	vector<vector<int> > changeCoeffs;
	
	
	private:
	vector<vector<int> > propCoeffs;
	vector<double> propVector;
	
};


#endif