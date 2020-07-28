#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char** argv){
	
	ifstream inFile(argv[1]);
	int inHold;
	inFile>>inHold;
	int size(inHold);
	vector<int> inVec(size,0);
	for(int i=0;i<size;i++){
		inFile>>inVec[i];
		cout<<inVec[i]<<" ";
	}
	inFile.close();
	

	return 0;
}