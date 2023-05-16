#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

void generateOutputFolder(vector<int> inTimes, string outFile);


int main() {


    string outputFolder = "InputFolder";
    mkdir(outputFolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    outputFolder += "//";

    ofstream scriptOutWrite("runScript.sh");
    scriptOutWrite << "#!/bin/sh" << endl;
    scriptOutWrite << "#SBATCH -N 3" << endl;
    scriptOutWrite << "#SBATCH --ntasks=60" << endl;
    scriptOutWrite << "#SBATCH --cpus-per-task=1" << endl;

    int numOfThreads(20);

    string baseString("mpirun -np " + to_string(numOfThreads) + " ./SIRPSO_MPI.exe \"" + to_string(numOfThreads) + "\" \"");

    vector<int> times = {0, 10, 20};
    int runCounter(0);
    generateOutputFolder(times, outputFolder + "outputOne.txt");
    scriptOutWrite << (baseString + outputFolder + "outputOne\" \"" + to_string(runCounter) + "\" &") << endl;
    runCounter++;
    vector<int> times2 = {0, 40, 80};
    generateOutputFolder(times2, outputFolder + "outputTwo.txt");
    scriptOutWrite << (baseString + outputFolder + "outputTwo\" \"" + to_string(runCounter) + "\" &") << endl;
    runCounter++;
    vector<int> times3 = {0, 10, 20, 40, 60, 80};
    generateOutputFolder(times3, outputFolder + "outputThree.txt");
    scriptOutWrite << (baseString + outputFolder + "outputThree\" \"" + to_string(runCounter) + "\" &") << endl;
    scriptOutWrite << "wait" << endl;


    scriptOutWrite.close();


    return 0;
}

void generateOutputFolder(vector<int> inTimes, string outFile) {
    ofstream outStream(outFile);
    outStream << inTimes.size() << endl;
    for (int i = 0; i < (int)inTimes.size(); i++) {
        outStream << inTimes[i] << " ";
    }
    outStream.close();
}