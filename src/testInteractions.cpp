#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <chrono>
#include <filesystem>

using namespace std;

#include <boost/random/mersenne_twister.hpp>
boost::mt19937 generator;

#include "GillespieFunctions.h"

int main() {

    generator.seed(time(NULL));
    string localPath = std::filesystem::current_path();
    string outputFolder = localPath + "\\interactionFolder";
    std::filesystem::create_directory(outputFolder);
    outputFolder += "\\";
    vector<int> intSpeciesReset = {1, 0, 0, 0, 1, 0, 0, 0, 0};
    vector<int> specNum = intSpeciesReset;
    vector<double> reportTimes = {0, 10, 20, 40, 80};

    Gillespie ReactionObject(outputFolder + "bindingTwoStateCoeffs");
    ReactionObject.initializeData(outputFolder + "bindingTwoStateConsts", ReactionObject.reactConsts, intSpeciesReset);

    if (intSpeciesReset.size() != ReactionObject.changeCoeffs[0].size()) {
        cout << "Reset Species is wrong" << endl;
        return 0;
    }

    vector<double> resetConsts = ReactionObject.reactConsts;

    const int numOfRuns(2500);

    vector<vector<vector<int> > > outData(reportTimes.size() - 1, vector<vector<int> > (specNum.size(), vector<int>(numOfRuns, 0)));

    for (int run = 0; run < numOfRuns; run++) {

        specNum = intSpeciesReset;
        ReactionObject.reactConsts = resetConsts;
        double currentTime(0);
        int reportIndex(0);
        do {
            ReactionObject.reactConsts[6];//=resetConsts[6]*(10.)/(10.+specNum[3]);
            ReactionObject.reactConsts[8];//=resetConsts[8]*(10.)/(10.+specNum[3]);
            tuple<int, double> hold = ReactionObject.PerformTimeStep2(specNum);
            currentTime += get<1>(hold);


            if (get<1>(hold) < 0) {
                for (int index = reportIndex; index < (int)reportTimes.size(); index++) {
                    for (int i = 0; i < (int)specNum.size(); i++) {
                        outData[index][i][run] = specNum[i];
                    }
                }
                reportIndex = (int)reportTimes.size();
            } else {
                ReactionObject.specChange(specNum, get<0>(hold), ReactionObject.changeCoeffs);
                if (currentTime > reportTimes[reportIndex + 1]) {
                    for (int i = 0; i < (int)specNum.size(); i++) {
                        outData[reportIndex][i][run] = specNum[i];
                    }
                    reportIndex++;
                }
            }
        } while (reportIndex < (int)reportTimes.size() - 1);

    }

    string idString(outputFolder + "bindingTwoStateGenes_base");
    ofstream outTest(idString + ".txt");
    for (int i = 1; i < (int)reportTimes.size(); i++) {
        for (int j = 0; j < 2; j++) {
            outTest << reportTimes[i] << " ";
        }
    }
    outTest << endl;
    for (int i = 0; i < numOfRuns; i++) {
        for (int j = 0; j < (int)reportTimes.size() - 1; j++) {
            outTest << outData[j][3][i] << " " << outData[j][8][i] << " ";
        }
        outTest << endl;
    }
    outTest.close();

    return 0;
}