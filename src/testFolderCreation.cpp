#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

int main(int argc, char** argv) {
    string inputString(argv[1]);
    string outputFolder = "testFolder_" + inputString;
    mkdir(outputFolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    outputFolder += "//";
    ofstream outFile(outputFolder + "test.txt");
    outFile << "This is from program :" << argv[1];
    outFile.close();

    return 0;
}