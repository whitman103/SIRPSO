#include <fstream>
#include <sys/types.h>

int main(int argc, char** argv){
    string outputFolder="testFolder_"+to_string(argv[1]);
	mkdir(outputFolder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    outputFolder+="//";
    ofstream outFile(outputFolder+"test.txt");
    outFile<<"This is from program :"<<argv[1];
    outFile.close();

    return 0;
}