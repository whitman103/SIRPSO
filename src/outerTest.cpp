#include <iostream>
#include <vector>

using namespace std;

int main() {

    vector<double> vec = {3, 1, 5};
    vector<vector<double> > outer(vec.size(), vector<double>(vec.size(), 0));

    for (int i = 0; i < (int)vec.size(); i++) {
        for (int j = 0; j < (int)vec.size(); j++) {
            cout << vec[i]*vec[j] << " ";
        }
        cout << endl;
    }

    return 0;
}