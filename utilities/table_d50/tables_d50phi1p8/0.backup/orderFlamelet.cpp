#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{
    size_t nTable = 0;
    vector<int> gasNum, liqNum;

    while (true)
    {
        ifstream inTableFile("noOrdered/flameletTable_" + std::to_string(nTable) + ".csv");
        if (!inTableFile)
            break;
        
        string line, str;
        getline(inTableFile, line); // The first line
        getline(inTableFile, line); // The second line
        istringstream buffer(line);
        getline(buffer, str, ',');
        // cout << "Z[0] = " << stod(str) << endl;
        if (stod(str) > 0.9){
            gasNum.emplace_back(nTable);
        } 
        else{
            liqNum.emplace_back(nTable);
        }

        nTable++;
    }
    // cout << "total number of tables = " << nTable << endl;
    // cout << "gas flamelet number = ";
    // for (int i = 0; i < gasNum.size(); i++){
    //     cout << gasNum[i] << ",";
    // }
    // cout << endl;

    // ======== rename/reorder the flamelets, liquid first ===========
    // 1. liquid flamelets
    int n_liquid = 0;
    for (auto& i : liqNum){        
        cout << i << ",";
        string oldName = "noOrdered/flameletTable_" + std::to_string(i) + ".csv";
        string newName = "ordered/" + std::to_string(n_liquid) + ".csv";
        if (!rename(oldName.c_str(), newName.c_str()))
        {
            std::cout << "rename success "<< std::endl;
        }
        n_liquid++;
    }
    // 2. gas flamelets
    int n_gas = n_liquid;
    for (auto& i : gasNum){        
        cout << i << ",";
        string oldName = "noOrdered/flameletTable_" + std::to_string(i) + ".csv";
        string newName = "ordered/" + std::to_string(n_gas) + ".csv";
        if (!rename(oldName.c_str(), newName.c_str()))
        {
            std::cout << "rename success "<< std::endl;
        }
        n_gas++;
    }
    
    // 3. rename again
    for (int i = 0; i < nTable; i++){
        string oldName = "ordered/" + std::to_string(i) + ".csv";
        string newName = "ordered/flameletTable_" + std::to_string(i) + ".csv";
        if (!rename(oldName.c_str(), newName.c_str()))
        {
            cout << i << "," << " 2nd rename success "<< std::endl;
        }
    }
    

    return(0);
}