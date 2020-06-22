#include "InvDistWeighted.h"
#include <iostream>
#include <stdexcept>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::to_string;
using std::runtime_error;
using std::min;
using std::max;
using std::pow;


InvDistWeighted::InvDistWeighted(std::string zetaIndex)
    : pw_(-3)
{
    collect(zetaIndex);
    if (ntables(0) <= 0 || ntables(1) <= 0) throw runtime_error("flameletTable not found");
    else Y_.resize(tables_[0][0].nsp()); // arbitrary table is OK
}


void InvDistWeighted::find(double sourceEvap, double z, double c, double h)
{
    int stateFlag = 0; // gaseous flamelets
    if (sourceEvap > 5) // liquid flamelets
    {
        stateFlag = 1;
    }
    
    // 1-find in each table and calculate minmax
    double cmin = numlimGreat;
    double cmax = -numlimGreat;
    double hmin = numlimGreat;
    double hmax = -numlimGreat;
    for (int i = 0; i < ntables(stateFlag); i++) {
        tables_[stateFlag][i].find(stateFlag, z);
        for (int j = 0; j < tables_[stateFlag][i].ncoeffs(); j++) {
            cmin = min(tables_[stateFlag][i].data[j].C, cmin);
            cmax = max(tables_[stateFlag][i].data[j].C, cmax);
            hmin = min(tables_[stateFlag][i].data[j].H, hmin);
            hmax = max(tables_[stateFlag][i].data[j].H, hmax);
        }
    }

    // cout << "cmin|cmax|hmin|hmax: " << cmin << '|'
    //      << cmax << '|' << hmin << '|' << hmax << endl;

    // 2-normalize and calculate distance
    // dpow[i] = (distance[i])^pw_
    // denominator = \sum{dpow[i]}
    c = max(cmin, min(c, cmax));
    h = max(hmin, min(h, hmax));
    const double cnorm = normalize(c, cmin, cmax);
    const double hnorm = normalize(h, hmin, hmax);
    // cout << "cnorm|hnorm: " << cnorm << '|' << hnorm << endl;
    double denominator = 0.0;
    for (int i = 0; i < ntables(stateFlag); i++) {
        for (int j = 0; j < tables_[stateFlag][i].ncoeffs(); j++) {
            double cj = normalize(tables_[stateFlag][i].data[j].C, cmin, cmax);
            double hj = normalize(tables_[stateFlag][i].data[j].H, hmin, hmax);
            double d = max(distance(cnorm, hnorm, cj, hj), numlimSmall);
            double dpow = pow(d, pw_);
            tables_[stateFlag][i].data[j].idpw = dpow;
            denominator += dpow;
        }
    }
    denominator = max(denominator, numlimSmall);

    // 3-compute lookup values
    omegaYc_ = 0.0; T_ = 0.0;
    int nsp = Y_.size();
    for (int k = 0; k < nsp; k++)
        Y_[k] = 0.0;
    for (int i = 0; i < ntables(stateFlag); i++) {
        for (int j = 0; j < tables_[stateFlag][i].ncoeffs(); j++) {
            double weight = tables_[stateFlag][i].data[j].idpw / denominator;
            omegaYc_ += weight * tables_[stateFlag][i].data[j].O;
            T_ += weight * tables_[stateFlag][i].data[j].T;
            for (int k = 0; k < nsp; k++)
                Y_[k] += weight * tables_[stateFlag][i].data[j].Y[k];
        }
    }
}

void InvDistWeighted::collect(const std::string zetaIndex)
{
    std::vector<string> stateTable = {"/gas", "/liquid"};
    std::vector<Table> table1;
    for (int i = 0; i < 2; i++)
    {
        int n = 0;
        table1.clear();
        while (true)
        {
            string filename = ("tables/" + zetaIndex + stateTable[i] + "/flameletTable_" + to_string(n) + ".csv");
            // std::cout << filename << std::endl;
            //dir + "/flameletTable_" + to_string(n) + ".csv";
            
            try
            {
                table1.emplace_back(filename);
            }
            catch (runtime_error err)
            {
                cout << "zetaIndex = " << zetaIndex << ", Number of tables: " << table1.size() << endl;
                break;
            }
            ++n;
        }
        
        tables_.push_back(table1);
    }
}

