#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include "Table.h"
using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::getline;
using std::count;
using std::stod;
using std::max;
using std::min;
using std::runtime_error;


Table::Table(const string filename)
    : zCutOff_(1e-6)
{
    read(filename);
    data.reserve(2);
}

void Table::find(int stateFlag, double Z)
{
    data.clear();
    if (stateFlag == 0) // gas
    {
        Z = max(0.0, min(Z, 1.0));
    }else
    {
        Z = max(0.0, min(Z, *max_element(z_.begin(),z_.end()) ));
    }
    
    
    
    int len = z_.size();

    if (Z <= zCutOff_) {
        int i = 0;
        while(i < len && z_[i] > zCutOff_)
            ++i;
        if (i < len) {
            data.emplace_back(nsp_);
            data.back().lb = i;
            data.back().ub = i;
            data.back().lw = 0.5;
            data.back().uw = 0.5;
        }
    } else if (Z >= 1.0 - zCutOff_) {
        int i = 0;
        while (i < len && z_[i] < 1.0 - zCutOff_)
            ++i;
        if (i < len) {
            data.emplace_back(nsp_);
            data.back().lb = i;
            data.back().ub = i;
            data.back().lw = 0.5;
            data.back().uw = 0.5;
        }
    } else {
        for (int i = 1; i < len; i++) {
            if (Z > z_[i-1]) {
                if (Z > z_[i]) continue;
                else if (Z < z_[i]) {
                    data.emplace_back(nsp_);
                    data.back().lb = i-1;
                    data.back().ub = i;
                    data.back().lw = (z_[i] - Z)/(z_[i] - z_[i-1]);
                    data.back().uw = (Z - z_[i-1])/(z_[i] - z_[i-1]);
                } else {  // Z == z_[i]
                    data.emplace_back(nsp_);
                    data.back().lb = i;
                    data.back().ub = i;
                    data.back().lw = 0.5;
                    data.back().uw = 0.5;
                }
            } else if (Z < z_[i-1]) {
                if (Z < z_[i]) continue;
                else if (Z > z_[i]) {
                    data.emplace_back(nsp_);
                    data.back().lb = i;
                    data.back().ub = i-1;
                    data.back().lw = (z_[i-1] - Z)/(z_[i-1] - z_[i]);
                    data.back().uw = (Z - z_[i])/(z_[i-1] - z_[i]);
                } else {  // Z == z_[i]
                    data.emplace_back(nsp_);
                    data.back().lb = i;
                    data.back().ub = i;
                    data.back().lw = 0.5;
                    data.back().uw = 0.5;
                }
            } else if (i == 1) {  // Z == z_[i-1] && i == 1
                data.emplace_back(nsp_);
                data.back().lb = i-1;
                data.back().ub = i-1;
                data.back().lw = 0.5;
                data.back().uw = 0.5;
            } else continue;  // Z == z_[i-1] && i != 1
        }
    }

    len = data.size();
    for (int i = 0; i < len; i++) {
        data[i].C = data[i].linterp(c_);
        data[i].H = data[i].linterp(h_);
        data[i].O = data[i].linterp(omegaC_);
        data[i].T = data[i].linterp(t_);
        for (int k = 0; k < nsp_; k++)
            data[i].Y[k] = data[i].linterp(y_[k]);
    }
}


// Input format:
// Z,C,H,Omega,T,Y1,...,Yn
void Table::read(const string filename)
{
    ifstream in;
    in.open(filename);
    if (!in)
        throw runtime_error("Couldn't open: " + filename);
    string line, str;
    getline(in, line);
    nsp_ = (count(line.begin(), line.end(), ',') - 4);
    y_.resize(nsp_);  // nsp
    // cout << "Number of species: " << nsp_ << endl;
    while(getline(in, line)) {
        int n = 0;
        istringstream buffer(line);
        while(getline(buffer, str, ',')) {
            switch (n) {
            case 0:
                z_.push_back(stold(str));
                break;
            case 1:
                c_.push_back(stold(str));
                break;
            case 2:
                h_.push_back(stold(str));
                break;
            case 3:
                omegaC_.push_back(stold(str));
                break;
            case 4:
                t_.push_back(stold(str));
                break;
            default:
                y_[n-5].push_back(stold(str));
                break;
            }
            ++n;
        }
    }
    in.close();
}

