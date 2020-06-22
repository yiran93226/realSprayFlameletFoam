#include "FlameletLib.h"
#include <iostream>
#include <dirent.h>
#include <sys/types.h>
#include <vector>
#include <algorithm>
#include <string>

FlameletLib::FlameletLib(std::string type)
{
    std::vector<std::string> foldNames = list_dir("tables/");
    std::sort(foldNames.begin(), foldNames.end());   
    foldNames.erase(foldNames.begin(), foldNames.begin() + 2); // remove folder "." and ".."

    for (size_t i=0; i < foldNames.size(); i++)
    {
        std::cout << foldNames[i] << " ";
    }
    std::cout << std::endl;
    
    num_ = foldNames.size();
    
    InvDistWeighteds_.resize(num_);
    Zetas_.resize(num_);
    for (size_t i=0; i<num_; i++)
    {        
        
        Zetas_[i] = stod(foldNames[i]);  
        //nvDistWeighteds_.emplace_back(i);
        InvDistWeighted* tTableS = new InvDistWeighted(foldNames[i]);
        
        std::cout << "#FLAMELETTABLE Constructed" << std::endl;
        InvDistWeighteds_[i] = tTableS;
    }
    minZeta_ = Zetas_[0];
    maxZeta_ = Zetas_[num_-1];
}

void FlameletLib::find(double sourceEvap, double Z, double Zeta, double Yc, double h)
{
    if ( (Zeta > minZeta_) && (Zeta < maxZeta_) )
    {
        interZeta_ = Zeta;
        for (size_t i=0; i<num_; i++)
        {
            if (Zetas_[i] > interZeta_)
            {
                positionL_ = i-1;
                positionH_ = i;
                break;
            }
        }
        weightL_ = (Zetas_[positionH_] - interZeta_) / (Zetas_[positionH_] - Zetas_[positionL_]);
        weightH_ = (interZeta_ - Zetas_[positionL_]) / (Zetas_[positionH_] - Zetas_[positionL_]);
    }
    else if(Zeta <= minZeta_)
    {
        interZeta_ = minZeta_;
        positionL_ = 0;
        positionH_ = 0;
        weightL_ = 0.5;
        weightH_ = 0.5;
    }
    else
    {
        interZeta_= maxZeta_;
        positionL_ = num_-1;
        positionH_ = num_-1;
        weightL_ = 0.5;
        weightH_ = 0.5;
    }
    
    InvDistWeighteds_[positionL_]->find(sourceEvap, Z, Yc, h);
    InvDistWeighteds_[positionH_]->find(sourceEvap, Z, Yc, h);

    // if (Z>1e-6 && Z<1-1e-6)
    // {
    //     InvDistWeighteds_[positionL_].find(Z, Yc, h);
    //     InvDistWeighteds_[positionH_].find(Z, Yc, h);
    // }
    // else
    // {
    //     InvDistWeighteds_[positionL_].bdSet(Z);
    //     InvDistWeighteds_[positionH_].bdSet(Z);
    // }
}

// Private member function
std::vector<std::string> FlameletLib::list_dir(const char *path)
{
    struct dirent *entry;
    DIR *dir = opendir(path);
    std::vector<std::string> foldNames;
    int sum =0;

    if (dir == NULL)
    {
        return foldNames;
        std::cout << "no folder existed" << std::endl;
    }
    while ((entry = readdir(dir)) != NULL)
    {
        // cout << entry->d_name << endl;
        foldNames.emplace_back(entry->d_name);
        sum++;
    }
    closedir(dir);
    std::cout << "sum of the folders = " << sum << std::endl;
    return foldNames;
}