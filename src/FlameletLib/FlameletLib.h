#ifndef FLAMELETLIB_H_
#define FLAMELETLIB_H_

#include "InvDistWeighted.h"

class FlameletLib
{
public:
    // Constructor
    explicit FlameletLib(const std::string type);
    FlameletLib(const FlameletLib&) = delete;
    FlameletLib& operator=(const FlameletLib&) = delete;
    ~FlameletLib() = default;

    // Member functions
    void find(double sourceEvap, double Z, double Zeta, double Yc, double h);

    // size_t nsp() const
    // {
    //     return tableSolvers_[0]->nsp();
    // }

    double lookupT() const
    {
        return InvDistWeighteds_[positionH_]->lookupT() * weightH_ +
                InvDistWeighteds_[positionL_]->lookupT() * weightL_;
    }

    double lookupOmegaYc() const
    {
        return InvDistWeighteds_[positionH_]->lookupOmegaYc() * weightH_ +
                InvDistWeighteds_[positionL_]->lookupOmegaYc() * weightL_;
    }

    double lookupY(size_t i) const
    {
        return InvDistWeighteds_[positionH_]->lookupY(i) * weightH_ + 
                InvDistWeighteds_[positionL_]->lookupY(i) * weightL_;
    }

private:
 
    // Private member data
    double weightL_;
    double weightH_;
    size_t positionL_;
    size_t positionH_;
    double Z_;
    double Yc_;
    double minZeta_;
    double maxZeta_;
    double interZeta_;
    size_t num_;
    std::vector<InvDistWeighted*> InvDistWeighteds_;
    std::vector<double> Zetas_;

    //Private member function
    std::vector<std::string> list_dir(const char *path);
};

#endif