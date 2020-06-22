#ifndef TABLE_H
#define TABLE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

class table
{
public:
    // Constructor
    table(std::ifstream& tableFile);

    // Getters
    std::vector<double> getZ() const
    {
        return Z_;
    }

    std::vector<double> getYc() const
    {
        return Yc_;
    }

    std::vector<double> getH() const
    {
        return h_;
    }

    std::vector<double> getT() const
    {
        return T_;
    }

    std::vector<double> getOmegaYc() const
    {
        return omegaYc_;
    }

    std::vector<std::vector<double> > getY() const
    {
        return Y_;
    }
    double getNsp() const
    {
        return Y_.size();
    }

    std::string getFirstLine() const
    {
        return firstLine_;
    }

    // Member functions
    double lookupOmegaYc() const
    {
        return omegaYc_[positionH_] * weightH_ + omegaYc_[positionL_] * weightL_;
    }

    double lookupT() const
    {
        return T_[positionH_] * weightH_ + T_[positionL_] * weightL_;
    }

    double lookupY(size_t i) const
    {
        return Y_[i][positionH_] * weightH_ + Y_[i][positionL_] * weightL_;
    }

    double lookupYc_() const
    {
        return interYc_;
    }

    void find(double Z);


private:
    // Private member functions
    void read();

    // Private member data
    std::vector<double> Z_;
    std::vector<double> Yc_;
    std::vector<double> h_;
    std::vector<double> omegaYc_;
    std::vector<double> T_;
    std::vector<std::vector<double> > Y_;

    double maxZ_;
    double minZ_;
    double interZ_;
    double interYc_;
    double weightL_;
    double weightH_;
    size_t positionL_;
    size_t positionH_;

    size_t lenZ_;
    size_t nColumn_;
    std::string firstLine_;
    std::ifstream& tableFile_;
};

#endif