#ifndef INVDISTWEIGHTED_H_
#define INVDISTWEIGHTED_H_

#include <vector>
#include <string>
#include <cmath>
#include "Table.h"

class InvDistWeighted {
public:
    explicit InvDistWeighted(std::string zetaIndex);
    InvDistWeighted(const InvDistWeighted&) = delete;
    InvDistWeighted& operator=(const InvDistWeighted&) = delete;
    ~InvDistWeighted() = default;

    void find(double, double, double, double);

    // access after find
    double lookupOmegaYc() const {
        return omegaYc_;
    }

    double lookupT() const {
        return T_;
    }

    double lookupY(int k) const {
        return Y_[k];
    }

    int ntables(int stateFlag) const { return tables_[stateFlag].size(); }

private:
    void collect(const std::string zetaIndex);

    double normalize(const double x, const double min,
                     const double max) const {
        return (x - min) / std::max(max - min, numlimSmall);
    }

    double distance(const double x1, const double y1,
                    const double x2, const double y2) const {
        return std::sqrt( std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2) );
    }

    const double numlimSmall = 1e-12;
    const double numlimGreat = 1e12;

    std::vector<std::vector<Table> > tables_; // tables_[0] denotes gaseous flamelets, and tables_[1] denotes liquid flamelets

    double pw_;

    double omegaYc_;
    double T_;
    std::vector<double> Y_;
};

#endif  // INVDISTWEIGHTED_H_
