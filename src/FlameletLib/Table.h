#ifndef TABLE_H_
#define TABLE_H_

#include <vector>
#include <string>

struct coeffs {
    int lb;
    int ub;
    double lw;
    double uw;
    double C;
    double H;
    double O;
    double T;
    std::vector<double> Y;

    double idpw = 0.0;  // inverse distance pow

    // constructor
    explicit coeffs(int nsp) { Y.resize(nsp); }

    double linterp(const std::vector<double>& vec) const {
        return vec[lb] * lw + vec[ub] * uw;
    }
};


class Table {
public:
    explicit Table(const std::string);
    Table(const Table&) = default;
    Table& operator=(const Table&) = default;
    ~Table() = default;

    void find(int, double);

    int ncoeffs() const { return data.size(); }

    int nsp() const { return nsp_; }

    // discrete data on C-H plane
    std::vector<coeffs> data;

private:
    // private member functions
    void read(const std::string);


    // private member data
    std::vector<double> z_;  // mixture fraction [-] OFFSET-0
    std::vector<double> c_;  // progress variable [-] OFFSET-1
    std::vector<double> h_;  // enthalpy [J/kg] OFFSET-2
    std::vector<double> omegaC_;  // c_source [kg/m3 s] OFFSET-3
    std::vector<double> t_;  // temperature [K] OFFSET-4
    std::vector<std::vector<double>> y_;  // species mass fractions [-] OFFSET-5
    int nsp_;
    double zCutOff_;
};


#endif  // TABLE_H_
