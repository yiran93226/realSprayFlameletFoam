//***************** data preparation ******************//
const List<scalar> x = {1, 5, 10, 15, 20, 25, 30};
const List<scalar> r = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
scalar injectorD = 10.5e-3;

//***************************************//
std::vector<std::vector<scalar>> parcelInCell; // store parcel diameter 
std::vector<std::vector<scalar>> particleNumInCell; // store particle number in the parcel
std::vector<std::vector<scalar>> UxVol; // store particle velocity flux (Ux*parcel_volume)
std::vector<std::vector<scalar>> Vol; // store particle volume in a parcel

std::vector<scalar> cellVol; // store the cell volume of the specific position
// std::vector<scalar> sUx(x.size() * r.size(), 0); // particle velocity flux (Ux*parcel_volume) accumulation in a time region
// std::vector<scalar> sVol(x.size() * r.size(), 0); // parcel volume accumulation in a time region

parcelInCell.resize(x.size() * r.size());
particleNumInCell.resize(x.size() * r.size());
UxVol.resize(x.size() * r.size());
Vol.resize(x.size() * r.size());

//*************************THE FIRST LOG FILE************************************//
OFstream d32File
(
    //runTime.path()/runTime.timeName()/("logSummary." + runTime.timeName() + ".csv")
    runTime.path()/("d32." + runTime.timeName() + ".csv")
);

//***************************************//

d32File
    << "TimeStep" << ",";
for (int i = 0; i < x.size(); i++)
{
    for (int j = 0; j < r.size(); j++)
    {
        d32File
            << "D32(" << x[i] << "; " << r[j] << ")" << ",";
    }
}
d32File
    << endl;

//*************************THE SECOND LOG FILE************************************//
OFstream VelocityFile
(
    runTime.path()/("velocity." + runTime.timeName() + ".csv")
);

VelocityFile
    << "TimeStep" << ",";
for (int i = 0; i < x.size(); i++)
{
    for (int j = 0; j < r.size(); j++)
    {
        VelocityFile
            << "U(" << x[i] << "; " << r[j] << ")" << ",";
    }
}
VelocityFile
    << endl;

//*************************THE THIRD LOG FILE************************************//
OFstream VolFluxFile
(
    runTime.path()/("volFlux." + runTime.timeName() + ".csv")
);

VolFluxFile
    << "TimeStep" << ",";
for (int i = 0; i < x.size(); i++)
{
    for (int j = 0; j < r.size(); j++)
    {
        VolFluxFile
            << "vL(" << x[i] << "; " << r[j] << ")" << ",";
    }
}
VolFluxFile
    << endl;

//*************************THE FOURTH LOG FILE************************************//
OFstream VelocityFluctFile
(
    runTime.path()/("velocityFluct." + runTime.timeName() + ".csv")
);

VelocityFluctFile
    << "TimeStep" << ",";
for (int i = 0; i < x.size(); i++)
{
    for (int j = 0; j < r.size(); j++)
    {
        VelocityFluctFile
            << "U'(" << x[i] << "; " << r[j] << ")" << ",";
    }
}
VelocityFluctFile
    << endl;

//*************************THE FIFTH LOG FILE************************************//
// output the diameter and particle number at the first postion (x0,r0)
OFstream particleInfoFile
(
    runTime.path()/("particleInfo_x0r0." + runTime.timeName() + ".csv")
);

particleInfoFile
    << "TimeStep" << "," << "d" << "," << "n" << endl;