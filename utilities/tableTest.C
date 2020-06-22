#include <iostream>
#include "scalar.H"
#include "List.H"
#include "OFstream.H"
#include "Gamma.h"
#include "table.H"
#include <bits/stdc++.h> 
#include <sys/stat.h> 
#include <sys/types.h> 

using namespace Foam;

void flipTable(List<List<scalar>> &, bool);
void betaPDFIntegration(const scalar &Zeta, const List<List<scalar>> &singleData_, List<List<scalar>> &integratedData_);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
   // const List<scalar> Zetas = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99};
   const List<scalar> Zetas = {0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.25, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99};
   //const List<scalar> Zetas = {0}; // test

   size_t nTable = 0;
   while (true)
   {
      // std::ifstream inTableFile("pdfTables_pilot_Yc/Zeta_0/flameletTable_" + std::to_string(nTable) + ".csv");
      std::ifstream inTableFile("table_d50/phi1.8+phi2.5/Zeta_0/flameletTable_" + std::to_string(nTable) + ".csv");
      if (!inTableFile)
         break;

      table t_0(inTableFile);
      std::vector<double> Z = t_0.getZ();
      std::vector<double> Yc = t_0.getYc();
      std::vector<double> T = t_0.getT();
      std::vector<double> h = t_0.getH();
      std::vector<double> omegaYc = t_0.getOmegaYc();
      std::vector<std::vector<double>> Y = t_0.getY();
      double nsp = t_0.getNsp();

      std::string firstLine = t_0.getFirstLine();

      //  Info<< "species number = " << nsp << endl;

      List<scalar> single(Z.size(), 1.0);
      List<List<scalar>> singleData_(5 + nsp, single);

      // data assignment from the csv file
      for (size_t j = 0; j < Z.size(); j++)
      {
         singleData_[0][j] = Z[j];
      }

      for (size_t j = 0; j < Yc.size(); j++)
      {
         singleData_[1][j] = Yc[j];
      }

      for (size_t j = 0; j < omegaYc.size(); j++)
      {
         singleData_[2][j] = h[j];
      }

      for (size_t j = 0; j < omegaYc.size(); j++)
      {
         singleData_[3][j] = omegaYc[j];
         // std::cout << singleData_[3][j] << std::endl;
      }

      for (size_t j = 0; j < T.size(); j++)
      {
         singleData_[4][j] = T[j];
      }

      for (size_t i = 0; i < nsp; i++)
      {
         for (size_t j = 0; j < Y[i].size(); j++)
         {
            singleData_[i + 5][j] = Y[i][j];
         }
      }

      bool flag = false;
      std::cout << "table number = " << nTable << ",singleData_[0][0] = " << singleData_[0][0] << std::endl;

      if (singleData_[0][0] > 0.9)
      {
         std::cout << "flip number = " << nTable << std::endl;
         flag = true;
      }
      // before pdf, flip the table if flag == false
      std::cout << "flag = " << flag << std::endl;

      flipTable(singleData_, flag);

      List<List<scalar>> integratedData_(singleData_);

      for (int i = 0; i < Zetas.size(); i++)
      {
         const scalar Zeta = Zetas[i];
         //**************** separate the spray table at Z_max****************
         if (flag == false) //for two-phase table
         {
            std::cout << "separate number = " << nTable << std::endl;
            auto ZmaxPtr = max_element(Z.begin(), Z.end());
            int pos1 = ZmaxPtr - Z.begin();
            int pos2 = Z.end() - ZmaxPtr;
            std::cout << "Z.size()= " << Z.size() << std::endl;
            std::cout << "Z.end() - ZmaxPtr = " << Z.end() - ZmaxPtr << std::endl;
            std::cout << "ZmaxPtr - Z.begin() = " << ZmaxPtr - Z.begin() << std::endl;
            std::cout << *ZmaxPtr << " at the position 1 = " << pos1  << ", position 2 = " << pos2 << std::endl;
            
            // List<List<scalar> > table_1, table_2, table_all;
            List<scalar> single_1(pos1 + 1);
            // single.resize(pos1 + 1);
            List<List<scalar>> table_1(integratedData_.size(), single_1);

            List<scalar> single_2(Z.size() - pos1 - 1);
            // single.resize(Z.size() - pos1 - 1);
            List<List<scalar>> table_2(integratedData_.size(), single_2);

            for (int i = 0; i < integratedData_.size(); i++)
            {
               for (int j = 0; j < pos1 + 1; j++)
               {                  
                  table_1[i][j] = integratedData_[i][j];
               }

               // for (int j = 0; j < Z.end() - ZmaxPtr; j++)
               for (int j = 0; j < pos2 - 1; j++)
               {
                  table_2[i][j] = integratedData_[i][Z.size() - j - 1];
               }
            }
            
            //pdf
            List<List<scalar>> tableCopy_1(table_1);
            betaPDFIntegration(Zeta, tableCopy_1, table_1);

            List<List<scalar>> tableCopy_2(table_2);
            betaPDFIntegration(Zeta, tableCopy_2, table_2);

            // output csv file
            string fileName = "table_d50/" + std::to_string(Zeta).substr(0, 4);
            mkdir(fileName.c_str(), 0777);
            string fileName_liquid = "table_d50/" + std::to_string(Zeta).substr(0, 4) + "/liquid";
            mkdir(fileName_liquid.c_str(), 0777);
            string fileName_gas = "table_d50/" + std::to_string(Zeta).substr(0, 4) + "/gas";
            mkdir(fileName_gas.c_str(), 0777);
            // if (mkdir(fileName.c_str(), 0777) == -1)
            // {
            //    // cerr << "Error: " << strerror(errno) << std::endl;
            //    cout << "Dirctory " << fileName << " has existed." << std::endl;
            // }

            // std::ofstream outTableFile1("pdfTables_pilot_Yc/Zeta_" + std::to_string(i + 1) + "/table_1.csv");
            // std::ofstream outTableFile1("tables_d50phi1p8/Zeta_" + std::to_string(i + 1) + "/table_1.csv");
            std::ofstream outTableFile1(fileName_liquid + "/flameletTable_" + std::to_string(nTable)+ ".csv");
            if (!outTableFile1)
            {
               std::cerr << "Unable to save data!\n";
            }
            else
            {
               outTableFile1 << firstLine << std::endl;
               for (int j = 0; j < table_1[0].size(); j++)
               {
                  for (int k = 0; k < table_1.size(); k++)
                  {
                     if (k != table_1.size() - 1)
                     {
                        outTableFile1 << table_1[k][j] << ',';
                     }
                     else
                     {
                        outTableFile1 << table_1[k][j];
                     }
                  }
                  outTableFile1 << std::endl;
               }
            }

            // output csv file
            // std::ofstream outTableFile2("pdfTables_pilot_Yc/Zeta_" + std::to_string(i + 1) + "/table_2.csv");
            // std::ofstream outTableFile2("tables_d50phi1p8/Zeta_" + std::to_string(i + 1) + "/table_2.csv");
            // fileName = "tables_d50phi1p8/" + std::to_string(Zeta).substr(0, 4) + "/gas";
            // if (mkdir(fileName.c_str(), 0777) == -1)
            // {
            //    // cerr << "Error: " << strerror(errno) << std::endl;
            //    cout << "Dirctory " << fileName << " has existed." << std::endl;
            // }
            std::ofstream outTableFile2(fileName_gas +  "/flameletTable_" + std::to_string(nTable)+ ".csv");            
            
            if (!outTableFile2)
            {
               std::cerr << "Unable to save data!\n";
            }
            else
            {
               outTableFile2 << firstLine << std::endl;
               for (int j = 0; j < table_2[0].size(); j++)
               {
                  for (int k = 0; k < table_2.size(); k++)
                  {
                     if (k != table_2.size() - 1)
                     {
                        outTableFile2 << table_2[k][j] << ',';
                     }
                     else
                     {
                        outTableFile2 << table_2[k][j];
                     }
                  }
                  outTableFile2 << std::endl;
               }
            }

            //*********************** merge the table *********************
            // for (int i = 0; i < integratedData_.size(); i++)
            // {
            //    for (int j = 0; j < table_1[0].size(); j++)
            //    {
            //       integratedData_[i][j] = table_1[i][j];
            //    }

            //    // table_2 should be flipped again
            //    for (int j = 0; j < table_2[0].size(); j++)
            //    {
            //       integratedData_[i][j + table_1[0].size()] = table_2[i][table_2[0].size() - j - 1];
            //       // integratedData_[i][j + table_1[0].size()] = table_2[i][j];
            //    }
            // }          

         }

         if (flag == true)  // flag == true for one-phase table
         {
            betaPDFIntegration(Zeta, singleData_, integratedData_);

            // output csv file
            // std::ofstream outTableFile("pdfTables_pilot_Yc/Zeta_" + std::to_string(i + 1) + "/flameletTable_" + std::to_string(nTable) + ".csv");
            string fileName = "table_d50/" + std::to_string(Zeta).substr(0, 4);
            mkdir(fileName.c_str(), 0777);
            string fileName_gas = "table_d50/" + std::to_string(Zeta).substr(0, 4) + "/gas";
            mkdir(fileName_gas.c_str(), 0777);

            std::ofstream outTableFile(fileName_gas + "/flameletTable_" + std::to_string(nTable) + ".csv");
            if (!outTableFile)
            {
               std::cerr << "Unable to save data!\n";
            }
            else
            {
               outTableFile << firstLine << std::endl;
               for (int j = 0; j < integratedData_[0].size(); j++)
               {
                  for (int k = 0; k < integratedData_.size(); k++)
                  {
                     if (k != integratedData_.size() - 1)
                     {
                        outTableFile << integratedData_[k][j] << ',';
                     }
                     else
                     {
                        outTableFile << integratedData_[k][j];
                     }
                  }
                  outTableFile << std::endl;
               }
            }
         }
         
      }

      nTable++;
   }
   Info << "total number of tables = " << nTable << endl;
}

//============================= function definitions =============================%
void flipTable(List<List<scalar>> &integratedData_, bool flag)
{
   // determine if flipping the table
   if (flag)
   {
      List<List<scalar>> copyData_(integratedData_);
      std::cout << "integratedData_.size() = " << integratedData_.size() << std::endl;
      std::cout << "integratedData_[0].size() = " << integratedData_[0].size() << std::endl;

      for (int i = 0; i < integratedData_.size(); i++)
      {
         for (int j = 0; j < integratedData_[i].size(); j++)
            integratedData_[i][j] = copyData_[i][integratedData_[i].size() - 1 - j];
      }
   }
}

void betaPDFIntegration(const scalar &Zeta, const List<List<scalar>> &singleData_, List<List<scalar>> &integratedData_)
{
   if (Zeta != 0)
   {
      List<scalar> Z_(integratedData_[0]);
      List<scalar> varZ_(integratedData_[0].size(), 0.0);
      List<scalar> pdfAlpha(integratedData_[0].size(), 0.0);
      List<scalar> pdfBeta(integratedData_[0].size(), 0.0);
      List<scalar> PDF(integratedData_[0].size(), 0.0);

      // define the maximum value of Z, used in allocating Z
      // Z in table must be increasingly ordered
      scalar maxZ_ = Z_[Z_.size() - 1];
      // scalar maxZ_ = 1;
      std::cout << "maxZ_ = " << maxZ_ << std::endl;

      for (int i = 0; i < integratedData_[0].size(); i++)
      {
         varZ_[i] = sqr(Zeta) * (Z_[i] * (1.0 - Z_[i]));

         if (varZ_[i] > 1e-7)
         {            
            pdfAlpha[i] = Z_[i] * ((Z_[i] * (1.0 - Z_[i])) / varZ_[i] - 1);
            pdfBeta[i] = (1.0 - Z_[i]) * ((Z_[i] * (1.0 - Z_[i])) / varZ_[i] - 1);
            // pdfAlpha[i] = Z_[i] * (1.0/sqr(Zeta) - 1);
            // pdfBeta[i] = (1.0 - Z_[i]) * pdfAlpha[i] / Z_[i];

            // Limit alpha and beta but keep their ratio
            if (pdfAlpha[i] > 500)
            {
               pdfBeta[i] = pdfBeta[i] / pdfAlpha[i] * 500;
               pdfAlpha[i] = 500;
            }
            if (pdfBeta[i] > 500)
            {
               pdfAlpha[i] = pdfAlpha[i] / pdfBeta[i] * 500;
               pdfBeta[i] = 500;
            }

            int gridPoints = 500;
            List<long double> hZ_(gridPoints, 0.0);

            for (int j = 0; j < gridPoints; j++)
            {
               hZ_[j] = 1e-4 + j * (maxZ_ - 1e-4) / (gridPoints - 1);
            }         

            if ((pdfAlpha[i] > 1) && (pdfBeta[i] > 1))
            {               
               // Calculate BetaPDF
               PDF.clear();
               PDF.resize(hZ_.size(), 0.0);

               for (int j = 0; j < hZ_.size(); j++)
               {
                  PDF[j] = (std::pow(hZ_[j], (pdfAlpha[i] - 1e0))) * std::pow((1e0 - hZ_[j]), (pdfBeta[i] - 1e0)) * Gamma(pdfAlpha[i] + pdfBeta[i]) / (min(Gamma(pdfAlpha[i]), 1e17) * min(Gamma(pdfBeta[i]), 1e17));
               }
            }

            else if ((pdfAlpha[i] <= 1) && (pdfBeta[i] > 1))
            {
               // PDF Singularity at Z = 0

               // Calculate BetaPDF
               PDF.clear();
               PDF.resize(hZ_.size(), 0.0);
               for (int j = 1; j < hZ_.size(); j++)
               {
                  PDF[j] = std::pow(hZ_[j], (pdfAlpha[i] - 1e0)) * std::pow((1e0 - hZ_[j]), (pdfBeta[i] - 1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]), 1e17) / (min(Gamma(pdfAlpha[i]), 1e17) * min(Gamma(pdfBeta[i]), 1e17));
               }
               PDF[0] = 1.5 * PDF[1] / pdfAlpha[i];
            }

            else if ((pdfAlpha[i] > 1) && (pdfBeta[i] <= 1))
            {
               // PDF Singularity at Z = 1
               // Allocation of Z for PDF integration
               
               // Calculate BetaPDF
               PDF.clear();
               PDF.resize(hZ_.size(), 0.0);

               for (int j = 0; j < hZ_.size() - 1; j++)
               {
                  PDF[j] = std::pow((hZ_[j]), (pdfAlpha[i] - 1e0)) * std::pow((1e0 - hZ_[j]), (pdfBeta[i] - 1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]), 1e17) / (min(Gamma(pdfAlpha[i]), 1e17) * min(Gamma(pdfBeta[i]), 1e17));
               }
               PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta[i];
            }

            else if ((pdfAlpha[i] <= 1) && (pdfBeta[i] <= 1))
            {
               // PDF Singularity at Z = 1 and Z = 0

               // Calculate BetaPDF
               PDF.clear();
               PDF.resize(hZ_.size(), 0.0);

               for (int j = 1; j < hZ_.size() - 1; j++)
               {
                  PDF[j] = std::pow(hZ_[j], (pdfAlpha[i] - 1e0)) * std::pow((1e0 - hZ_[j]), (pdfBeta[i] - 1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]), 1e17) / (min(Gamma(pdfAlpha[i]), 1e17) * min(Gamma(pdfBeta[i]), 1e17));
                  // std::cout << "pdfBeta = " << pdfBeta[i] << "\n" << std::endl;
               }
               PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta[i];
               PDF[0] = 1.5 * PDF[1] / pdfAlpha[i];
            }
            

            // ============================ Calculate the area of the PDF for scaling =================================
            scalar intPDF = 0;
            for (int j = 1; j < hZ_.size(); j++)
            {
               intPDF += (hZ_[j - 1] - hZ_[j]) * (PDF[j - 1] + PDF[j]) / 2;
            }

            // Interpolate singleData entries to the new mixture fraction space
            List<scalar> hY_(gridPoints, 0.0);
            scalar intY = 0;

            for (int j = 0; j < integratedData_.size(); j++)
            {
               hY_ = 0.0;
               hY_[0] = singleData_[j][0];
               hY_[hY_.size() - 1] = singleData_[j][singleData_[j].size() - 1];
               
               intY = 0;

               for (int k = 1; k < hZ_.size() - 1; k++)
               {
                  int ubZ = 0;
                  for (int l = 0; l < Z_.size(); l++)
                  {
                     ubZ = l;
                     if (hZ_[k] < Z_[l])
                     {
                        break;
                     }                           
                  }
                  int lbZ = ubZ - 1;

                  // Interpolation to hZ space
                  hY_[k] = (singleData_[j][ubZ] - singleData_[j][lbZ]) / max(Z_[ubZ] - Z_[lbZ], SMALL) * (hZ_[k] - Z_[lbZ]) + singleData_[j][lbZ];
                  
                  // PDF Integration using the trapezoidal rule
                  intY += (hZ_[k - 1] - hZ_[k]) * (hY_[k - 1] * PDF[k - 1] + hY_[k] * PDF[k]) / (2.0 * intPDF);
               }
               
               // Special treatment for the boundaries
               intY += (hZ_[hZ_.size() - 2] - hZ_[hZ_.size() - 1]) * (hY_[hZ_.size() - 2] * PDF[hZ_.size() - 2] + hY_[hZ_.size() - 1] * PDF[hZ_.size() - 1]) / (2.0 * intPDF);

               if (singleData_[0][singleData_[0].size() - 1] > 0.5) //gas phase
               {
                  if (i != 0 && i != integratedData_[0].size() - 1 && j > 2) //no transfer of Z, Yc and Ha columns
                  // if (i != 0 && i != integratedData_[0].size() - 1)
                     integratedData_[j][i] = intY;
               }
               else
               {
                  if (i != 0 && j > 2) //no transfer of Z, Yc and Ha columns
                     integratedData_[j][i] = intY;
               }
            }
         }
      }
   }
}
