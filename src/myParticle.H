// track all positions of parcels and extract the parcels in a certain cell

// Lsit x and r are defined in startSummary.H

int timeNumber = 10; // 

if (tcount % timeNumber == 0) // to ensure there are no same parcels occured in the cell at the previous time step
{
     int recordNumber = 5; // write to file every n records
     int deltaTimeStep = recordNumber*timeNumber; 

     if (tcount % deltaTimeStep == 0)
     {
          d32File()
              << tcount/timeNumber << ","; // instead of tcount
          VelocityFile()
              << tcount/timeNumber << ",";
          VolFluxFile()
              << tcount/timeNumber << ",";
          VelocityFluctFile()
              << tcount/timeNumber << ",";
     }

     //*********** store all parcels in an accumulative time region (deltaTimeStep) *************//
     for (int i = 0; i < x.size(); i++)
     {
          for (int j = 0; j < r.size(); j++)
          {
               int num = r.size() * i + j;
               //Info << "i = " << i << ", j = " << j << ", num = " << num << endl;

               vector position(x[i] * injectorD, r[j] * injectorD, 0);
               label cell_i = mesh.findCell(position);

               //************* record all parcels in the certain cell ***********************//
               scalar sUxVol_oneTimeStep = 0.0;
               scalar sVol_oneTimeStep = 0.0;

               forAllIter(basicSprayCloud, parcels, pIter)
               {
                    label cellI = pIter().cell();
                    if (cellI == cell_i)
                    {
                         if (pIter().d() > 5 * 1.0e-6) //only size larger than 10um can be stored
                         {
                              parcelInCell[num].push_back(pIter().d());
                              particleNumInCell[num].push_back(pIter().nParticle());

                              sUxVol_oneTimeStep += pIter().U()[0] * pIter().nParticle() * pIter().volume();
                              sVol_oneTimeStep += pIter().nParticle() * pIter().volume();

                              // sUx[num] += pIter().U()[0] * pIter().nParticle() * pIter().volume();
                              // sVol[num] += pIter().nParticle() * pIter().volume();
                         }
                    }
               }
               UxVol[num].push_back(sUxVol_oneTimeStep);
               Vol[num].push_back(sVol_oneTimeStep);
          }
     }

     //************* calculate SMD, Ux, volFlux *********************//
     if (tcount % deltaTimeStep == 0)
     {
          //============ calculate the cell volume in each position ================//
          for (int i = 0; i < x.size(); i++)
          {
               for (int j = 0; j < r.size(); j++)
               {
                    vector position(x[i] * injectorD, r[j] * injectorD, 0);
                    label cell_i = mesh.findCell(position);
                    cellVol.push_back(mesh.V()[cell_i]);
               }
          }

          //************* calculate SMD, Ux, volFlux *********************//
          for (size_t i = 0; i < parcelInCell.size(); i++) //different position
          {
               scalar s3 = 0.0;
               scalar s2 = 0.0;
               scalar sUx = 0.0;
               scalar sVol = 0.0;
               std::vector<scalar> transUx;
               int count = 0;

               for (size_t n = 0; n < parcelInCell[i].size(); n++)
               {
                    // s3 += particleNumInCell[i][n] * pow(parcelInCell[i][n], 3);
                    // s2 += particleNumInCell[i][n] * pow(parcelInCell[i][n], 2);
                    s3 += pow(parcelInCell[i][n], 3);
                    s2 += pow(parcelInCell[i][n], 2);

                    count++;
               }

               for (size_t n = 0; n < UxVol[i].size(); n++)
               {
                    sUx += UxVol[i][n];
                    sVol += Vol[i][n];
                    transUx.push_back(UxVol[i][n] / max(Vol[i][n], vSmall));
                    if (i == 0)
                    {
                         Info << "transient Ux = " << UxVol[i][n] / max(Vol[i][n], vSmall) << nl << endl;
                    }
               }

               if (count)
               {
                    scalar D32 = 1.0e+6 * s3 / max(s2, vSmall);

                    //Info << "timestep = " << tcount << ", position (x,r) = (" << x[i] << "," << r[j] << "), "
                    //<< "D32 = " << D32 << endl;
                    //Info << "timestep = " << tcount << ", position number = " << i << ", D32 = "
                    // << D32 << endl;

                    // scalar Ux = sUx[i] / max(sVol[i], vSmall);
                    scalar Ux = sUx / max(sVol, vSmall);
                    // alpha needs to be divided by deltaTimeStep to represent the average value in each time step
                    // scalar alpha = sVol[i] / max(cellVol[i]*deltaTimeStep, vSmall);
                    scalar alpha = sVol / max(cellVol[i] * recordNumber, vSmall);
                    scalar volFlux = alpha * Ux;

                    scalar sum_UxFluc = 0.0;
                    int countN = 0;
                    for (size_t n = 0; n < transUx.size(); n++)
                    {
                         if (transUx[n] > 0)
                         {
                              sum_UxFluc += std::abs(transUx[n] - Ux);
                              countN++;
                         }
                    }
                    scalar UxFluc = sum_UxFluc / max(countN, 1);

                    d32File()
                        << D32 << ","; //mu
                    VelocityFile()
                        << Ux << ","; //m/s
                    VolFluxFile()
                        << volFlux * 100 << ","; // cm/s
                    VelocityFluctFile()
                        << UxFluc << ","; // m/s
               }
               else
               {
                    Info << "No parcel in the cell " << nl << endl;
                    d32File()
                        << 0 << ",";
                    VelocityFile()
                        << 0 << ",";
                    VolFluxFile()
                        << 0 << ",";
                    VelocityFluctFile()
                        << 0 << ",";
               }
          }
          d32File()
              << endl;
          VelocityFile()
              << endl;
          VolFluxFile()
              << endl;
          VelocityFluctFile()
              << endl;

          //============ write the particle diameter and number at the postion (x0,r0) ================//
          particleInfoFile()
               << tcount/timeNumber << "," << endl;
          for (size_t i = 0; i < parcelInCell[0].size(); i++)
          {
               particleInfoFile()
                    << "," << parcelInCell[0][i] << "," << particleNumInCell[0][i] << endl;
          }
          particleInfoFile() // separate each time recording
               << endl;

          //============== clear the data for the next accumulative time region ===============
          parcelInCell.clear();
          particleNumInCell.clear();
          cellVol.clear();
          UxVol.clear();
          Vol.clear();

          parcelInCell.resize(x.size() * r.size());
          particleNumInCell.resize(x.size() * r.size());
          UxVol.resize(x.size() * r.size());
          Vol.resize(x.size() * r.size());
     }
}
