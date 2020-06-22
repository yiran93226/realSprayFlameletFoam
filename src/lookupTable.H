{
    volScalarField& he = thermo.he();

    // calculated from transport equations
    const scalarField& ZCells = Z.internalField();
    const scalarField& varZCells = varZ.internalField();
    const scalarField& YcCells = Yc.internalField();
    const scalarField& HaCells = Ha.internalField();
    const scalarField& pCells = p.internalField();
    const scalarField& sourceEvapCell = sourceEvap.internalField();

    // lookup-ed values
    scalarField& heCells = he.primitiveFieldRef();
    scalarField& TCells = T.primitiveFieldRef();
    scalarField& omegaYcCells = omegaYc.primitiveFieldRef();

    scalar Zeta = 0.0;
    // double y = 0.0;

    volScalarField Yt(0.0*Y[0]);
    forAll(ZCells, celli)
    {
        Zeta = Foam::sqrt(varZCells[celli]/Foam::max(ZCells[celli]*(1 - ZCells[celli]), SMALL));
        flameletTable.find(sourceEvapCell[celli], ZCells[celli], Zeta, YcCells[celli], HaCells[celli]);
        heCells[celli] = 0.0;
        TCells[celli] = flameletTable.lookupT();
        omegaYcCells[celli] = flameletTable.lookupOmegaYc();
    
        forAll(Y, i)
        {
            scalarField& YCells = Y[i].primitiveFieldRef();
            YCells[celli] = flameletTable.lookupY(i);
            Yt[celli] += Y[i][celli];
        }
        
        forAll(Y, i)
        {
            scalarField& YCells = Y[i].primitiveFieldRef();
            YCells[celli] /= Yt[celli];
            heCells[celli] += YCells[celli]*thermo.composition().Hs(i, pCells[celli], TCells[celli]);
        }

        //constrain Yc value by C*sum(Y_CO2 + Y_H2O)
        // scalar sumY = Y[CO2Index][celli] + Y[H2OIndex][celli];
        // if ( (YcCells[celli] - sumY)/sumY > 0.5 )
        // {
        //     YcCells[celli] = 1.5*sumY;
        // }
        
    }

    forAll(Z.boundaryField(), patchI)
    {        
        const fvPatchScalarField& pZ = Z.boundaryField()[patchI];
        const fvPatchScalarField& pvarZ = varZ.boundaryField()[patchI];
        const fvPatchScalarField& pYc = Yc.boundaryField()[patchI];
        const fvPatchScalarField& pHa = Ha.boundaryField()[patchI];
        const fvPatchScalarField& pp = p.boundaryField()[patchI];
        const fvPatchScalarField& pSourceEvap = sourceEvap.boundaryField()[patchI];

        fvPatchScalarField& phe = he.boundaryFieldRef()[patchI];
        fvPatchScalarField& pT = T.boundaryFieldRef()[patchI];
        fvPatchScalarField& pYt = Yt.boundaryFieldRef()[patchI];

        forAll(pZ, facei)
        {
            Zeta = Foam::sqrt(pvarZ[facei]/Foam::max(pZ[facei]*(1 - pZ[facei]), SMALL));
            flameletTable.find(pSourceEvap[facei], pZ[facei], Zeta, pYc[facei], pHa[facei]);
            phe[facei] = 0.0;
            pT[facei] = flameletTable.lookupT();
            
            forAll(Y, i)
            {
                fvPatchScalarField& pY = Y[i].boundaryFieldRef()[patchI];
                pY[facei] = flameletTable.lookupY(i);
                pYt[facei] += pY[facei];
            }
            
            forAll(Y, i)
            {
                fvPatchScalarField& pY = Y[i].boundaryFieldRef()[patchI];
                pY[facei] /= pYt[facei];
                phe[facei] += pY[facei]*thermo.composition().Hs(i, pp[facei], pT[facei]);
            }
        }
        //constrain Yc value by C*sum(Y_CO2 + Y_H2O)
        // scalar sumY = Y[CO2Index][facei] + Y[H2OIndex][facei];
        // if ( (pYc[facei] - sumY)/sumY > 0.5 )
        // {
        //     pYc[facei] = 1.5*sumY;
        // }
    }
    
    thermo.correct();
    // radiation->correct();

    Info<< "T gas min/max   " << min(T).value() << ", "
        << max(T).value() << endl;
    Info<< "p min/max : " << min(p).value() << ", "
        << max(p).value() << endl;
    Info<< "Y oxidizer min/max : " << min(Y[oxidizerIndex]).value() << ", "
        << max(Y[oxidizerIndex]).value() << endl;
    Info<< "Sum Y min/max : " << min(Yt).value() << ", "
        << max(Yt).value() << endl;
    Info<< "OmegaYc min/max : " << min(omegaYc).value() << ", "
        << max(omegaYc).value() << endl;

    flameletTable.find(0, 0.1, 0, 0.261, -516480.0);
    Info<< "Find Z0.1 TEMPERATURE: " << flameletTable.lookupT() << endl;
    // Info<< "Find Z0.1 omegaYc: " << flameletTable.lookupOmegaYc() << endl;
}
