{
    fvScalarMatrix HaEqn
    (
        fvm::ddt(rho, Ha) + mvConvection->fvmDiv(phi, Ha)
      - fvm::laplacian(turbulence->alphaEff(), Ha)
     ==
        parcels.hsTrans() / (mesh.V() * runTime.deltaT())
      + parcels.Srho() * dimensionedScalar("h0", dimEnergy/dimMass, -5.1098e06)
      + fvOptions(rho, Ha)
    );

    HaEqn.relax();

    fvOptions.constrain(HaEqn);

    HaEqn.solve();

    fvOptions.correct(Ha);


    Info<< "==== Ha min/max   " << min(Ha).value() << ", "
        << max(Ha).value() << endl;
}

