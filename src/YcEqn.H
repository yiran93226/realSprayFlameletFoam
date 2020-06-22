{
    fvScalarMatrix YcEqn
    (
        fvm::ddt(rho, Yc) + mvConvection->fvmDiv(phi, Yc)
      - fvm::laplacian(turbulence->muEff(), Yc)
     ==
        fvOptions(rho, Yc)
      + omegaYc
    ); 

    YcEqn.relax();

    fvOptions.constrain(YcEqn);

    YcEqn.solve(mesh.solver("Yc"));

    fvOptions.correct(Yc);

    Yc.max(0.0);


    Info<< "Z min/max : " << min(Z).value() << ", "
        << max(Z).value() << endl;

    Info<< "varZ min/max : " << min(varZ).value() << ", "
        << max(varZ).value() << endl;

    Info<< "Yc min/max : " << min(Yc).value() << ", "
        << max(Yc).value() << endl;
}