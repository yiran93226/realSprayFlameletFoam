tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
       mesh,
       fields,
       phi,
       mesh.divScheme("div(phi,Z_Yc)")
    )
);

{
    fvScalarMatrix ZEqn
    (
        fvm::ddt(rho, Z) + mvConvection->fvmDiv(phi, Z)
      - fvm::laplacian(turbulence->muEff(), Z)
     ==
        fvOptions(rho, Z)
      + parcels.Srho()
    );

    ZEqn.relax();

    fvOptions.constrain(ZEqn);
    
    ZEqn.solve(mesh.solver("Z"));
    
    fvOptions.correct(Z);
    
    Z.max(0.0);


    if (turbulenceModelType == "RAS")
    {
        fvScalarMatrix varZEqn
        (
            fvm::ddt(rho, varZ) + mvConvection->fvmDiv(phi, varZ)
          - fvm::laplacian(turbulence->muEff(), varZ)
          - 2.0*turbulence->muEff()*magSqr(fvc::grad(Z))
          + 2.0*rho*turbulence->epsilon()/turbulence->k()*varZ
          - 0.5*varZ.internalField()*parcels.Srho()/(Z.internalField()+SMALL)
        );

        varZEqn.relax();
        
        varZEqn.solve(mesh.solver("varZ"));
        
        varZ.max(0.0);
    }

    else if (turbulenceModelType == "LES")
    {
        //varZ = 0.15 * sqr(turbulence->delta()) * magSqr(fvc::grad(Z));
	varZ = 0;
    }
    

    else if (turbulenceModelType == "laminar")
    {
        // Info<< "DZ min/max : "
        //     << min(turbulence->muEff()/rho).value() << ", "
        //     << max(turbulence->muEff()/rho).value() << endl;
    }
    else 
    {
        FatalErrorInFunction
            << "Turbulence model " << turbulenceModelType
            << " is not implemented in the current flamelet-based solver"
            << abort(FatalError);
    }
}
