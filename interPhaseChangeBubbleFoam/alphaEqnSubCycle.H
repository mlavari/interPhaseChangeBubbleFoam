{
    // Standard face-flux compression coefficient
    surfaceScalarField phic(interface.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (interface.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    volScalarField divU
    (
        mesh.moving()
      ? fvc::div(phi + mesh.phi())
      : fvc::div(phi)
    );

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum
        (
            mesh.newIOobject("rhoPhiSum"),
            mesh,
            dimensionedScalar(rhoPhi.dimensions(), Zero)
        );

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
    }
    else
    {
        #include "alphaEqn.H"
    }

    //- ML: update the gamma value from multiplication of alpha1_*beta1_
    mixture->update_gamma();
    
    //- ML: updating the gamma
    rho == gamma1*rho1 + gamma2*rho2;
}
