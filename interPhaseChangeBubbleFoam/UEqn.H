    
    beta2 = 1.0 - beta1;
    word betaScheme("div(phi,beta)"); //- ML : modify the flux of beta
    // Need to store rho2Beta2 for ddt(rho, U)
    volScalarField rho2Beta2
    (
        IOobject
        (
            "rhoDisperse",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        beta2*rho2
    );

    tmp<surfaceScalarField> beta2Phi //- ML : Find the contribution of beta on rhoPhi
    (
        fvc::flux
        (
            phi,
            neg(beta1 - 1.0)*beta2, //- ML : Work on lagrangian Cells only
            betaScheme
        )
    );

    surfaceScalarField rho2PhiBeta2
    (
        IOobject
        (
            "rho2PhiBeta2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        dimensionedScalar(dimMass/dimTime, Zero)
    );

    rho2PhiBeta2 = rho2*beta2Phi;
    
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
      - fvm::ddt(rho2Beta2, U) //-ML: subtract effects of 1-beta1
      - fvm::div(rho2PhiBeta2, U) //-ML: subtract effects of 1-beta1
      + turbulence->divDevRhoReff(rho, U)
      + parcels.SU(U) //-ML: Return the momentum source from the lagrangian framework
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    interface.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }
