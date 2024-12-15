{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            interface.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    Pair<tmp<volScalarField>> vDotP = mixture->vDotP();
    const volScalarField& vDotcP = vDotP[0]();
    const volScalarField& vDotvP = vDotP[1]();
    
    //-ML: reciprocal of beta1 and rho1
    dimensionedScalar invRho1 = 1.0 / rho1;
    volScalarField invBeta = (1.0 / beta1);
    
    //-ML: multiplication of invRho1 and invBeta
    volScalarField invRho1Beta = invBeta * invRho1;

    volScalarField vDotcvP = vDotvP - vDotcP;
    volScalarField EuvDotcvP = pos0(beta1 - 1.0)*vDotcvP;

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
          - (EuvDotcvP)*(mixture->pSat() - rho*gh)  //-ML: p_rgh + rho*gh
          + fvm::Sp(EuvDotcvP, p_rgh)
          + invRho1Beta*parcels.SU_mDotBubble() //-ML: Rate of change of bubble mass
          + invBeta*parcels.Dbetal() //-ML: Total derivative of betal
        );

        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
