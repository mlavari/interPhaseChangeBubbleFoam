Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Creating hybridPhaseChangeTwoPhaseMixture\n" << endl;
autoPtr<hybridPhaseChangeTwoPhaseMixture> mixture =
    hybridPhaseChangeTwoPhaseMixture::New(U, phi);

volScalarField& alpha1(mixture->alpha1());
volScalarField& alpha2(mixture->alpha2());

//- ML: Adding the beta and gamma
volScalarField& beta1(mixture->beta1());
volScalarField& beta2(mixture->beta2());
volScalarField& gamma1(mixture->gamma1());
volScalarField& gamma2(mixture->gamma2());

//- ML: Double check values
alpha2 = 1.0 - alpha1;
beta2 = 1.0 - beta1;

const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    gamma1*rho1 + gamma2*rho2
);
rho.oldTime();

// Construct interface from alpha1 distribution
interfaceProperties interface(alpha1, U, mixture());

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture())
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

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

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

//- ML: Add setFlux for beta and gamma here
mesh.setFluxRequired(beta1.name());
mesh.setFluxRequired(gamma1.name());

#include "createFvOptions.H"

IOobject alphaPhi10Header
(
    IOobject::groupName("alphaPhi0", alpha1.group()),
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

// MULES flux from previous time-step
surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);


//-ML: Creating elements important for lagrangian bubbles
volScalarField muc
(
    IOobject
    (
        "muc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->nu()*rho
);

//- ML: Creating clouds
const word kinematicBubbleCloudName
(
    args.getOrDefault<word>("cloud", "kinematicBubbleCloud")
);

Info<< "Constructing Cloud: " << kinematicBubbleCloudName << endl;

basicKinematicBubbleCloud parcels
(
    kinematicBubbleCloudName,
    rho,
    U,
    muc,
    p,//-ML:  It is added for bubble Dynamics-Rayleigh Plesset Equation
    alpha1,//-ML: It is added to modify liquid volume fraction
    beta1,//-ML: It is added to modify liquid volume fraction
    g
);


//-ML: Set the Euler to Lagrang controllers
Info<< "Set Eulerian-to-Lagrangian Controllers " <<  endl;
IOdictionary EulerToLagrang
(
    IOobject
    (
        "kinematicBubbleCloudProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dictionary readEulerToLagrang(EulerToLagrang.subDict("Transition_Algorithm"));

//- ML: Flag for activation of EulerianToLagrangian
label EulerToLagrang_activation(readEulerToLagrang.lookupOrDefault<bool>("EulerianToLagrangian_activation", 0));
Info << "Set the EulerToLagrang_activation: "  <<  EulerToLagrang_activation  << endl;

//- ML: Number of cells that trigger the transtion from Eulerian to Lagrangian
scalar EL_cellThreshold(readEulerToLagrang.lookupOrDefault<scalar>("EL_cellThreshold", 0));
Info << "Set the EL_cellThreshold: "  <<  EL_cellThreshold  << endl;

//- ML: Set threshold for tracking structures with alpha below this Threshold  
scalar EL_alphaThreshold(readEulerToLagrang.lookupOrDefault<scalar>("EL_alphaThreshold", 0));
Info << "Set the EL_alphaThreshold: "  <<  EL_alphaThreshold  << endl;
