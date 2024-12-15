/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SchnerrSauer.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace hybridPhaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(SchnerrSauer, 0);
    addToRunTimeSelectionTable
    (
        hybridPhaseChangeTwoPhaseMixture,
        SchnerrSauer,
        components
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hybridPhaseChangeTwoPhaseMixtures::SchnerrSauer::SchnerrSauer
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    hybridPhaseChangeTwoPhaseMixture(typeName, U, phi),

    n_("n", dimless/dimVolume, hybridPhaseChangeTwoPhaseMixtureCoeffs_),
    dNuc_("dNuc", dimLength, hybridPhaseChangeTwoPhaseMixtureCoeffs_),
    Cc_("Cc", dimless, hybridPhaseChangeTwoPhaseMixtureCoeffs_),
    Cv_("Cv", dimless, hybridPhaseChangeTwoPhaseMixtureCoeffs_),

    p0_(pSat().dimensions(), Zero)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::hybridPhaseChangeTwoPhaseMixtures::SchnerrSauer::rRb
(
    const volScalarField& limitedAlpha1
) const
{
    return pow
    (
        ((4*constant::mathematical::pi*n_)/3)
       *limitedAlpha1/(1.0 + alphaNuc() - limitedAlpha1),
        1.0/3.0
    );
}


Foam::dimensionedScalar
Foam::hybridPhaseChangeTwoPhaseMixtures::SchnerrSauer::alphaNuc() const
{
    dimensionedScalar Vnuc = n_*constant::mathematical::pi*pow3(dNuc_)/6;
    return Vnuc/(1 + Vnuc);
}


Foam::tmp<Foam::volScalarField>
Foam::hybridPhaseChangeTwoPhaseMixtures::SchnerrSauer::pCoeff
(
    const volScalarField& p
) const
{
    volScalarField limitedAlpha1(clamp(alpha1_, zero_one{}));
    volScalarField rho
    (
        limitedAlpha1*rho1() + (scalar(1) - limitedAlpha1)*rho2()
    );

    return
        (3*rho1()*rho2())*sqrt(2/(3*rho1()))
       *rRb(limitedAlpha1)/(rho*sqrt(mag(p - pSat()) + 0.01*pSat()));
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::hybridPhaseChangeTwoPhaseMixtures::SchnerrSauer::mDotAlphal() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(clamp(alpha1_, zero_one{}));

    return Pair<tmp<volScalarField>>
    (
        Cc_*limitedAlpha1*pCoeff*max(p - pSat(), p0_),

        Cv_*(1.0 + alphaNuc() - limitedAlpha1)*pCoeff*min(p - pSat(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::hybridPhaseChangeTwoPhaseMixtures::SchnerrSauer::mDotP() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(clamp(alpha1_, zero_one{}));
    volScalarField apCoeff(limitedAlpha1*pCoeff);

    return Pair<tmp<volScalarField>>
    (
        Cc_*(1.0 - limitedAlpha1)*pos0(p - pSat())*apCoeff,

        (-Cv_)*(1.0 + alphaNuc() - limitedAlpha1)*neg(p - pSat())*apCoeff
    );
}


void Foam::hybridPhaseChangeTwoPhaseMixtures::SchnerrSauer::correct()
{}


bool Foam::hybridPhaseChangeTwoPhaseMixtures::SchnerrSauer::read()
{
    if (hybridPhaseChangeTwoPhaseMixture::read())
    {
        hybridPhaseChangeTwoPhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");

        hybridPhaseChangeTwoPhaseMixtureCoeffs_.readEntry("n", n_);
        hybridPhaseChangeTwoPhaseMixtureCoeffs_.readEntry("dNuc", dNuc_);
        hybridPhaseChangeTwoPhaseMixtureCoeffs_.readEntry("Cc", Cc_);
        hybridPhaseChangeTwoPhaseMixtureCoeffs_.readEntry("Cv", Cv_);

        return true;
    }

    return false;
}


// ************************************************************************* //
