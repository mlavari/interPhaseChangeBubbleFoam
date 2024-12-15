/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "hybridPhaseChangeTwoPhaseMixture.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hybridPhaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(hybridPhaseChangeTwoPhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hybridPhaseChangeTwoPhaseMixture::hybridPhaseChangeTwoPhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    hybridIncompressibleTwoPhaseMixture(U, phi),
    hybridPhaseChangeTwoPhaseMixtureCoeffs_(optionalSubDict(type + "Coeffs")),
    pSat_("pSat", dimPressure, *this)
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::hybridPhaseChangeTwoPhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/rho1() - alpha1_*(1.0/rho1() - 1.0/rho2()));
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::hybridPhaseChangeTwoPhaseMixture::vDotP() const
{
    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}


bool Foam::hybridPhaseChangeTwoPhaseMixture::read()
{
    if (hybridIncompressibleTwoPhaseMixture::read())
    {
        hybridPhaseChangeTwoPhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");
        readEntry("pSat", pSat_);

        return true;
    }

    return false;
}


// ************************************************************************* //
