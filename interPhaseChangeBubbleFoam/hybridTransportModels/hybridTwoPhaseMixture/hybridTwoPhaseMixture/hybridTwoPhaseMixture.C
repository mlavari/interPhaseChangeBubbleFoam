/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "hybridTwoPhaseMixture.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hybridTwoPhaseMixture::hybridTwoPhaseMixture
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    phase1Name_(dict.get<wordList>("phases")[0]),
    phase2Name_(dict.get<wordList>("phases")[1]),
    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - alpha1_
    ),
    alpha1_BC_(alpha1_.boundaryField().types()),
    beta1_
    (
        IOobject
        (
            IOobject::groupName("beta", phase1Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1.0),
        alpha1_BC()
    ),
    beta2_
    (
        IOobject
        (
            IOobject::groupName("beta", phase2Name_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - beta1_
    ),
    gamma1_
    (
        IOobject
        (
            IOobject::groupName("gamma", phase1Name_),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1_ * beta1_,
        alpha1_BC()
    ),
    gamma2_
    (
        IOobject
        (
            IOobject::groupName("gamma", phase2Name_),
            mesh.time().timeName(),
            mesh
        ),
        1.0 - gamma1_
    )
{
}


// ************************************************************************* //
