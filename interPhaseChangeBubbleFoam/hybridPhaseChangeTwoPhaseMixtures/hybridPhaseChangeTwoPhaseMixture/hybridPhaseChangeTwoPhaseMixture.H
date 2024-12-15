/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::hybridPhaseChangeTwoPhaseMixture

Description

SourceFiles
    hybridPhaseChangeTwoPhaseMixture.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef hybridPhaseChangeTwoPhaseMixture_H
#define hybridPhaseChangeTwoPhaseMixture_H

#include "hybridIncompressibleTwoPhaseMixture.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class hybridPhaseChangeTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class hybridPhaseChangeTwoPhaseMixture
:
    public hybridIncompressibleTwoPhaseMixture
{

protected:

    // Protected data

        dictionary hybridPhaseChangeTwoPhaseMixtureCoeffs_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;


    // Private Member Functions

        //- No copy construct
        hybridPhaseChangeTwoPhaseMixture(const hybridPhaseChangeTwoPhaseMixture&) = delete;

        //- No copy assignment
        void operator=(const hybridPhaseChangeTwoPhaseMixture&) = delete;


public:

    //- Runtime type information
    TypeName("hybridPhaseChangeTwoPhaseMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            hybridPhaseChangeTwoPhaseMixture,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<hybridPhaseChangeTwoPhaseMixture> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        hybridPhaseChangeTwoPhaseMixture
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~hybridPhaseChangeTwoPhaseMixture() = default;


    // Member Functions

        //- Return const-access to the saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField>> mDotP() const = 0;

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        Pair<tmp<volScalarField>> vDotAlphal() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients to multiply (p - pSat)
        Pair<tmp<volScalarField>> vDotP() const;

        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
