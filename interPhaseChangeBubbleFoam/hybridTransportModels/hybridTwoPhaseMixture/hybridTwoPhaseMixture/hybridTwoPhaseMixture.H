/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::hybridTwoPhaseMixture

Description
    A two-phase mixture model

SourceFiles
    hybridTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef hybridTwoPhaseMixture_H
#define hybridTwoPhaseMixture_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hybridTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class hybridTwoPhaseMixture
{
protected:

    // Protected data

        word phase1Name_;
        word phase2Name_;
        //- The Eulerian phase-fraction
        volScalarField alpha1_;
        volScalarField alpha2_;
        //-ML: The Boundary Condition of alpha1
        wordList alpha1_BC_;
        //-ML: The Lagrangian phase-fraction
        volScalarField beta1_;
        volScalarField beta2_;
        //-ML: The Hybrid phase-fraction
        volScalarField gamma1_;
        volScalarField gamma2_;

public:

    // Constructors

        //- Construct from components
        hybridTwoPhaseMixture
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~hybridTwoPhaseMixture()
    {}


    // Member Functions

        const word& phase1Name() const
        {
            return phase1Name_;
        }

        const word& phase2Name() const
        {
            return phase2Name_;
        }

        //- Return Eulerian phase-fraction of phase 1
        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        //- Return Eulerian phase-fraction of phase 1
        volScalarField& alpha1()
        {
            return alpha1_;
        }

        //- Return Eulerian phase-fraction of phase 2
        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        //- Return Eulerian phase-fraction of phase 2
        volScalarField& alpha2()
        {
            return alpha2_;
        }

        //-ML: Return Boundary Condition
        const wordList& alpha1_BC() const
        {
            return alpha1_BC_;
        }

        //-ML: Return Boundary Condition
        wordList& alpha1_BC()
        {
            return alpha1_BC_;
        }

        //-ML: Return Lagrangian phase-fraction of phase 1
        const volScalarField& beta1() const
        {
            return beta1_;
        }

        //-ML: Return Lagrangian phase-fraction of phase 1
        volScalarField& beta1()
        {
            return beta1_;
        }

        //-ML: Return Lagrangian phase-fraction of phase 2
        const volScalarField& beta2() const
        {
            return beta2_;
        }

        //-ML: Return Lagrangian phase-fraction of phase 2
        volScalarField& beta2()
        {
            return beta2_;
        }

        //-ML: Return Hybrid phase-fraction of phase 1
        const volScalarField& gamma1() const
        {
            return gamma1_;
        }

        //-ML: Return Hybrid phase-fraction of phase 1
        volScalarField& gamma1()
        {
            return gamma1_;
        }

        //-ML: Return Hybrid phase-fraction of phase 2
        const volScalarField& gamma2() const
        {
            return gamma2_;
        }

        //-ML: Return Hybrid phase-fraction of phase 2
        volScalarField& gamma2()
        {
            return gamma2_;
        }

        //-ML: Return Hybrid phase-fraction of phase 2
        void update_gamma()
        {
            gamma1_ = alpha1_ * beta1_;
            gamma2_ = 1.0 - gamma1_;
        }
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
