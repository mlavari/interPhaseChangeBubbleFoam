/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "basicKinematicBubbleCollidingCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// KinematicBubble
#include "makeParcelForces.H"
#include "makeParcelDispersionModels.H"
#include "makeParcelInjectionModels.H"
#include "makeParcelCollisionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelStochasticCollisionModels.H"
#include "makeParcelSurfaceFilmModels.H"

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(basicKinematicBubbleCollidingCloud);

// KinematicBubble sub-models
makeParcelForces(basicKinematicBubbleCollidingCloud);
makeParcelDispersionModels(basicKinematicBubbleCollidingCloud);
makeParcelInjectionModels(basicKinematicBubbleCollidingCloud);
makeParcelCollisionModels(basicKinematicBubbleCollidingCloud);
makeParcelPatchInteractionModels(basicKinematicBubbleCollidingCloud);
makeParcelStochasticCollisionModels(basicKinematicBubbleCollidingCloud);
makeParcelSurfaceFilmModels(basicKinematicBubbleCollidingCloud);

// MPPIC sub-models
makeMPPICParcelDampingModels(basicKinematicBubbleCollidingCloud);
makeMPPICParcelIsotropyModels(basicKinematicBubbleCollidingCloud);
makeMPPICParcelPackingModels(basicKinematicBubbleCollidingCloud);

// ************************************************************************* //
