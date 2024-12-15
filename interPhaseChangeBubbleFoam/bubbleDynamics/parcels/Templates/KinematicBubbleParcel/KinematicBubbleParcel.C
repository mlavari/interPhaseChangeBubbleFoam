/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "KinematicBubbleParcel.H"
#include "forceSuSp.H"
#include "integrationScheme.H"
#include "meshTools.H"
#include "cloudSolution.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::label Foam::KinematicBubbleParcel<ParcelType>::maxTrackAttempts = 1;


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicBubbleParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar distance
)
{
    //-ML: Initialise finding fields
    vector posB= this->position();
    scalar R_sphere = distance * (this->d()/2); // Define the radius of the sphere
    //-ML: Number of divisions in theta (azimuthal angle) and phi (polar angle)
    label numTheta = 10;
    label numPhi = 10;
    vector USum(0, 0, 0);
    scalar rhoSum = 0.0;
    scalar muSum = 0.0;
    scalar pLSum = 0.0;
    label validSamples = 0; //-ML: Counter for valid samples

    //-ML: Loop over the spherical surface with uniform sampling
    for (label i = 0; i < numTheta; i++)
    {
        scalar thetaBubble = (2 * constant::mathematical::pi * i) / numTheta; //-ML: Azimuthal angle
        
        for (label j = 0; j < numPhi; j++)
        {
            scalar phiBubble = (constant::mathematical::pi * j) / (numPhi - 1); //-ML: Polar angle

            //-ML: Compute the unit vector components explicitly
            scalar x = sin(phiBubble) * cos(thetaBubble);
            scalar y = sin(phiBubble) * sin(thetaBubble);
            scalar z = cos(phiBubble);

            //-ML: Create the vector using the computed components
            vector unitVec(x, y, z);

            //-ML: Scale by the sphere radius and shift by the center
            vector samplePos = posB + R_sphere * unitVec;

            //-ML: Find the cell that contains this position
            label cellSample = cloud.mesh().findCell(samplePos);

            //-ML: Check if a valid cell was found
            if (cellSample != -1) 
            {
                //-ML: Interpolate values at the sample position
                vector USample = td.UInterp().interpolate(samplePos, cellSample);
                scalar rhoSample = td.rhoInterp().interpolate(samplePos, cellSample);
                scalar muSample = td.muInterp().interpolate(samplePos, cellSample);
                scalar pLSample = td.pInterp().interpolate(samplePos, cellSample);

                //-ML: Accumulate the values
                USum += USample;
                rhoSum += rhoSample;
                muSum += muSample;
                pLSum += pLSample;

                //-ML: Increment valid sample counter
                validSamples++;
            }
        }
    }

    //-ML: Calculate the average values over the surface
    scalar avgRho = (validSamples > 0) ? (rhoSum / validSamples) : 0.0;
    vector avgU = (validSamples > 0) ? (USum / validSamples) : Foam::Vector<double>::zero;
    scalar avgMu = (validSamples > 0) ? (muSum / validSamples) : 0.0;
    scalar avgPL = (validSamples > 0) ? (pLSum / validSamples) : 0.0;

    td.rhoc() = avgRho;
    if (td.rhoc() < cloud.constProps().rhoMin())
    {
        if (debug)
        {
            WarningInFunction
                << "Limiting observed density in cell " << this->cell()
                << " to " << cloud.constProps().rhoMin() <<  nl << endl;
        }

        td.rhoc() = cloud.constProps().rhoMin();
    }
    td.Uc() = avgU;
    td.muc() = avgMu;
    td.pc() = avgPL; 
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicBubbleParcel<ParcelType>::calcDispersion
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    td.Uc() = cloud.dispersion().update
    (
        dt,
        this->cell(),
        U_,
        td.Uc(),
        UTurb_,
        tTurb_
    );
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicBubbleParcel<ParcelType>::calcUCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typename TrackCloudType::parcelType& p =
        static_cast<typename TrackCloudType::parcelType&>(*this);

    this->UCorrect_ = Zero;

    this->UCorrect_ =
        cloud.dampingModel().velocityCorrection(p, dt);

    this->UCorrect_ +=
        cloud.packingModel().velocityCorrection(p, dt);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicBubbleParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    td.Uc() += cloud.UTrans()[this->cell()]/massCell(td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicBubbleParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = nParticle_;
    const scalar mass0 = mass();
    
    // Reynolds number
    const scalar Re = this->Re(td);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        calcVelocity(cloud, td, dt, Re, td.muc(), mass0, Su, dUTrans, Spu);

    this->U_ += this->UCorrect_;

    //-ML: Mass generation in cell due to RP equation
    scalar mDotB = (mass() - mass0);

    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (cloud.solution().coupled())
    {

    //-ML: Find position and radius of bubble
    vector posB = this->position();
    scalar radB = this->d() / 2;
    label startCell = this->cell();

    //-ML: DynamicList to store cells and their distances
    DynamicList<std::pair<label, scalar>> cellsWithDistances;

    //-ML: Add the start cell to the list 
    vector cellCenter = this->mesh().cellCentres()[startCell];
    cellsWithDistances.append(std::make_pair(startCell, mag(cellCenter - posB)));  //-ML: Add the first cell with distance

    //-ML: Create a set to avoid re-visiting cells
    labelHashSet visitedCells;
    visitedCells.insert(startCell);

    //-ML: Dynamic list to store all neighbors for further exploration
    DynamicList<label> cells(0);
    cells.append(this->mesh().cellCells()[startCell]);


    //-ML: Start outer loop to explore neighbors
    while (cells.size() > 0)
    {
        //-ML: New dynamic list to store newly found neighbors within the radius
        DynamicList<label> newNeighbours;

        //-ML: Loop over all current neighbors in the cells list
        forAll(cells, i)
        {
            label cellI = cells[i];
            if (!visitedCells.found(cellI))  //-ML: Check if cellI is already visited
            {
                cellCenter = this->mesh().cellCentres()[cellI];

                //-ML: Check if this neighbor cell's center is within the radius
                scalar distance = mag(cellCenter - posB);  //-ML: Calculate distance to posB

                if (distance <= (radB))
                {
                    cellsWithDistances.append(std::make_pair(cellI, distance));  //-ML: Store the cell and its distance
                    newNeighbours.append(this->mesh().cellCells()[cellI]);  //-ML: Add neighbors to the list
                }

                visitedCells.insert(cellI);  //-ML: Mark cell as visited
            }
        }

        //-ML: Sort the cells based on their distance to posB
        std::sort(cellsWithDistances.begin(), cellsWithDistances.end(),
                [](const std::pair<label, scalar>& a, const std::pair<label, scalar>& b)
                {
                    return a.second < b.second;  //-ML: Compare distances
                });

        //-ML: Continue to the next iteration with newNeighbours
        cells.clear();
        cells.transfer(newNeighbours);
    }

    scalar sigma = radB / this->LE_deviation();  //-ML: Gaussian standard deviation
    scalar totalWeightWithinRadius = 0.0;

    //-ML: Calculate the Gaussian weight for each cell and accumulate total weight within radius
    List<scalar> gaussianWeights(cellsWithDistances.size());

    forAll(cellsWithDistances, j) 
    {

        scalar gaussFactor;

        if (cellsWithDistances.size() == 1) {
            gaussFactor = 1.0;  //-ML: Set gaussFactor to 1 if there is only one cell
        } else {
            scalar distance = cellsWithDistances[j].second;  //-ML: Get distance from pair
            //-ML: Calculate Gaussian weight for this cell
            gaussFactor = exp(-0.5 * pow(distance / sigma, 2.0));
        }
        gaussianWeights[j] = gaussFactor;

        //-ML: Accumulate total weight
        totalWeightWithinRadius += gaussFactor;
    }

    //-ML: Distribute the bubble volume using normalized Gaussian weights
    forAll(cellsWithDistances, j) {

        label cellJ = cellsWithDistances[j].first;

        //-ML: Normalize the Gaussian weight and compute `vi_j` for this cell
        scalar normalizedWeight = gaussianWeights[j] / totalWeightWithinRadius;
        //-ML: Update momentum transfer
        cloud.UTrans()[cellJ] += normalizedWeight*np0*dUTrans;

        //-ML: Update momentum transfer coefficient
        cloud.UCoeff()[cellJ] += normalizedWeight*np0*Spu;

        //-ML: Update Return rate of change of mass
        cloud.mDotBubble()[cellJ] += normalizedWeight*np0*mDotB;

    }    
    }
}


template<class ParcelType>
template<class TrackCloudType>
const Foam::vector Foam::KinematicBubbleParcel<ParcelType>::calcVelocity
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar Re,
    const scalar mu,
    const scalar mass,
    const vector& Su,
    vector& dUTrans,
    scalar& Spu
)
{
    const typename TrackCloudType::parcelType& p =
        static_cast<const typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    const typename TrackCloudType::forceType& forces = cloud.forces();

    // Momentum source due to particle forces
    const forceSuSp Fcp = forces.calcCoupled(p, ttd, dt, mass, Re, mu); 
    const forceSuSp Fncp = forces.calcNonCoupled(p, ttd, dt, mass, Re, mu);
    const scalar massEff = forces.massEff(p, ttd, mass);//-ML: in case mass transfer is added

    /*
    // Proper splitting ...
    // Calculate the integration coefficients
    const vector acp = (Fcp.Sp()*td.Uc() + Fcp.Su())/massEff;
    const vector ancp = (Fncp.Sp()*td.Uc() + Fncp.Su() + Su)/massEff;
    const scalar bcp = Fcp.Sp()/massEff;
    const scalar bncp = Fncp.Sp()/massEff;

    // Integrate to find the new parcel velocity
    const vector deltaUcp =
        cloud.UIntegrator().partialDelta
        (
            U_, dt, acp + ancp, bcp + bncp, acp, bcp
        );
    const vector deltaUncp =
        cloud.UIntegrator().partialDelta
        (
            U_, dt, acp + ancp, bcp + bncp, ancp, bncp
        );
    const vector deltaT = deltaUcp + deltaUncp;
    */

    // Shortcut splitting assuming no implicit non-coupled force ...
    // Calculate the integration coefficients
    const vector acp = (Fcp.Sp()*td.Uc() + Fcp.Su())/massEff;
    const vector ancp = (Fncp.Su() + Su)/massEff;
    const scalar bcp = Fcp.Sp()/massEff;

    // Integrate to find the new parcel velocity
    const vector deltaU = cloud.UIntegrator().delta(U_, dt, acp + ancp, bcp); // -ML: deltaU = Unew - Uold
    const vector deltaUncp = ancp*dt;
    const vector deltaUcp = deltaU - deltaUncp;

    // Calculate the new velocity and the momentum transfer terms
    vector Unew = U_ + deltaU;

    //- ML: Include bubble dynamics
    if(this->bubble_activation())
    {
        #include "bubbleDynamics.H"
    }

    //dUTrans -= massEff*deltaUcp;
    vector Unew2 = U_ + deltaUcp;
    scalar mass2 = this->mass();
    const scalar massEff2 = forces.massEff(p, ttd, mass2);//-ML: In case virtual mass transfer is added
    dUTrans -= ((massEff2*Unew2)-(massEff*U_)); //-ML: Modify the source term to make it according to the equations

    Spu = dt*Fcp.Sp();

    // Apply correction to velocity and dUTrans for reduced-D cases
    const polyMesh& mesh = cloud.pMesh();
    meshTools::constrainDirection(mesh, mesh.solutionD(), Unew);
    meshTools::constrainDirection(mesh, mesh.solutionD(), dUTrans);

    return Unew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicBubbleParcel<ParcelType>::KinematicBubbleParcel
(
    const KinematicBubbleParcel<ParcelType>& p
)
:
    ParcelType(p),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    //-ML: Growth rate of radius [m/s]
    R_dot_(p.R_dot_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    UCorrect_(p.UCorrect_),
    //- ML: flag for activation of bubble dynamics
    bubble_activation_(p.bubble_activation_),
    //- ML: surface tension for bubble dynamics
    bubbleSigma_(p.bubbleSigma_),
    //- ML: initial pressure for bubble dynamics
    p0_(p.p0_),
    //- ML: initial radius for bubble dynamics            
    R0_(p.R0_),
    //- ML: vapour pressure for bubble dynamics            
    pv_(p.pv_),
    //- ML: time step for each loop of RP solver           
    RPdT_(p.RPdT_),
    //- ML: kappa isotropic value for bubble dynamics            
    bubbleKappa_(p.bubbleKappa_),
    //- ML: Averaging distance for bubble dynamics
    averagingDistance_(p.averagingDistance_),
    //- ML: ODE solver type for bubble dynamics
    ODESolverType_(p.ODESolverType_),
    //- ML: flag for activation of EulerianToLagrangian
    EulerianToLagrangian_activation_(p.EulerianToLagrangian_activation_),
    //- ML: flag for activation of LagrangianToEulerian
    LagrangianToEulerian_activation_(p.LagrangianToEulerian_activation_),
    //- ML:sigma = bubbleRadius / deviation
    LE_deviation_(p.LE_deviation_),        
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    LE_minCellOccupancy_(p.LE_minCellOccupancy_),          
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    LE_cellThreshold_(p.LE_cellThreshold_),          
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    LE_alphaThreshold_(p.LE_alphaThreshold_),         
    //- ML: Minimum bubble radius threshold that will be tracked
    LE_bubbleSizeThreshold_(p.LE_bubbleSizeThreshold_),         
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    LE_boxCheckEnabled_(p.LE_boxCheckEnabled_),
    //- ML: the top-left corner coordinates
    LE_boxTopLeftCorner_(p.LE_boxTopLeftCorner_),  
    //- ML: the bottom-right corner coordinates
    LE_boxBottomRightCorner_(p.LE_boxBottomRightCorner_)
{}


template<class ParcelType>
Foam::KinematicBubbleParcel<ParcelType>::KinematicBubbleParcel
(
    const KinematicBubbleParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    active_(p.active_),
    typeId_(p.typeId_),
    nParticle_(p.nParticle_),
    d_(p.d_),
    //- Growth rate of radius [m/s]
    R_dot_(p.R_dot_),
    dTarget_(p.dTarget_),
    U_(p.U_),
    rho_(p.rho_),
    age_(p.age_),
    tTurb_(p.tTurb_),
    UTurb_(p.UTurb_),
    UCorrect_(p.UCorrect_),
    //- ML: flag for activation of bubble dynamics
    bubble_activation_(p.bubble_activation_),
    //- ML: surface tension for bubble dynamics
    bubbleSigma_(p.bubbleSigma_),
    //- ML: initial pressure for bubble dynamics
    p0_(p.p0_),
    //- ML: initial radius for bubble dynamics            
    R0_(p.R0_),
    //- ML: vapour pressure for bubble dynamics            
    pv_(p.pv_),   
    //- ML: time step for each loop of RP solver           
    RPdT_(p.RPdT_),
    //- ML: kappa isotropic value for bubble dynamics            
    bubbleKappa_(p.bubbleKappa_),
    //- ML: Averaging distance for bubble dynamics
    averagingDistance_(p.averagingDistance_),
    //- ML: ODE solver type for bubble dynamics
    ODESolverType_(p.ODESolverType_),
    //- ML: flag for activation of EulerianToLagrangian
    EulerianToLagrangian_activation_(p.EulerianToLagrangian_activation_),
    //- ML: flag for activation of LagrangianToEulerian
    LagrangianToEulerian_activation_(p.LagrangianToEulerian_activation_),
    //- ML:sigma = bubbleRadius / deviation
    LE_deviation_(p.LE_deviation_),        
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    LE_minCellOccupancy_(p.LE_minCellOccupancy_),          
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    LE_cellThreshold_(p.LE_cellThreshold_),          
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    LE_alphaThreshold_(p.LE_alphaThreshold_),         
    //- ML: Minimum bubble radius threshold that will be tracked
    LE_bubbleSizeThreshold_(p.LE_bubbleSizeThreshold_),         
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    LE_boxCheckEnabled_(p.LE_boxCheckEnabled_),
    //- ML: the top-left corner coordinates
    LE_boxTopLeftCorner_(p.LE_boxTopLeftCorner_),  
    //- ML: the bottom-right corner coordinates
    LE_boxBottomRightCorner_(p.LE_boxBottomRightCorner_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
bool Foam::KinematicBubbleParcel<ParcelType>::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    auto& p = static_cast<typename TrackCloudType::parcelType&>(*this);
    auto& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    ttd.switchProcessor = false;
    ttd.keepParticle = true;

    const cloudSolution& solution = cloud.solution();
    const scalarField& cellLengthScale = cloud.cellLengthScale();
    const scalar maxCo = solution.maxCo();

    //- ML: Reinitialize the R_dot before calculating the bubble
    if(this->bubble_activation())
    {
        p.R_dot_ = 0;
    }

    while (ttd.keepParticle && !ttd.switchProcessor && p.stepFraction() < 1)
    {
        // Cache the current position, cell and step-fraction
        const point start = p.position();
        const scalar sfrac = p.stepFraction();

        // Total displacement over the time-step
        const vector s = trackTime*U_;

        // Cell length scale
        const scalar l = cellLengthScale[p.cell()];

        // Deviation from the mesh centre for reduced-D cases
        const vector d = p.deviationFromMeshCentre();

        // Fraction of the displacement to track in this loop. This is limited
        // to ensure that the both the time and distance tracked is less than
        // maxCo times the total value.
        scalar f = 1 - p.stepFraction();
        f = min(f, maxCo);
        f = min(f, maxCo*l/max(SMALL*l, mag(s)));
        if (p.active())
        {
            // Track to the next face
            p.trackToFace(f*s - d, f);
        }
        else
        {
            // At present the only thing that sets active_ to false is a stick
            // wall interaction. We want the position of the particle to remain
            // the same relative to the face that it is on. The local
            // coordinates therefore do not change. We still advance in time and
            // perform the relevant interactions with the fixed particle.
            p.stepFraction() += f;
        }

        const scalar dt = (p.stepFraction() - sfrac)*trackTime;

        // Avoid problems with extremely small timesteps
        if (dt > ROOTVSMALL)
        {
            // Update cell based properties
            p.setCellValues(cloud, ttd);

            p.calcDispersion(cloud, ttd, dt);

            if (solution.cellValueSourceCorrection())
            {
                p.cellValueSourceCorrection(cloud, ttd, dt);
            }

            p.calcUCorrection(cloud, ttd, dt);

            p.calc(cloud, ttd, dt); //-ML: new velocity is calculated based on the forces and source terms
        }

        p.age() += dt;

        if (p.active() && p.onFace())
        {
            ttd.keepParticle = cloud.functions().postFace(p, ttd);
        }

        ttd.keepParticle = cloud.functions().postMove(p, dt, start, ttd);

        if (p.active() && p.onFace() && ttd.keepParticle)
        {
            p.hitFace(s, cloud, ttd);
        }
    }

    //- ML: Include transition from Lagrangian to Eulerian
    if(this->LagrangianToEulerian_activation())
    {
        #include "LagrangianToEulerian.H"
    }

    return ttd.keepParticle;
}


template<class ParcelType>
template<class TrackCloudType>
bool Foam::KinematicBubbleParcel<ParcelType>::hitPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    auto& p = static_cast<typename TrackCloudType::parcelType&>(*this);
    auto& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);

    const polyPatch& pp = p.mesh().boundaryMesh()[p.patch()];

    // Invoke post-processing model
    td.keepParticle = cloud.functions().postPatch(p, pp, ttd);

    if (isA<processorPolyPatch>(pp))
    {
        // Skip processor patches
        return false;
    }
    else if (cloud.surfaceFilm().transferParcel(p, pp, td.keepParticle))
    {
        // Surface film model consumes the interaction, i.e. all
        // interactions done
        return true;
    }
    else
    {
        // This does not take into account the wall interaction model
        // Just the polyPatch type. Then, a patch type which has 'rebound'
        // interaction model will count as escaped parcel while it is not
        if (!isA<wallPolyPatch>(pp) && !polyPatch::constraintType(pp.type()))
        {
            cloud.patchInteraction().addToEscapedParcels(nParticle_*mass());
        }

        // Invoke patch interaction model
        return cloud.patchInteraction().correct(p, pp, td.keepParticle);
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicBubbleParcel<ParcelType>::hitProcessorPatch
(
    TrackCloudType&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::KinematicBubbleParcel<ParcelType>::hitWallPatch
(
    TrackCloudType&,
    trackingData&
)
{
    // wall interactions are handled by the generic hitPatch method
}


template<class ParcelType>
void Foam::KinematicBubbleParcel<ParcelType>::transformProperties(const tensor& T)
{
    ParcelType::transformProperties(T);

    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::KinematicBubbleParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    ParcelType::transformProperties(separation);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "KinematicBubbleParcelIO.C"

// ************************************************************************* //
