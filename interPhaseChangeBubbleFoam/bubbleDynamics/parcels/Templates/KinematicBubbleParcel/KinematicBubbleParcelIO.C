/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::KinematicBubbleParcel<ParcelType>::propertyList_ =
    Foam::KinematicBubbleParcel<ParcelType>::propertyList();


template<class ParcelType>
const std::size_t Foam::KinematicBubbleParcel<ParcelType>::sizeofFields
(
    sizeof(KinematicBubbleParcel<ParcelType>)
  - offsetof(KinematicBubbleParcel<ParcelType>, active_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicBubbleParcel<ParcelType>::KinematicBubbleParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    active_(false),
    typeId_(0),
    nParticle_(0.0),
    d_(0.0),
    //- ML: Growth rate of radius [m/s]
    R_dot_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    UCorrect_(Zero),
    //- ML: flag for activation of bubble dynamics
    bubble_activation_(false),
    //- ML: surface tension for bubble dynamics
    bubbleSigma_(0.0),
    //- ML: initial pressure for bubble dynamics
    p0_(0.0),
    //- ML: initial radius for bubble dynamics            
    R0_(0.0),
    //- ML: vapour pressure for bubble dynamics            
    pv_(0.0),
    //- ML: time step for each loop of RP solver           
    RPdT_(0.0),
    //- ML: kappa isotropic value for bubble dynamics            
    bubbleKappa_(0.0),
    //- ML: Averaging distance for bubble dynamics
    averagingDistance_(0.0),
    //- ML: ODE solver type for bubble dynamics
    ODESolverType_("defaultSolver"),
    //- ML: flag for activation of EulerianToLagrangian
    EulerianToLagrangian_activation_(false),
    //- ML: flag for activation of LagrangianToEulerian
    LagrangianToEulerian_activation_(false),
    //- ML:sigma = bubbleRadius / deviation
    LE_deviation_(0.0),        
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    LE_minCellOccupancy_(0.0),          
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    LE_cellThreshold_(0.0),          
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    LE_alphaThreshold_(0.0),         
    //- ML: Minimum bubble radius threshold that will be tracked
    LE_bubbleSizeThreshold_(0.0),         
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    LE_boxCheckEnabled_(false),
    //- ML: the top-left corner coordinates
    LE_boxTopLeftCorner_(Zero),  
    //- ML: the bottom-right corner coordinates
    LE_boxBottomRightCorner_(Zero)
{
    //- ML: ODE solver type for bubble dynamics from dictionary
    IOdictionary readBubbleProperties_
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
    dictionary readODEProperties_(readBubbleProperties_.subDict("bubbleProperties"));
    word test_(readODEProperties_.lookupOrDefault<word>("ODESolverType", "defaultSolver"));
    ODESolverType(test_);

    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is  >> active_
                >> typeId_
                >> nParticle_
                >> d_
                //- ML: Growth rate of radius [m/s]
                >> R_dot_
                >> dTarget_
                >> U_
                >> rho_
                >> age_
                >> tTurb_
                >> UTurb_
                >> UCorrect_
                //- ML: flag for activation of bubble dynamics
                >> bubble_activation_
                //- ML: surface tension for bubble dynamics
                >> bubbleSigma_
                //- ML: initial pressure for bubble dynamics
                >> p0_
                //- ML: initial radius for bubble dynamics            
                >> R0_
                //- ML: vapour pressure for bubble dynamics            
                >> pv_
                //- ML: time step for each loop of RP solver          
                >> RPdT_
                //- ML: kappa isotropic value for bubble dynamics            
                >> bubbleKappa_
                //- ML: Averaging distance for bubble dynamics
                >> averagingDistance_
                //- ML: ODE solver type for bubble dynamics
                //>> ODESolverType_;
                //- ML: flag for activation of EulerianToLagrangian
                >> EulerianToLagrangian_activation_
                //- ML: flag for activation of LagrangianToEulerian
                >> LagrangianToEulerian_activation_
                //- ML:sigma = bubbleRadius / deviation
                >> LE_deviation_
                //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
                >> LE_minCellOccupancy_   
                //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
                >> LE_cellThreshold_        
                //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
                >> LE_alphaThreshold_ 
                //- ML: Minimum bubble radius threshold that will be tracked
                >> LE_bubbleSizeThreshold_        
                //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
                >> LE_boxCheckEnabled_
                //- ML: the top-left corner coordinates
                >> LE_boxTopLeftCorner_
                //- ML: the bottom-right corner coordinates
                >> LE_boxBottomRightCorner_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawLabel(is, &active_);
            readRawLabel(is, &typeId_);
            readRawScalar(is, &nParticle_);
            readRawScalar(is, &d_);
            //- ML: Growth rate of radius [m/s]
            readRawScalar(is, &R_dot_);
            readRawScalar(is, &dTarget_);
            readRawScalar(is, U_.data(), vector::nComponents);
            readRawScalar(is, &rho_);
            readRawScalar(is, &age_);
            readRawScalar(is, &tTurb_);
            readRawScalar(is, UTurb_.data(), vector::nComponents);
            readRawScalar(is, UCorrect_.data(), vector::nComponents);
            //- ML: flag for activation of bubble dynamics
            readRawLabel(is, &bubble_activation_);
            //- ML: surface tension for bubble dynamics
            readRawScalar(is, &bubbleSigma_);
            //- ML: initial pressure for bubble dynamics
            readRawScalar(is, &p0_);
            //- ML: initial radius for bubble dynamics            
            readRawScalar(is, &R0_);
            //- ML: vapour pressure for bubble dynamics            
            readRawScalar(is, &pv_);
            //- ML: time step for each loop of RP solver            
            readRawScalar(is, &RPdT_);
            //- ML: kappa isotropic value for bubble dynamics            
            readRawScalar(is, &bubbleKappa_);
            //- ML: Averaging distance for bubble dynamics
            readRawScalar(is, &averagingDistance_);
            //- ML: ODE solver type for bubble dynamics
            //readRaw(is, &ODESolverType_);
            //- ML: flag for activation of EulerianToLagrangian
            readRawLabel(is, &EulerianToLagrangian_activation_);
            //- ML: flag for activation of LagrangianToEulerian
            readRawLabel(is, &LagrangianToEulerian_activation_);
            //- ML:sigma = bubbleRadius / deviation
            readRawScalar(is, &LE_deviation_);        
            //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
            readRawScalar(is, &LE_minCellOccupancy_);          
            //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
            readRawScalar(is, &LE_cellThreshold_);          
            //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
            readRawScalar(is, &LE_alphaThreshold_);         
            //- ML: Minimum bubble radius threshold that will be tracked
            readRawScalar(is, &LE_bubbleSizeThreshold_);         
            //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
            readRawLabel(is, &LE_boxCheckEnabled_);
            //- ML: the top-left corner coordinates
            readRawScalar(is, LE_boxTopLeftCorner_.data(), vector::nComponents); 
            //- ML: the bottom-right corner coordinates
            readRawScalar(is, LE_boxBottomRightCorner_.data(), vector::nComponents);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&active_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicBubbleParcel<ParcelType>::readFields(CloudType& c)
{
    const bool readOnProc = c.size();

    ParcelType::readFields(c);

    IOField<label> active
    (
        c.newIOobject("active", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, active);

    IOField<label> typeId
    (
        c.newIOobject("typeId", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, typeId);

    IOField<scalar> nParticle
    (
        c.newIOobject("nParticle", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d
    (
        c.newIOobject("d", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, d);

    //- ML: Growth rate of radius [m/s]
    IOField<scalar> R_dot
    (
        c.newIOobject("R_dot", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, R_dot);

    IOField<scalar> dTarget
    (
        c.newIOobject("dTarget", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, dTarget);

    IOField<vector> U
    (
        c.newIOobject("U", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, U);

    IOField<scalar> rho
    (
        c.newIOobject("rho", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age
    (
        c.newIOobject("age", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, age);

    IOField<scalar> tTurb
    (
        c.newIOobject("tTurb", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb
    (
        c.newIOobject("UTurb", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, UTurb);

    IOField<vector> UCorrect
    (
        c.newIOobject("UCorrect", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, UCorrect);

    //- ML: flag for activation of bubble dynamics
    IOField<label> bubble_activation
    (
        c.newIOobject("bubble_activation", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, bubble_activation);

    //- ML: surface tension for bubble dynamics
    IOField<scalar> bubbleSigma
    (
        c.newIOobject("bubbleSigma", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, bubbleSigma);

    //- ML: initial pressure for bubble dynamics
    IOField<scalar> p0
    (
        c.newIOobject("p0", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, p0);

    //- ML: initial radius for bubble dynamics 
    IOField<scalar> R0
    (
        c.newIOobject("R0", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, R0);    

    //- ML: vapour pressure for bubble dynamics   
    IOField<scalar> pv
    (
        c.newIOobject("pv", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, pv);        

    //- ML: time step for each loop of RP solver
    IOField<scalar> RPdT
    (
        c.newIOobject("RPdT", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, RPdT);           

    //- ML: kappa isotropic value for bubble dynamics   
    IOField<scalar> bubbleKappa
    (
        c.newIOobject("bubbleKappa", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, bubbleKappa);         

    //- ML: Averaging distance for bubble dynamics
    IOField<scalar> averagingDistance
    (
        c.newIOobject("averagingDistance", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, averagingDistance);     

    //- ML: ODE solver type for bubble dynamics
    /*
    IOField<word> ODESolverType
    (
        c.newIOobject("ODESolverType", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, ODESolverType);   
    */

    //- ML: flag for activation of EulerianToLagrangian
    IOField<label> EulerianToLagrangian_activation
    (
        c.newIOobject("EulerianToLagrangian_activation", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, EulerianToLagrangian_activation);

    //- ML: flag for activation of LagrangianToEulerian
    IOField<label> LagrangianToEulerian_activation
    (
        c.newIOobject("LagrangianToEulerian_activation", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, LagrangianToEulerian_activation);

    //- ML: sigma = bubbleRadius / deviation
    IOField<scalar> LE_deviation
    (
        c.newIOobject("LE_deviation", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, LE_deviation);           
    //- ML: Minimum liquid volume fraction value that lagrangian cell can occupy
    IOField<scalar> LE_minCellOccupancy
    (
        c.newIOobject("LE_minCellOccupancy", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, LE_minCellOccupancy);             
    //- ML: Number of cells that trigger the transtion from Lagrangian to Eulerian
    IOField<scalar> LE_cellThreshold
    (
        c.newIOobject("LE_cellThreshold", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, LE_cellThreshold);                 
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    IOField<scalar> LE_alphaThreshold
    (
        c.newIOobject("LE_alphaThreshold", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, LE_alphaThreshold);              
    //- ML: Minimum bubble radius threshold that will be tracked
    IOField<scalar> LE_bubbleSizeThreshold
    (
        c.newIOobject("LE_bubbleSizeThreshold", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, LE_bubbleSizeThreshold);        
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    IOField<label> LE_boxCheckEnabled
    (
        c.newIOobject("LE_boxCheckEnabled", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, LE_boxCheckEnabled);
    //- ML: the top-left corner coordinates
    IOField<vector> LE_boxTopLeftCorner
    (
        c.newIOobject("LE_boxTopLeftCorner", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, LE_boxTopLeftCorner); 
    //- ML: the bottom-right corner coordinates
    IOField<vector> LE_boxBottomRightCorner
    (
        c.newIOobject("LE_boxBottomRightCorner", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, LE_boxBottomRightCorner);


    label i = 0;

    for (KinematicBubbleParcel<ParcelType>& p : c)
    {
        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        //- ML: Growth rate of radius [m/s]
        p.R_dot_ = R_dot[i];
        p.dTarget_ = dTarget[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];
        p.UCorrect_ = UCorrect[i];
        //- ML: flag for activation of bubble dynamics
        p.bubble_activation_ = bubble_activation[i];
        //- ML: surface tension for bubble dynamics
        p.bubbleSigma_ = bubbleSigma[i];
        //- ML: initial pressure for bubble dynamics
        p.p0_ = p0[i];
        //- ML: initial radius for bubble dynamics            
        p.R0_ = R0[i];
        //- ML: vapour pressure for bubble dynamics            
        p.pv_ = pv[i];
        //- ML: time step for each loop of RP solver           
        p.RPdT_ = RPdT[i];
        //- ML: kappa isotropic value for bubble dynamics            
        p.bubbleKappa_ = bubbleKappa[i];
        //- ML: Averaging distance for bubble dynamics
        p.averagingDistance_ = averagingDistance[i];
        //- ML: ODE solver type for bubble dynamics
        //p.ODESolverType_ = ODESolverType[i];
        //- ML: flag for activation of EulerianToLagrangian
        p.EulerianToLagrangian_activation_ = EulerianToLagrangian_activation[i];
        //- ML: flag for activation of LagrangianToEulerian
        p.LagrangianToEulerian_activation_ = LagrangianToEulerian_activation[i];
        //- ML:sigma = bubbleRadius / deviation
        p.LE_deviation_ = LE_deviation[i];       
        //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
        p.LE_minCellOccupancy_ = LE_minCellOccupancy[i];        
        //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
        p.LE_cellThreshold_ = LE_cellThreshold[i];        
        //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
        p.LE_alphaThreshold_ = LE_alphaThreshold[i];       
        //- ML: Minimum bubble radius threshold that will be tracked
        p.LE_bubbleSizeThreshold_ = LE_bubbleSizeThreshold[i];         
        //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
        p.LE_boxCheckEnabled_ = LE_boxCheckEnabled[i];
        //- ML: the top-left corner coordinates
        p.LE_boxTopLeftCorner_ = LE_boxTopLeftCorner[i]; 
        //- ML: the bottom-right corner coordinates
        p.LE_boxBottomRightCorner_ = LE_boxBottomRightCorner[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicBubbleParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    const label np = c.size();
    const bool writeOnProc = c.size();

    IOField<label> active(c.newIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.newIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> nParticle
    (
        c.newIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.newIOobject("d", IOobject::NO_READ), np);
    //- ML: Growth rate of radius [m/s]
    IOField<scalar> R_dot(c.newIOobject("R_dot", IOobject::NO_READ), np);
    IOField<scalar> dTarget(c.newIOobject("dTarget", IOobject::NO_READ), np);
    IOField<vector> U(c.newIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> rho(c.newIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.newIOobject("age", IOobject::NO_READ), np);
    IOField<scalar> tTurb(c.newIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.newIOobject("UTurb", IOobject::NO_READ), np);
    IOField<vector> UCorrect(c.newIOobject("UCorrect", IOobject::NO_READ), np);
    //- ML: flag for activation of bubble dynamics
    IOField<label> bubble_activation(c.newIOobject("bubble_activation", IOobject::NO_READ), np);
    //- ML: surface tension for bubble dynamics
    IOField<scalar> bubbleSigma(c.newIOobject("bubbleSigma", IOobject::NO_READ), np);
    //- ML: initial pressure for bubble dynamics
    IOField<scalar> p0(c.newIOobject("p0", IOobject::NO_READ), np);
    //- ML: initial radius for bubble dynamics            
    IOField<scalar> R0(c.newIOobject("R0", IOobject::NO_READ), np);
    //- ML: vapour pressure for bubble dynamics            
    IOField<scalar> pv(c.newIOobject("pv", IOobject::NO_READ), np);
    //- ML: time step for each loop of RP solver            
    IOField<scalar> RPdT(c.newIOobject("RPdT", IOobject::NO_READ), np);
    //- ML: kappa isotropic value for bubble dynamics            
    IOField<scalar> bubbleKappa(c.newIOobject("bubbleKappa", IOobject::NO_READ), np);
    //- ML: Averaging distance for bubble dynamics
    IOField<scalar> averagingDistance(c.newIOobject("averagingDistance", IOobject::NO_READ), np);
    //- ML: ODE solver type for bubble dynamics
    //IOField<word> ODESolverType(c.newIOobject("ODESolverType", IOobject::NO_READ), np);
    //- ML: flag for activation of EulerianToLagrangian
    IOField<label> EulerianToLagrangian_activation(c.newIOobject("EulerianToLagrangian_activation", IOobject::NO_READ), np);
    //- ML: flag for activation of LagrangianToEulerian
    IOField<label> LagrangianToEulerian_activation(c.newIOobject("LagrangianToEulerian_activation", IOobject::NO_READ), np);
    //- ML:sigma = bubbleRadius / deviation
    IOField<scalar> LE_deviation(c.newIOobject("LE_deviation", IOobject::NO_READ), np);     
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    IOField<scalar> LE_minCellOccupancy(c.newIOobject("LE_minCellOccupancy", IOobject::NO_READ), np);        
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    IOField<scalar> LE_cellThreshold(c.newIOobject("LE_cellThreshold", IOobject::NO_READ), np);           
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    IOField<scalar> LE_alphaThreshold(c.newIOobject("LE_alphaThreshold", IOobject::NO_READ), np);         
    //- ML: Minimum bubble radius threshold that will be tracked
    IOField<scalar> LE_bubbleSizeThreshold(c.newIOobject("LE_bubbleSizeThreshold", IOobject::NO_READ), np);           
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    IOField<label> LE_boxCheckEnabled(c.newIOobject("LE_boxCheckEnabled", IOobject::NO_READ), np);  
    //- ML: the top-left corner coordinates
    IOField<vector> LE_boxTopLeftCorner(c.newIOobject("LE_boxTopLeftCorner", IOobject::NO_READ), np);   
    //- ML: the bottom-right corner coordinates
    IOField<vector> LE_boxBottomRightCorner(c.newIOobject("LE_boxBottomRightCorner", IOobject::NO_READ), np);  

    label i = 0;

    for (const KinematicBubbleParcel<ParcelType>& p : c)
    {
        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        //- ML: Growth rate of radius [m/s]
        R_dot[i] = p.R_dot();
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();
        UCorrect[i] = p.UCorrect();
        //- ML: flag for activation of bubble dynamics
        bubble_activation[i] = p.bubble_activation();
        //- ML: surface tension for bubble dynamics
        bubbleSigma[i] = p.bubbleSigma();
        //- ML: initial pressure for bubble dynamics
        p0[i] = p.p0();
        //- ML: initial radius for bubble dynamics            
        R0[i] = p.R0();
        //- ML: vapour pressure for bubble dynamics            
        pv[i] = p.pv();
        //- ML: time step for each loop of RP solver         
        RPdT[i] = p.RPdT();
        //- ML: kappa isotropic value for bubble dynamics            
        bubbleKappa[i] = p.bubbleKappa();
        //- ML: Averaging distance for bubble dynamics
        averagingDistance[i] = p.averagingDistance();
        //- ML: ODE solver type for bubble dynamics
        //ODESolverType[i] = p.ODESolverType();
        //- ML: flag for activation of EulerianToLagrangian
        EulerianToLagrangian_activation[i] = p.EulerianToLagrangian_activation();
        //- ML: flag for activation of LagrangianToEulerian
        LagrangianToEulerian_activation[i] = p.LagrangianToEulerian_activation();
        //- ML:sigma = bubbleRadius / deviation
        LE_deviation[i] = p.LE_deviation();        
        //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
        LE_minCellOccupancy[i] = p.LE_minCellOccupancy();          
        //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
        LE_cellThreshold[i] = p.LE_cellThreshold();         
        //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
        LE_alphaThreshold[i] = p.LE_alphaThreshold();         
        //- ML: Minimum bubble radius threshold that will be tracked
        LE_bubbleSizeThreshold[i] = p.LE_bubbleSizeThreshold();         
        //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
        LE_boxCheckEnabled[i] = p.LE_boxCheckEnabled();
        //- ML: the top-left corner coordinates
        LE_boxTopLeftCorner[i] = p.LE_boxTopLeftCorner();  
        //- ML: the bottom-right corner coordinates
        LE_boxBottomRightCorner[i] = p.LE_boxBottomRightCorner();

        ++i;
    }

    active.write(writeOnProc);
    typeId.write(writeOnProc);
    nParticle.write(writeOnProc);
    d.write(writeOnProc);
    //- ML: Growth rate of radius [m/s]
    R_dot.write(writeOnProc);
    dTarget.write(writeOnProc);
    U.write(writeOnProc);
    rho.write(writeOnProc);
    age.write(writeOnProc);
    tTurb.write(writeOnProc);
    UTurb.write(writeOnProc);
    UCorrect.write(writeOnProc);
    //- ML: flag for activation of bubble dynamics
    bubble_activation.write(writeOnProc);
    //- ML: surface tension for bubble dynamics
    bubbleSigma.write(writeOnProc);
    //- ML: initial pressure for bubble dynamics
    p0.write(writeOnProc);
    //- ML: initial radius for bubble dynamics            
    R0.write(writeOnProc);
    //- ML: vapour pressure for bubble dynamics            
    pv.write(writeOnProc);
    //- ML: time step for each loop of RP solver         
    RPdT.write(writeOnProc);
    //- ML: kappa isotropic value for bubble dynamics            
    bubbleKappa.write(writeOnProc);
    //- ML: Averaging distance for bubble dynamics
    averagingDistance.write(writeOnProc);
    //- ML: ODE solver type for bubble dynamics
    //ODESolverType.write(writeOnProc);
    //- ML: flag for activation of EulerianToLagrangian
    EulerianToLagrangian_activation.write(writeOnProc);
    //- ML: flag for activation of LagrangianToEulerian
    LagrangianToEulerian_activation.write(writeOnProc);
    //- ML:sigma = bubbleRadius / deviation
    LE_deviation.write(writeOnProc);    
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    LE_minCellOccupancy.write(writeOnProc);         
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    LE_cellThreshold.write(writeOnProc);        
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    LE_alphaThreshold.write(writeOnProc);      
    //- ML: Minimum bubble radius threshold that will be tracked
    LE_bubbleSizeThreshold.write(writeOnProc);      
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    LE_boxCheckEnabled.write(writeOnProc);
    //- ML: the top-left corner coordinates
    LE_boxTopLeftCorner.write(writeOnProc);
    //- ML: the bottom-right corner coordinates
    LE_boxBottomRightCorner.write(writeOnProc);

}


template<class ParcelType>
void Foam::KinematicBubbleParcel<ParcelType>::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    ParcelType::writeProperties(os, filters, delim, namesOnly);

    #undef  writeProp
    #define writeProp(Name, Value)                                            \
        ParcelType::writeProperty(os, Name, Value, namesOnly, delim, filters)

    writeProp("active", active_);
    writeProp("typeId", typeId_);
    writeProp("nParticle", nParticle_);
    writeProp("d", d_);
    //- ML: Growth rate of radius [m/s]
    writeProp("R_dot", R_dot_);
    writeProp("dTarget", dTarget_);
    writeProp("U", U_);
    writeProp("rho", rho_);
    writeProp("age", age_);
    writeProp("tTurb", tTurb_);
    writeProp("UTurb", UTurb_);
    writeProp("UCorrect", UCorrect_);
    //- ML: flag for activation of bubble dynamics
    writeProp("bubble_activation", bubble_activation_);
    //- ML: surface tension for bubble dynamics
    writeProp("bubbleSigma", bubbleSigma_);
    //- ML: initial pressure for bubble dynamics
    writeProp("p0", p0_);
    //- ML: initial radius for bubble dynamics            
    writeProp("R0", R0_);
    //- ML: vapour pressure for bubble dynamics            
    writeProp("pv", pv_);
    //- ML: time step for each loop of RP solver          
    writeProp("RPdT", RPdT_);
    //- ML: kappa isotropic value for bubble dynamics            
    writeProp("bubbleKappa", bubbleKappa_);
    //- ML: Averaging distance for bubble dynamics
    writeProp("averagingDistance", averagingDistance_);
    //- ML: ODE solver type for bubble dynamics
    //writeProp("ODESolverType", ODESolverType_);
    //- ML: flag for activation of EulerianToLagrangian
    writeProp("EulerianToLagrangian_activation", EulerianToLagrangian_activation_);
    //- ML: flag for activation of LagrangianToEulerian
    writeProp("LagrangianToEulerian_activation", LagrangianToEulerian_activation_);
    //- ML:sigma = bubbleRadius / deviation
    writeProp("LE_deviation", LE_deviation_);       
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    writeProp("LE_minCellOccupancy", LE_minCellOccupancy_);          
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    writeProp("LE_cellThreshold", LE_cellThreshold_);          
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    writeProp("LE_alphaThreshold", LE_alphaThreshold_);          
    //- ML: Minimum bubble radius threshold that will be tracked
    writeProp("LE_bubbleSizeThreshold", LE_bubbleSizeThreshold_);        
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    writeProp("LE_boxCheckEnabled", LE_boxCheckEnabled_);
    //- ML: the top-left corner coordinates
    writeProp("LE_boxTopLeftCorner", LE_boxTopLeftCorner_); 
    //- ML: the bottom-right corner coordinates
    writeProp("LE_boxBottomRightCorner", LE_boxBottomRightCorner_);

    #undef writeProp
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicBubbleParcel<ParcelType>::readObjects
(
    CloudType& c,
    const objectRegistry& obr
)
{
    ParcelType::readObjects(c, obr);

    if (!c.size()) return;

    const auto& active = cloud::lookupIOField<label>("active", obr);
    const auto& typeId = cloud::lookupIOField<label>("typeId", obr);
    const auto& nParticle = cloud::lookupIOField<scalar>("nParticle", obr);
    const auto& d = cloud::lookupIOField<scalar>("d", obr);
    //- ML: Growth rate of radius [m/s]
    const auto& R_dot = cloud::lookupIOField<scalar>("R_dot", obr);
    const auto& dTarget = cloud::lookupIOField<scalar>("dTarget", obr);
    const auto& U = cloud::lookupIOField<vector>("U", obr);
    const auto& rho = cloud::lookupIOField<scalar>("rho", obr);
    const auto& age = cloud::lookupIOField<scalar>("age", obr);
    const auto& tTurb = cloud::lookupIOField<scalar>("tTurb", obr);
    const auto& UTurb = cloud::lookupIOField<vector>("UTurb", obr);
    const auto& UCorrect = cloud::lookupIOField<vector>("UCorrect", obr);
    //- ML: flag for activation of bubble dynamics
    const auto& bubble_activation = cloud::lookupIOField<label>("bubble_activation", obr);
    //- ML: surface tension for bubble dynamics
    const auto& bubbleSigma = cloud::lookupIOField<scalar>("bubbleSigma", obr);
    //- ML: initial pressure for bubble dynamics
    const auto& p0 = cloud::lookupIOField<scalar>("p0", obr);
    //- ML: initial radius for bubble dynamics            
    const auto& R0 = cloud::lookupIOField<scalar>("R0", obr);
    //- ML: vapour pressure for bubble dynamics            
    const auto& pv = cloud::lookupIOField<scalar>("pv", obr);
    //- ML: time step for each loop of RP solver           
    const auto& RPdT = cloud::lookupIOField<scalar>("RPdT", obr);
    //- ML: kappa isotropic value for bubble dynamics            
    const auto& bubbleKappa = cloud::lookupIOField<scalar>("bubbleKappa", obr);
    //- ML: Averaging distance for bubble dynamics
    const auto& averagingDistance = cloud::lookupIOField<scalar>("averagingDistance", obr);
    //- ML: ODE solver type for bubble dynamics
    //const auto& ODESolverType = cloud::lookupIOField<word>("ODESolverType", obr);
    //- ML: flag for activation of EulerianToLagrangian
    const auto& EulerianToLagrangian_activation = cloud::lookupIOField<label>("EulerianToLagrangian_activation", obr);
    //- ML: flag for activation of LagrangianToEulerian
    const auto& LagrangianToEulerian_activation = cloud::lookupIOField<label>("LagrangianToEulerian_activation", obr);
    //- ML:sigma = bubbleRadius / deviation
    const auto& LE_deviation = cloud::lookupIOField<scalar>("LE_deviation", obr);       
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    const auto& LE_minCellOccupancy = cloud::lookupIOField<scalar>("LE_minCellOccupancy", obr);          
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    const auto& LE_cellThreshold = cloud::lookupIOField<scalar>("LE_cellThreshold", obr);        
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    const auto& LE_alphaThreshold = cloud::lookupIOField<scalar>("LE_alphaThreshold", obr);      
    //- ML: Minimum bubble radius threshold that will be tracked
    const auto& LE_bubbleSizeThreshold = cloud::lookupIOField<scalar>("LE_bubbleSizeThreshold", obr);        
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    const auto& LE_boxCheckEnabled = cloud::lookupIOField<label>("LE_boxCheckEnabled", obr);
    //- ML: the top-left corner coordinates
    const auto& LE_boxTopLeftCorner = cloud::lookupIOField<vector>("LE_boxTopLeftCorner", obr); 
    //- ML: the bottom-right corner coordinates
    const auto& LE_boxBottomRightCorner = cloud::lookupIOField<vector>("LE_boxBottomRightCorner", obr);


    label i = 0;

    for (KinematicBubbleParcel<ParcelType>& p : c)
    {
        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        //- ML: Growth rate of radius [m/s]
        p.R_dot_ = R_dot[i];
        p.dTarget_ = dTarget[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];
        p.UCorrect_ = UCorrect[i];
        //- ML: flag for activation of bubble dynamics
        p.bubble_activation_ = bubble_activation[i];
        //- ML: surface tension for bubble dynamics
        p.bubbleSigma_ = bubbleSigma[i];
        //- ML: initial pressure for bubble dynamics
        p.p0_ = p0[i];
        //- ML: initial radius for bubble dynamics            
        p.R0_ = R0[i];
        //- ML: vapour pressure for bubble dynamics            
        p.pv_ = pv[i];
        //- ML: time step for each loop of RP solver          
        p.RPdT_ = RPdT[i];
        //- ML: kappa isotropic value for bubble dynamics            
        p.bubbleKappa_ = bubbleKappa[i];
        //- ML: Averaging distance for bubble dynamics
        p.averagingDistance_ = averagingDistance[i];
        //- ML: ODE solver type for bubble dynamics
        //p.ODESolverType_ = ODESolverType[i];
        //- ML: flag for activation of bubble dynamics
        p.EulerianToLagrangian_activation_ = EulerianToLagrangian_activation[i];
        //- ML: flag for activation of bubble dynamics
        p.LagrangianToEulerian_activation_ = LagrangianToEulerian_activation[i];
        //- ML:sigma = bubbleRadius / deviation
        p.LE_deviation_ = LE_deviation[i];     
        //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
        p.LE_minCellOccupancy_ = LE_minCellOccupancy[i];        
        //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
        p.LE_cellThreshold_ = LE_cellThreshold[i];         
        //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
        p.LE_alphaThreshold_ = LE_alphaThreshold[i];        
        //- ML: Minimum bubble radius threshold that will be tracked
        p.LE_bubbleSizeThreshold_ = LE_bubbleSizeThreshold[i];         
        //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
        p.LE_boxCheckEnabled_ = LE_boxCheckEnabled[i];
        //- ML: the top-left corner coordinates
        p.LE_boxTopLeftCorner_ = LE_boxTopLeftCorner[i]; 
        //- ML: the bottom-right corner coordinates
        p.LE_boxBottomRightCorner_ = LE_boxBottomRightCorner[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicBubbleParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    const label np = c.size();

    auto& active = cloud::createIOField<label>("active", np, obr);
    auto& typeId = cloud::createIOField<label>("typeId", np, obr);
    auto& nParticle = cloud::createIOField<scalar>("nParticle", np, obr);
    auto& d = cloud::createIOField<scalar>("d", np, obr);
    //- ML: Growth rate of radius [m/s]
    auto& R_dot = cloud::createIOField<scalar>("R_dot", np, obr);
    auto& dTarget = cloud::createIOField<scalar>("dTarget", np, obr);
    auto& U = cloud::createIOField<vector>("U", np, obr);
    auto& rho = cloud::createIOField<scalar>("rho", np, obr);
    auto& age = cloud::createIOField<scalar>("age", np, obr);
    auto& tTurb = cloud::createIOField<scalar>("tTurb", np, obr);
    auto&& UTurb = cloud::createIOField<vector>("UTurb", np, obr);
    auto&& UCorrect = cloud::createIOField<vector>("UCorrect", np, obr);
    //- ML: flag for activation of bubble dynamics
    auto& bubble_activation = cloud::createIOField<label>("bubble_activation", np, obr);
    //- ML: surface tension for bubble dynamics
    auto& bubbleSigma = cloud::createIOField<scalar>("bubbleSigma", np, obr);
    //- ML: initial pressure for bubble dynamics
    auto& p0 = cloud::createIOField<scalar>("p0", np, obr);
    //- ML: initial radius for bubble dynamics            
    auto& R0 = cloud::createIOField<scalar>("R0", np, obr);
    //- ML: vapour pressure for bubble dynamics            
    auto& pv = cloud::createIOField<scalar>("pv", np, obr);
    //- ML: time step for each loop of RP solver           
    auto& RPdT = cloud::createIOField<scalar>("RPdT", np, obr);
    //- ML: kappa isotropic value for bubble dynamics            
    auto& bubbleKappa = cloud::createIOField<scalar>("bubbleKappa", np, obr);
    //- ML: Averaging distance for bubble dynamics
    auto& averagingDistance = cloud::createIOField<scalar>("averagingDistance", np, obr);
    //- ML: ODE solver type for bubble dynamics
    //auto& ODESolverType = cloud::createIOField<word>("ODESolverType", np, obr); 
    //- ML: flag for activation of EulerianToLagrangian
    auto& EulerianToLagrangian_activation = cloud::createIOField<label>("EulerianToLagrangian_activation", np, obr);
    //- ML: flag for activation of LagrangianToEulerian
    auto& LagrangianToEulerian_activation = cloud::createIOField<label>("LagrangianToEulerian_activation", np, obr);
    //- ML:sigma = bubbleRadius / deviation
    auto& LE_deviation = cloud::createIOField<scalar>("LE_deviation", np, obr);        
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    auto& LE_minCellOccupancy = cloud::createIOField<scalar>("LE_minCellOccupancy", np, obr);          
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    auto& LE_cellThreshold = cloud::createIOField<scalar>("LE_cellThreshold", np, obr);         
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    auto& LE_alphaThreshold = cloud::createIOField<scalar>("LE_alphaThreshold", np, obr);       
    //- ML: Minimum bubble radius threshold that will be tracked
    auto& LE_bubbleSizeThreshold = cloud::createIOField<scalar>("LE_bubbleSizeThreshold", np, obr);       
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    auto& LE_boxCheckEnabled = cloud::createIOField<label>("LE_boxCheckEnabled", np, obr);
    //- ML: the top-left corner coordinates
    auto& LE_boxTopLeftCorner = cloud::createIOField<vector>("LE_boxTopLeftCorner", np, obr); 
    //- ML: the bottom-right corner coordinates
    auto& LE_boxBottomRightCorner = cloud::createIOField<vector>("LE_boxBottomRightCorner", np, obr);

    label i = 0;

    for (const KinematicBubbleParcel<ParcelType>& p : c)
    {
        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        //- ML: Growth rate of radius [m/s]
        R_dot[i] = p.R_dot();
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();
        UCorrect[i] = p.UCorrect();
        //- ML: flag for activation of bubble dynamics
        bubble_activation[i] = p.bubble_activation();
        //- ML: surface tension for bubble dynamics
        bubbleSigma[i] = p.bubbleSigma();
        //- ML: initial pressure for bubble dynamics
        p0[i] = p.p0();
        //- ML: initial radius for bubble dynamics            
        R0[i] = p.R0();
        //- ML: vapour pressure for bubble dynamics            
        pv[i] = p.pv();
        //- ML: time step for each loop of RP solver           
        RPdT[i] = p.RPdT();
        //- ML: kappa isotropic value for bubble dynamics            
        bubbleKappa[i] = p.bubbleKappa();
        //- ML: Averaging distance for bubble dynamics
        averagingDistance[i] = p.averagingDistance();
        //- ML: ODE solver type for bubble dynamics
        //ODESolverType[i] = p.ODESolverType();
        //- ML: flag for activation of EulerianToLagrangian
        EulerianToLagrangian_activation[i] = p.EulerianToLagrangian_activation();
        //- ML: flag for activation of LagrangianToEulerian
        LagrangianToEulerian_activation[i] = p.LagrangianToEulerian_activation();
        //- ML:sigma = bubbleRadius / deviation
        LE_deviation[i] = p.LE_deviation();       
        //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
        LE_minCellOccupancy[i] = p.LE_minCellOccupancy();         
        //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
        LE_cellThreshold[i] = p.LE_cellThreshold();         
        //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
        LE_alphaThreshold[i] = p.LE_alphaThreshold();       
        //- ML: Minimum bubble radius threshold that will be tracked
        LE_bubbleSizeThreshold[i] = p.LE_bubbleSizeThreshold();       
        //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
        LE_boxCheckEnabled[i] = p.LE_boxCheckEnabled();
        //- ML: the top-left corner coordinates
        LE_boxTopLeftCorner[i] = p.LE_boxTopLeftCorner(); 
        //- ML: the bottom-right corner coordinates
        LE_boxBottomRightCorner[i] = p.LE_boxBottomRightCorner();

        ++i;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const KinematicBubbleParcel<ParcelType>& p
)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << bool(p.active())
            << token::SPACE << p.typeId()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            //- ML: Growth rate of radius [m/s]
            << token::SPACE << p.R_dot()
            << token::SPACE << p.dTarget()
            << token::SPACE << p.U()
            << token::SPACE << p.rho()
            << token::SPACE << p.age()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb()
            << token::SPACE << p.UCorrect()
            //- ML: flag for activation of bubble dynamics
            << token::SPACE << bool(p.bubble_activation())
            //- ML: surface tension for bubble dynamics
            << token::SPACE << p.bubbleSigma()
            //- ML: initial pressure for bubble dynamics
            << token::SPACE << p.p0()
            //- ML: initial radius for bubble dynamics            
            << token::SPACE << p.R0()
            //- ML: vapour pressure for bubble dynamics            
            << token::SPACE << p.pv()
            //- ML: time step for each loop of RP solver          
            << token::SPACE << p.RPdT()
            //- ML: kappa isotropic value for bubble dynamics            
            << token::SPACE << p.bubbleKappa()
            //- ML: Averaging distance for bubble dynamics
            << token::SPACE << p.averagingDistance()
            //- ML: ODE solver type for bubble dynamics
            //<< token::SPACE << p.ODESolverType();
            //- ML: flag for activation of EulerianToLagrangian
            << token::SPACE << bool(p.EulerianToLagrangian_activation())
            //- ML: flag for activation of LagrangianToEulerian
            << token::SPACE << bool(p.LagrangianToEulerian_activation())
            //- ML:sigma = bubbleRadius / deviation
            << token::SPACE << p.LE_deviation()  
            //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
            << token::SPACE << p.LE_minCellOccupancy()         
            //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
            << token::SPACE << p.LE_cellThreshold()  
            //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
            << token::SPACE << p.LE_alphaThreshold()       
            //- ML: Minimum bubble radius threshold that will be tracked
            << token::SPACE << p.LE_bubbleSizeThreshold()        
            //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
            << token::SPACE << bool(p.LE_boxCheckEnabled())
            //- ML: the top-left corner coordinates
            << token::SPACE << p.LE_boxTopLeftCorner() 
            //- ML: the bottom-right corner coordinates
            << token::SPACE << p.LE_boxBottomRightCorner();

    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            KinematicBubbleParcel<ParcelType>::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
