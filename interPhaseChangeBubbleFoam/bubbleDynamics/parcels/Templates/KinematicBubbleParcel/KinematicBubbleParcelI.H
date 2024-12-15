/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline
Foam::KinematicBubbleParcel<ParcelType>::constantProperties::constantProperties()
:
    dict_(),
    parcelTypeId_(dict_, -1),
    rhoMin_(dict_, 0.0),
    rho0_(dict_, 0.0),
    minParcelMass_(dict_, 0.0)
{}


template<class ParcelType>
inline Foam::KinematicBubbleParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    dict_(cp.dict_),
    parcelTypeId_(cp.parcelTypeId_),
    rhoMin_(cp.rhoMin_),
    rho0_(cp.rho0_),
    minParcelMass_(cp.minParcelMass_)
{}


template<class ParcelType>
inline Foam::KinematicBubbleParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    dict_(parentDict.subOrEmptyDict("constantProperties")),
    parcelTypeId_(dict_, "parcelTypeId", -1),
    rhoMin_(dict_, "rhoMin", 1e-15),
    rho0_(dict_, "rho0", -1),
    minParcelMass_(dict_, "minParcelMass", 1e-15)
{}


template<class ParcelType>
inline Foam::KinematicBubbleParcel<ParcelType>::KinematicBubbleParcel
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(owner, coordinates, celli, tetFacei, tetPti),
    active_(true),
    typeId_(-1),
    nParticle_(0),
    d_(0.0),
    //- ML: Growth rate of radius [m/s]
    R_dot_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    //- ML: Read cloud parameters from dictionary
    cloudProperties_(owner.lookupObject<IOdictionary>("kinematicBubbleCloudProperties")),
    //- ML: Read bubble parameters from dictionary    
    bubbleProperties_(cloudProperties_.subDict("bubbleProperties")),
    //- ML: flag for activation of bubble dynamics
    bubble_activation_(bubbleProperties_.lookupOrDefault<bool>("bubble_activation", 0)),
    //- ML: surface tension for bubble dynamics
    bubbleSigma_(bubbleProperties_.lookupOrDefault<scalar>("bubbleSigma", 0.0)),
    //- ML: initial pressure for bubble dynamics
    p0_(bubbleProperties_.lookupOrDefault<scalar>("p0", 0.0)),
    //- ML: initial radius for bubble dynamics            
    R0_(bubbleProperties_.lookupOrDefault<scalar>("R0", 0.0)),
    //- ML: vapour pressure for bubble dynamics            
    pv_(bubbleProperties_.lookupOrDefault<scalar>("pv", 0.0)), 
    //- ML: time step for each loop of RP solver            
    RPdT_(bubbleProperties_.lookupOrDefault<scalar>("RPdT", 0.0)),
    //- ML: kappa isotropic value for bubble dynamics            
    bubbleKappa_(bubbleProperties_.lookupOrDefault<scalar>("bubbleKappa", 0.0)),
    //- ML: Averaging distance for bubble dynamics
    averagingDistance_(bubbleProperties_.lookupOrDefault<scalar>("averagingDistance", 0.0)),
    //- ML: ODE solver type for bubble dynamics
    ODESolverType_(bubbleProperties_.lookupOrDefault<word>("ODESolverType", "defaultSolver")),
    //- ML: Read Transition_Algorithm from dictionary    
    Transition_Algorithm_(cloudProperties_.subDict("Transition_Algorithm")),
    //- ML: flag for activation of EulerianToLagrangian
    EulerianToLagrangian_activation_(Transition_Algorithm_.lookupOrDefault<bool>("EulerianToLagrangian_activation", 0)),
    //- ML: flag for activation of LagrangianToEulerian
    LagrangianToEulerian_activation_(Transition_Algorithm_.lookupOrDefault<bool>("LagrangianToEulerian_activation", 0)),
    //- ML:sigma = bubbleRadius / deviation
    LE_deviation_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_deviation", 0.0)),        
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    LE_minCellOccupancy_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_minCellOccupancy", 0.0)),          
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    LE_cellThreshold_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_cellThreshold", 0.0)),          
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    LE_alphaThreshold_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_alphaThreshold", 0.0)),         
    //- ML: Minimum bubble radius threshold that will be tracked
    LE_bubbleSizeThreshold_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_bubbleSizeThreshold", 0.0)),         
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    LE_boxCheckEnabled_(Transition_Algorithm_.lookupOrDefault<bool>("LE_boxCheckEnabled", 0)),
    //- ML: the top-left corner coordinates
    LE_boxTopLeftCorner_(Transition_Algorithm_.lookupOrDefault<vector>("LE_boxTopLeftCorner", Zero)),  
    //- ML: the bottom-right corner coordinates
    LE_boxBottomRightCorner_(Transition_Algorithm_.lookupOrDefault<vector>("LE_boxBottomRightCorner", Zero))
{}


template<class ParcelType>
inline Foam::KinematicBubbleParcel<ParcelType>::KinematicBubbleParcel
(
    const polyMesh& owner,
    const vector& position,
    const label celli
)
:
    ParcelType(owner, position, celli),
    active_(true),
    typeId_(-1),
    nParticle_(0),
    d_(0.0),
    //- Growth rate of radius [m/s]
    R_dot_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    //- ML: Read cloud parameters from dictionary
    cloudProperties_(owner.lookupObject<IOdictionary>("kinematicBubbleCloudProperties")),
    //- ML: Read bubble parameters from dictionary    
    bubbleProperties_(cloudProperties_.subDict("bubbleProperties")),
    //- ML: flag for activation of bubble dynamics
    bubble_activation_(bubbleProperties_.lookupOrDefault<bool>("bubble_activation", 0)),
    //- ML: surface tension for bubble dynamics
    bubbleSigma_(bubbleProperties_.lookupOrDefault<scalar>("bubbleSigma", 0.0)),
    //- ML: initial pressure for bubble dynamics
    p0_(bubbleProperties_.lookupOrDefault<scalar>("p0", 0.0)),
    //- ML: initial radius for bubble dynamics            
    R0_(bubbleProperties_.lookupOrDefault<scalar>("R0", 0.0)),
    //- ML: vapour pressure for bubble dynamics            
    pv_(bubbleProperties_.lookupOrDefault<scalar>("pv", 0.0)), 
    //- ML: time step for each loop of RP solver            
    RPdT_(bubbleProperties_.lookupOrDefault<scalar>("RPdT", 0.0)),
    //- ML: kappa isotropic value for bubble dynamics            
    bubbleKappa_(bubbleProperties_.lookupOrDefault<scalar>("bubbleKappa", 0.0)),
    //- ML: Averaging distance for bubble dynamics
    averagingDistance_(bubbleProperties_.lookupOrDefault<scalar>("averagingDistance", 0.0)),
    //- ML: ODE solver type for bubble dynamics
    ODESolverType_(bubbleProperties_.lookupOrDefault<word>("ODESolverType", "defaultSolver")),
    //- ML: Read Transition_Algorithm from dictionary    
    Transition_Algorithm_(cloudProperties_.subDict("Transition_Algorithm")),
    //- ML: flag for activation of EulerianToLagrangian
    EulerianToLagrangian_activation_(Transition_Algorithm_.lookupOrDefault<bool>("EulerianToLagrangian_activation", 0)),
    //- ML: flag for activation of LagrangianToEulerian
    LagrangianToEulerian_activation_(Transition_Algorithm_.lookupOrDefault<bool>("LagrangianToEulerian_activation", 0)),
    //- ML:sigma = bubbleRadius / deviation
    LE_deviation_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_deviation", 0.0)),        
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    LE_minCellOccupancy_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_minCellOccupancy", 0.0)),          
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    LE_cellThreshold_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_cellThreshold", 0.0)),          
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    LE_alphaThreshold_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_alphaThreshold", 0.0)),         
    //- ML: Minimum bubble radius threshold that will be tracked
    LE_bubbleSizeThreshold_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_bubbleSizeThreshold", 0.0)),         
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    LE_boxCheckEnabled_(Transition_Algorithm_.lookupOrDefault<bool>("LE_boxCheckEnabled", 0)),
    //- ML: the top-left corner coordinates
    LE_boxTopLeftCorner_(Transition_Algorithm_.lookupOrDefault<vector>("LE_boxTopLeftCorner", Zero)),  
    //- ML: the bottom-right corner coordinates
    LE_boxBottomRightCorner_(Transition_Algorithm_.lookupOrDefault<vector>("LE_boxBottomRightCorner", Zero))
{
}


template<class ParcelType>
inline Foam::KinematicBubbleParcel<ParcelType>::KinematicBubbleParcel
(
    const polyMesh& owner,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const constantProperties& constProps
)
:
    ParcelType(owner, coordinates, celli, tetFacei, tetPti),
    active_(true),
    typeId_(typeId),
    nParticle_(nParticle0),
    d_(d0),
    //- ML:  Growth rate of radius [m/s]
    R_dot_(0.0),
    dTarget_(dTarget0),
    U_(U0),
    rho_(constProps.rho0()),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    //- ML: Read cloud parameters from dictionary
    cloudProperties_(owner.lookupObject<IOdictionary>("kinematicBubbleCloudProperties")),
    //- ML: Read bubble parameters from dictionary    
    bubbleProperties_(cloudProperties_.subDict("bubbleProperties")),
    //- ML: flag for activation of bubble dynamics
    bubble_activation_(bubbleProperties_.lookupOrDefault<bool>("bubble_activation", 0)),
    //- ML: surface tension for bubble dynamics
    bubbleSigma_(bubbleProperties_.lookupOrDefault<scalar>("bubbleSigma", 0.0)),
    //- ML: initial pressure for bubble dynamics
    p0_(bubbleProperties_.lookupOrDefault<scalar>("p0", 0.0)),
    //- ML: initial radius for bubble dynamics            
    R0_(bubbleProperties_.lookupOrDefault<scalar>("R0", 0.0)),
    //- ML: vapour pressure for bubble dynamics            
    pv_(bubbleProperties_.lookupOrDefault<scalar>("pv", 0.0)), 
    //- ML: time step for each loop of RP solver            
    RPdT_(bubbleProperties_.lookupOrDefault<scalar>("RPdT", 0.0)),
    //- ML: kappa isotropic value for bubble dynamics            
    bubbleKappa_(bubbleProperties_.lookupOrDefault<scalar>("bubbleKappa", 0.0)),
    //- ML: Averaging distance for bubble dynamics
    averagingDistance_(bubbleProperties_.lookupOrDefault<scalar>("averagingDistance", 0.0)),
    //- ML: ODE solver type for bubble dynamics
    ODESolverType_(bubbleProperties_.lookupOrDefault<word>("ODESolverType", "defaultSolver")),
    //- ML: Read Transition_Algorithm from dictionary    
    Transition_Algorithm_(cloudProperties_.subDict("Transition_Algorithm")),
    //- ML: flag for activation of EulerianToLagrangian
    EulerianToLagrangian_activation_(Transition_Algorithm_.lookupOrDefault<bool>("EulerianToLagrangian_activation", 0)),
    //- ML: flag for activation of LagrangianToEulerian
    LagrangianToEulerian_activation_(Transition_Algorithm_.lookupOrDefault<bool>("LagrangianToEulerian_activation", 0)),
    //- ML:sigma = bubbleRadius / deviation
    LE_deviation_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_deviation", 0.0)),        
    //- ML:Minimum liquid volume fraction value that lagrangian cell can occupy
    LE_minCellOccupancy_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_minCellOccupancy", 0.0)),          
    //- ML:Number of cells that trigger the transtion from Lagrangian to Eulerian
    LE_cellThreshold_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_cellThreshold", 0.0)),          
    //- ML: Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
    LE_alphaThreshold_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_alphaThreshold", 0.0)),         
    //- ML: Minimum bubble radius threshold that will be tracked
    LE_bubbleSizeThreshold_(Transition_Algorithm_.lookupOrDefault<scalar>("LE_bubbleSizeThreshold", 0.0)),         
    //- ML: Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
    LE_boxCheckEnabled_(Transition_Algorithm_.lookupOrDefault<bool>("LE_boxCheckEnabled", 0)),
    //- ML: the top-left corner coordinates
    LE_boxTopLeftCorner_(Transition_Algorithm_.lookupOrDefault<vector>("LE_boxTopLeftCorner", Zero)),  
    //- ML: the bottom-right corner coordinates
    LE_boxBottomRightCorner_(Transition_Algorithm_.lookupOrDefault<vector>("LE_boxBottomRightCorner", Zero))

{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicBubbleParcel<ParcelType>::constantProperties::dict() const
{
    return dict_;
}


template<class ParcelType>
inline Foam::label
Foam::KinematicBubbleParcel<ParcelType>::constantProperties::parcelTypeId() const
{
    return parcelTypeId_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicBubbleParcel<ParcelType>::constantProperties::rhoMin() const
{
    return rhoMin_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicBubbleParcel<ParcelType>::constantProperties::rho0() const
{
    return rho0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::KinematicBubbleParcel<ParcelType>::constantProperties::minParcelMass() const
{
    return minParcelMass_.value();
}


// * * * * * * * KinematicBubbleParcel Member Functions  * * * * * * * //

template<class ParcelType>
inline bool Foam::KinematicBubbleParcel<ParcelType>::active() const
{
    return active_;
}


template<class ParcelType>
inline Foam::label Foam::KinematicBubbleParcel<ParcelType>::typeId() const
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::nParticle() const
{
    return nParticle_;
}

//- ML: Return const access to Growth rate of radius [m/s]
template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::R_dot() const
{
    return R_dot_;
}

template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::d() const
{
    return d_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::dTarget() const
{
    return dTarget_;
}


template<class ParcelType>
inline const Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::U() const
{
    return U_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::rho() const
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::age() const
{
    return age_;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::tTurb() const
{
    return tTurb_;
}


template<class ParcelType>
inline const Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::UTurb() const
{
    return UTurb_;
}

//- ML: Return bubble cloud dictionary
template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicBubbleParcel<ParcelType>::cloudProperties() const
{
    return cloudProperties_;
}

//- ML: Return bubble bubble dictionary
template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicBubbleParcel<ParcelType>::bubbleProperties() const
{
    return bubbleProperties_;
}

//- ML: Return const access to flag for activation of bubble dynamics
template<class ParcelType>
inline bool Foam::KinematicBubbleParcel<ParcelType>::bubble_activation() const
{
    return bubble_activation_;
}

//- ML: Return const access to surface tension for bubble dynamics
template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::bubbleSigma() const
{
    return bubbleSigma_;
}
            
//- ML: Return const access to initial pressure for bubble dynamics
template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::p0() const
{
    return p0_;
}

//- ML: Return const access to initial radius for bubble dynamics
template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::R0() const
{
    return R0_;
}

//- ML: Return const access to vapour pressure for bubble dynamics  
template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::pv() const
{
    return pv_;
}          


//- ML: Return const access to time step for each loop of RP solver          
template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::RPdT() const
{
    return RPdT_;
}          

//- ML: Return const access to kappa isotropic value for bubble dynamics 
template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::bubbleKappa() const
{
    return bubbleKappa_;
}         

//- ML: Return const access to Averaging distance for bubble dynamics
template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::averagingDistance() const
{
    return averagingDistance_;
}

//- ML: Return const access to ODE solver type for bubble dynamics
template<class ParcelType>
inline const Foam::word Foam::KinematicBubbleParcel<ParcelType>::ODESolverType() const
{
    return ODESolverType_;
}       


//- ML: Return Transition_Algorithm from dictionary   
template<class ParcelType>
inline const Foam::dictionary&
Foam::KinematicBubbleParcel<ParcelType>::Transition_Algorithm() const
{
    return Transition_Algorithm_;
}

//- ML: Return const access to flag for activation of EulerianToLagrangian
template<class ParcelType>
inline bool Foam::KinematicBubbleParcel<ParcelType>::EulerianToLagrangian_activation() const
{
    return EulerianToLagrangian_activation_;
}

//- ML: Return const access to flag for activation of LagrangianToEulerian
template<class ParcelType>
inline bool Foam::KinematicBubbleParcel<ParcelType>::LagrangianToEulerian_activation() const
{
    return LagrangianToEulerian_activation_;
}

//- ML: Return const access to sigma = bubbleRadius / deviation
template<class ParcelType>        
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::LE_deviation() const
{
    return LE_deviation_;
}       

//- ML: Return const access to Minimum liquid volume fraction value that lagrangian cell can occupy  
template<class ParcelType>      
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::LE_minCellOccupancy() const
{
    return LE_minCellOccupancy_;
}       

//- ML: Return const access to Number of cells that trigger the transtion from Lagrangian to Eulerian  
template<class ParcelType>      
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::LE_cellThreshold() const
{
    return LE_cellThreshold_;
}       

//- ML: Return const access to Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
template<class ParcelType>        
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::LE_alphaThreshold() const
{
    return LE_alphaThreshold_;
}       

//- ML: Return const access to Minimum bubble radius threshold that will be tracked       
template<class ParcelType> 
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::LE_bubbleSizeThreshold() const
{
    return LE_bubbleSizeThreshold_;
}       

//- ML: Return const access to Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
template<class ParcelType>
inline bool Foam::KinematicBubbleParcel<ParcelType>::LE_boxCheckEnabled() const
{
    return LE_boxCheckEnabled_;
}       

//- ML: Return const access to the top-left corner coordinates
template<class ParcelType>
inline const Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::LE_boxTopLeftCorner() const
{
    return LE_boxTopLeftCorner_;
}       

//- ML: Return const access to the bottom-right corner coordinates
template<class ParcelType>
inline const Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::LE_boxBottomRightCorner() const
{
    return LE_boxBottomRightCorner_;
}        

template<class ParcelType>
inline void Foam::KinematicBubbleParcel<ParcelType>::active(const bool state)
{
    active_ = state;
}


template<class ParcelType>
inline Foam::label& Foam::KinematicBubbleParcel<ParcelType>::typeId()
{
    return typeId_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::nParticle()
{
    return nParticle_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::d()
{
    return d_;
}

//- ML: Return access to Growth rate of radius [m/s]
template<class ParcelType>
inline Foam::scalar&  Foam::KinematicBubbleParcel<ParcelType>::R_dot()
{
    return R_dot_;
}

template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::dTarget()
{
    return dTarget_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::U()
{
    return U_;
}


template<class ParcelType>
inline const Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::UCorrect() const
{
    return UCorrect_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::UCorrect()
{
    return UCorrect_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::rho()
{
    return rho_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::age()
{
    return age_;
}


template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::tTurb()
{
    return tTurb_;
}


template<class ParcelType>
inline Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::UTurb()
{
    return UTurb_;
}

//- ML: Return const access to flag for activation of bubble dynamics
template<class ParcelType>
inline void Foam::KinematicBubbleParcel<ParcelType>::bubble_activation(const bool state)
{
    bubble_activation_ = state;
}

//- ML: Return const access to surface tension for bubble dynamics
template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::bubbleSigma()
{
    return bubbleSigma_;
}
            
//- ML: Return const access to initial pressure for bubble dynamics
template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::p0()
{
    return p0_;
}

//- ML: Return const access to initial radius for bubble dynamics
template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::R0()
{
    return R0_;
}

//- ML: Return const access to vapour pressure for bubble dynamics  
template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::pv()
{
    return pv_;
}          

//- ML: Return const access to time step for each loop of RP solver          
template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::RPdT()
{
    return RPdT_;
} 

//- ML: Return const access to kappa isotropic value for bubble dynamics 
template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::bubbleKappa()
{
    return bubbleKappa_;
}                

//- ML: Return const access to Averaging distance for bubble dynamics
template<class ParcelType>
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::averagingDistance()
{
    return averagingDistance_;
}

//- ML: Return const access to ODE solver type for bubble dynamics
template<class ParcelType>
inline void Foam::KinematicBubbleParcel<ParcelType>::ODESolverType(const word ODEtype)
{
    ODESolverType_ = ODEtype;
} 

//- ML: change flag for activation of EulerianToLagrangian
template<class ParcelType>
inline void Foam::KinematicBubbleParcel<ParcelType>::EulerianToLagrangian_activation(const bool state)
{
    EulerianToLagrangian_activation_ = state;
}

//- ML: change flag for activation of LagrangianToEulerian
template<class ParcelType>
inline void Foam::KinematicBubbleParcel<ParcelType>::LagrangianToEulerian_activation(const bool state)
{
    LagrangianToEulerian_activation_ = state;
}

//- ML: Return  access to sigma = bubbleRadius / deviation
template<class ParcelType>        
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::LE_deviation() 
{
    return LE_deviation_;
}       

//- ML: Return  access to Minimum liquid volume fraction value that lagrangian cell can occupy  
template<class ParcelType>      
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::LE_minCellOccupancy()
{
    return LE_minCellOccupancy_;
}       

//- ML: Return  access to Number of cells that trigger the transtion from Lagrangian to Eulerian  
template<class ParcelType>      
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::LE_cellThreshold() 
{
    return LE_cellThreshold_;
}       

//- ML: Return  access to Set threshold for interface proximity that trigger the transtion from Lagrangian to Eulerian
template<class ParcelType>        
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::LE_alphaThreshold() 
{
    return LE_alphaThreshold_;
}       

//- ML: Return  access to Minimum bubble radius threshold that will be tracked       
template<class ParcelType> 
inline Foam::scalar& Foam::KinematicBubbleParcel<ParcelType>::LE_bubbleSizeThreshold() 
{
    return LE_bubbleSizeThreshold_;
}       

//- ML: Return  access to Define the box boundaries for tracking- Only lagrangian inside this box will be tracked
template<class ParcelType>
inline void Foam::KinematicBubbleParcel<ParcelType>::LE_boxCheckEnabled(const bool state)
{
    return LE_boxCheckEnabled = state;
}       

//- ML: Return access to the top-left corner coordinates
template<class ParcelType>
inline Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::LE_boxTopLeftCorner() 
{
    return LE_boxTopLeftCorner_;
}       

//- ML: Return  access to the bottom-right corner coordinates
template<class ParcelType>
inline Foam::vector& Foam::KinematicBubbleParcel<ParcelType>::LE_boxBottomRightCorner() 
{
    return LE_boxBottomRightCorner_;
}        


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::massCell
(
    const trackingData& td
) const
{
    return td.rhoc()*this->mesh().cellVolumes()[this->cell()];
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::mass() const
{
    return rho_*volume();
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::momentOfInertia() const
{
    return 0.1*mass()*sqr(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::volume() const
{
    return volume(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::volume(const scalar d)
{
    return pi/6.0*pow3(d);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::areaP() const
{
    return areaP(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::areaP(const scalar d)
{
    return 0.25*areaS(d);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::areaS() const
{
    return areaS(d_);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::areaS(const scalar d)
{
    return pi*d*d;
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::Re
(
    const trackingData& td
) const
{
    return Re(td.rhoc(), U_, td.Uc(), d_, td.muc());
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::Re
(
    const scalar rhoc,
    const vector& U,
    const vector& Uc,
    const scalar d,
    const scalar muc
)
{
    return rhoc*mag(U - Uc)*d/max(muc, ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::We
(
    const trackingData& td,
    const scalar sigma
) const
{
    return We(td.rhoc(), U_, td.Uc(), d_, sigma);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::We
(
    const scalar rhoc,
    const vector& U,
    const vector& Uc,
    const scalar d,
    const scalar sigma
)
{
    return rhoc*magSqr(U - Uc)*d/max(sigma, ROOTVSMALL);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::Eo
(
    const trackingData& td,
    const scalar sigma
) const
{
    return Eo(td.g(), rho_, td.rhoc(), U_, d_, sigma);
}


template<class ParcelType>
inline Foam::scalar Foam::KinematicBubbleParcel<ParcelType>::Eo
(
    const vector& g,
    const scalar rho,
    const scalar rhoc,
    const vector& U,
    const scalar d,
    const scalar sigma
)
{
    const vector dir = U/max(mag(U), ROOTVSMALL);
    return mag(g & dir)*mag(rho - rhoc)*sqr(d)/max(sigma, ROOTVSMALL);
}


// ************************************************************************* //
