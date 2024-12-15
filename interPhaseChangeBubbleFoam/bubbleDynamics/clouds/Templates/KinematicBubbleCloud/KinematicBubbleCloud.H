/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::KinematicBubbleCloud

Group
    grpLagrangianIntermediateClouds

Description
    Templated base class for kinematicBubble cloud

    - cloud function objects

    - particle forces, e.g.
      - buoyancy
      - drag
      - pressure gradient
      - ...

    - sub-models:
      - dispersion model
      - injection model
      - patch interaction model
      - stochastic collision model
      - surface film model

SourceFiles
    KinematicBubbleCloudI.H
    KinematicBubbleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicBubbleCloud_H
#define KinematicBubbleCloud_H

#include "particle.H"
#include "Cloud.H"
#include "kinematicBubbleCloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "cloudSolution.H"

#include "ParticleForceList.H"
#include "CloudFunctionObjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class integrationScheme;

template<class CloudType>
class InjectionModelList;

template<class CloudType>
class DispersionModel;

template<class CloudType>
class PatchInteractionModel;

template<class CloudType>
class SurfaceFilmModel;

template<class CloudType>
class StochasticCollisionModel;

template<class CloudType>
class PackingModel;

template<class CloudType>
class DampingModel;

template<class CloudType>
class IsotropyModel;


/*---------------------------------------------------------------------------*\
                       Class KinematicBubbleCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class KinematicBubbleCloud
:
    public CloudType,
    public kinematicBubbleCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef KinematicBubbleCloud<CloudType> kinematicCloudType;

        //- Force models type
        typedef ParticleForceList<KinematicBubbleCloud<CloudType>> forceType;

        //- Function object type
        typedef CloudFunctionObjectList<KinematicBubbleCloud<CloudType>>
            functionType;


private:

    // Private Data

        //- Cloud copy pointer
        autoPtr<KinematicBubbleCloud<CloudType>> cloudCopyPtr_;


    // Private Member Functions

        //- No copy construct
        KinematicBubbleCloud(const KinematicBubbleCloud&) = delete;

        //- No copy assignment
        void operator=(const KinematicBubbleCloud&) = delete;


protected:

    // Protected Data

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- Dictionary of output properties
        IOdictionary outputProperties_;

        //- Solution properties
        cloudSolution solution_;

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;

        //- Sub-models dictionary
        const dictionary subModelProperties_;

        //- Random number generator - used by some injection routines
        mutable Random rndGen_;

        //- Cell occupancy information for each parcel, (demand driven)
        autoPtr<List<DynamicList<parcelType*>>> cellOccupancyPtr_;

        //- Cell length scale
        scalarField cellLengthScale_;


        // References to the carrier gas fields

            //- Density [kg/m3]
            const volScalarField& rho_;

            //- Velocity [m/s]
            const volVectorField& U_;

            //- Dynamic viscosity [Pa.s]
            const volScalarField& mu_;

            //- ML: Pressure [Pa]
            const volScalarField& p_;

            //- ML: Eulerian Liquid Volume Fraction
            volScalarField& alphal_;

            //- ML: Lagrangian Liquid Volume Fraction
            volScalarField& betal_;

            //- ML: Old Lagrangian Liquid Volume Fraction
            volScalarField betal_old_;


        // Environmental properties

            //- Gravity
            const dimensionedVector& g_;

            //- Averaged ambient domain pressure
            scalar pAmbient_;


        //- Optional particle forces
        forceType forces_;

        //- Optional cloud function objects
        functionType functions_;


        // References to the cloud sub-models

            //- Injector models
            InjectionModelList<KinematicBubbleCloud<CloudType>> injectors_;

            //- Dispersion model
            autoPtr<DispersionModel<KinematicBubbleCloud<CloudType>>>
                dispersionModel_;

            //- Patch interaction model
            autoPtr<PatchInteractionModel<KinematicBubbleCloud<CloudType>>>
                patchInteractionModel_;

            //- Stochastic collision model
            autoPtr<StochasticCollisionModel<KinematicBubbleCloud<CloudType>>>
                stochasticCollisionModel_;

            //- Surface film model
            autoPtr<SurfaceFilmModel<KinematicBubbleCloud<CloudType>>>
                surfaceFilmModel_;

            //- Packing model
            autoPtr<PackingModel<KinematicBubbleCloud<CloudType>>>
                packingModel_;

            //- Damping model
            autoPtr<DampingModel<KinematicBubbleCloud<CloudType>>>
                dampingModel_;

            //- Exchange model
            autoPtr<IsotropyModel<KinematicBubbleCloud<CloudType>>>
                isotropyModel_;


        // Reference to the particle integration schemes

            //- Velocity integration
            autoPtr<integrationScheme> UIntegrator_;


        // Sources

            //- Mass for kinematicBubble cloud
            autoPtr<volScalarField::Internal> rhokTrans_;

            //- Momentum
            autoPtr<volVectorField::Internal> UTrans_;

            //- Coefficient for carrier phase U equation
            autoPtr<volScalarField::Internal> UCoeff_;

            //- ML: Rate of change of mass
            autoPtr<volScalarField::Internal> mDotBubble_;

            //- ML: Source term for Lagrangian Liquid Volume Fraction
            autoPtr<volScalarField::Internal> Dbetal_;

        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Solve the cloud - calls all evolution functions
            template<class TrackCloudType>
            void solve
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );

            //- Build the cellOccupancy
            void buildCellOccupancy();

            //- Update (i.e. build) the cellOccupancy if it has
            //  already been used
            void updateCellOccupancy();

            //- Evolve the cloud
            template<class TrackCloudType>
            void evolveCloud
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );

            //- Post-evolve
            void postEvolve(const typename parcelType::trackingData& td);

            //- Reset state of cloud
            void cloudReset(KinematicBubbleCloud<CloudType>& c);


public:

    // Public Data

        //- Flag to write log into Info
        bool log;


    // Constructors

        //- ML: Construct given carrier gas fields
        KinematicBubbleCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const volScalarField& p,
            volScalarField& alphal, //- ML: Eulerian Liquid Volume Fraction
            volScalarField& betal, //- ML: Lagrangian Liquid Volume Fraction
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        KinematicBubbleCloud
        (
            KinematicBubbleCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        KinematicBubbleCloud
        (
            const fvMesh& mesh,
            const word& name,
            const KinematicBubbleCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new KinematicBubbleCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new KinematicBubbleCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~KinematicBubbleCloud() = default;


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const KinematicBubbleCloud& cloudCopy() const;


            // References to the mesh and databases

                //- Return reference to the mesh
                inline const fvMesh& mesh() const;

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

                //- Return output properties dictionary
                inline const IOdictionary& outputProperties() const;

                //- Return non-const access to the output properties dictionary
                inline IOdictionary& outputProperties();

                //- Return const access to the solution properties
                inline const cloudSolution& solution() const;

                //- Return access to the solution properties
                inline cloudSolution& solution();

                //- Return the constant properties
                inline const typename parcelType::constantProperties&
                    constProps() const;

                //- Return access to the constant properties
                inline typename parcelType::constantProperties& constProps();

                //- Return reference to the sub-models dictionary
                inline const dictionary& subModelProperties() const;


            // Cloud data

                //- Return reference to the random object
                inline Random& rndGen() const;

                //- Return the cell occupancy information for each
                //  parcel, non-const access, the caller is
                //  responsible for updating it for its own purposes
                //  if particles are removed or created.
                inline List<DynamicList<parcelType*>>& cellOccupancy();

                //- Return the cell length scale
                inline const scalarField& cellLengthScale() const;


            // References to the carrier gas fields

                //- Return carrier gas velocity
                inline const volVectorField& U() const;

                //- Return carrier gas density
                inline const volScalarField& rho() const;

                //- Return carrier gas dynamic viscosity
                inline const volScalarField& mu() const;

                //- ML: Return carrier gas pressure
                inline const volScalarField& p() const;

                //- ML: Return Eulerian Liquid Volume Fraction
                inline volScalarField& alphal();

                //- ML: Return Lagrangian Liquid Volume Fraction
                inline volScalarField& betal();

                //- ML: Return  old Lagrangian Liquid Volume Fraction
                inline volScalarField& betal_old();

            // Environmental properties

                //- Gravity
                inline const dimensionedVector& g() const;

                //- Return const-access to the ambient pressure
                inline scalar pAmbient() const;

                //- Return reference to the ambient pressure
                inline scalar& pAmbient();


            //- Optional particle forces
            inline const forceType& forces() const;

            //- Return the optional particle forces
            inline forceType& forces();

            //- Optional cloud function objects
            inline functionType& functions();


            // Sub-models

                //- Return const access to the injection model
                inline const InjectionModelList<KinematicBubbleCloud<CloudType>>&
                    injectors() const;

                //- Return reference to the injection model
                inline InjectionModelList<KinematicBubbleCloud<CloudType>>&
                    injectors();

                //- Return const-access to the dispersion model
                inline const DispersionModel<KinematicBubbleCloud<CloudType>>&
                    dispersion() const;

                //- Return reference to the dispersion model
                inline DispersionModel<KinematicBubbleCloud<CloudType>>&
                    dispersion();

                //- Return const-access to the patch interaction model
                inline const PatchInteractionModel<KinematicBubbleCloud<CloudType>>&
                    patchInteraction() const;

                //- Return reference to the patch interaction model
                inline PatchInteractionModel<KinematicBubbleCloud<CloudType>>&
                    patchInteraction();

                //- Return const-access to the stochastic collision model
                inline const
                    StochasticCollisionModel<KinematicBubbleCloud<CloudType>>&
                    stochasticCollision() const;

                //- Return reference to the stochastic collision model
                inline StochasticCollisionModel<KinematicBubbleCloud<CloudType>>&
                    stochasticCollision();

                //- Return const-access to the surface film model
                inline const SurfaceFilmModel<KinematicBubbleCloud<CloudType>>&
                    surfaceFilm() const;

                //- Return reference to the surface film model
                inline SurfaceFilmModel<KinematicBubbleCloud<CloudType>>&
                    surfaceFilm();


                //- Return const access to the packing model
                inline const PackingModel<KinematicBubbleCloud<CloudType>>&
                    packingModel() const;

                //- Return a reference to the packing model
                inline PackingModel<KinematicBubbleCloud<CloudType>>&
                    packingModel();

                //- Return const access to the damping model
                inline const DampingModel<KinematicBubbleCloud<CloudType>>&
                    dampingModel() const;

                //- Return a reference to the damping model
                inline DampingModel<KinematicBubbleCloud<CloudType>>&
                    dampingModel();

                //- Return const access to the isotropy model
                inline const IsotropyModel<KinematicBubbleCloud<CloudType>>&
                    isotropyModel() const;

                //- Return a reference to the isotropy model
                inline IsotropyModel<KinematicBubbleCloud<CloudType>>&
                    isotropyModel();


            // Integration schemes

                //-Return reference to velocity integration
                inline const integrationScheme& UIntegrator() const;


            // Sources

                //- Transfer the effect of parcel to the carrier phase
                inline void transferToCarrier
                (
                    const parcelType& p,
                    const typename parcelType::trackingData& td
                );


                // Momentum

                    //- Return reference to mass for kinematicBubble source
                    inline volScalarField::Internal& rhokTrans();

                    //- Return const reference to mass for kinematicBubble source
                    inline const volScalarField::Internal& rhokTrans() const;

                    //- Return reference to momentum source
                    inline volVectorField::Internal& UTrans();

                    //- Return const reference to momentum source
                    inline const volVectorField::Internal& UTrans() const;

                    //- Return coefficient for carrier phase U equation
                    inline volScalarField::Internal& UCoeff();

                    //- Return const coefficient for carrier phase U equation
                    inline const volScalarField::Internal& UCoeff() const;

                    //- ML: Return rate of change of mass
                    inline volScalarField::Internal& mDotBubble();

                    //- ML: Return rate of change of mass
                    inline const volScalarField::Internal& mDotBubble() const;

                    //- ML: Return source term for Lagrangian Liquid Volume Fraction
                    inline volScalarField::Internal& Dbetal();

                    //- ML: Return source term for Lagrangian Liquid Volume Fraction
                    inline const volScalarField::Internal& Dbetal() const;

                    //- Return tmp mass source for kinematicBubble
                    inline tmp<volScalarField::Internal> Srhok() const;

                    //- Return tmp momentum source term (compressible)
                    inline tmp<fvVectorMatrix> SU
                    (
                        volVectorField& U,
                        bool incompressible = false
                    ) const;

                    //- Return tmp mass rate
                    inline tmp<volScalarField::Internal> SU_mDotBubble() const;


        // Check

            //- Total number of parcels
            virtual label nParcels() const
            {
                return CloudType::nParcels();
            }

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Average particle per parcel
            inline scalar totalParticlePerParcel() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            //- Total rotational kinetic energy in the system
            inline scalar rotationalKineticEnergyOfSystem() const;

            //- Mean diameter Dij
            inline scalar Dij(const label i, const label j) const;

            //- Max diameter
            inline scalar Dmax() const;


            // Fields

                //- Volume swept rate of parcels per cell
                inline const tmp<volScalarField> vDotSweep() const;

                //- Return the particle volume fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> theta() const;

                //- Return the particle mass fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> alpha() const;

                //- Return the particle effective density field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> rhoEff() const;


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Relax field
            template<class Type>
            void relax
            (
                DimensionedField<Type, volMesh>& field,
                const DimensionedField<Type, volMesh>& field0,
                const word& name
            ) const;

            //- Scale field
            template<class Type>
            void scale
            (
                DimensionedField<Type, volMesh>& field,
                const word& name
            ) const;

            //- Apply relaxation to (steady state) cloud sources
            void relaxSources(const KinematicBubbleCloud<CloudType>& cloudOldTime);

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Pre-evolve
            void preEvolve
            (
                const typename parcelType::trackingData& td
            );

            //- Evolve the cloud
            void evolve();

            //- Particle motion
            template<class TrackCloudType>
            void motion
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );

            //- Calculate the patch normal and velocity to interact with,
            //  accounting for patch motion if required.
            void patchData
            (
                const parcelType& p,
                const polyPatch& pp,
                vector& normal,
                vector& Up
            ) const;


        // Mapping

            //- Update mesh
            void updateMesh();

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Print cloud information
            void info();

            //- Read particle fields from objects in the obr registry
            virtual void readObjects(const objectRegistry& obr);

            //- Write particle fields as objects into the obr registry
            virtual void writeObjects(objectRegistry& obr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KinematicBubbleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KinematicBubbleCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
