/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "fvmSup.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::KinematicBubbleCloud<CloudType>&
Foam::KinematicBubbleCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType>
inline const Foam::fvMesh& Foam::KinematicBubbleCloud<CloudType>::mesh() const
{
    return mesh_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicBubbleCloud<CloudType>::particleProperties() const
{
    return particleProperties_;
}


template<class CloudType>
inline const Foam::IOdictionary&
Foam::KinematicBubbleCloud<CloudType>::outputProperties() const
{
    return outputProperties_;
}


template<class CloudType>
inline Foam::IOdictionary& Foam::KinematicBubbleCloud<CloudType>::outputProperties()
{
    return outputProperties_;
}


template<class CloudType>
inline const Foam::cloudSolution&
Foam::KinematicBubbleCloud<CloudType>::solution() const
{
    return solution_;
}


template<class CloudType>
inline Foam::cloudSolution& Foam::KinematicBubbleCloud<CloudType>::solution()
{
    return solution_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::KinematicBubbleCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::KinematicBubbleCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::dictionary&
Foam::KinematicBubbleCloud<CloudType>::subModelProperties() const
{
    return subModelProperties_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicBubbleCloud<CloudType>::rho() const
{
    return rho_;
}


template<class CloudType>
inline const Foam::volVectorField& Foam::KinematicBubbleCloud<CloudType>::U() const
{
    return U_;
}


template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicBubbleCloud<CloudType>::mu() const
{
    return mu_;
}

//- ML: Return pressure of carrier
template<class CloudType>
inline const Foam::volScalarField& Foam::KinematicBubbleCloud<CloudType>::p() const
{
    return p_;
}

//- ML: Return Eulerian Liquid Volume Fraction
template<class CloudType>
inline Foam::volScalarField& Foam::KinematicBubbleCloud<CloudType>::alphal()
{
    return alphal_;
}

//- ML: Return Lagrangian Liquid Volume Fraction
template<class CloudType>
inline Foam::volScalarField& Foam::KinematicBubbleCloud<CloudType>::betal()
{
    return betal_;
}

//- ML: Return  old Lagrangian Liquid Volume Fraction
template<class CloudType>
inline Foam::volScalarField& Foam::KinematicBubbleCloud<CloudType>::betal_old()
{
    return betal_old_;
}

template<class CloudType>
inline const Foam::dimensionedVector& Foam::KinematicBubbleCloud<CloudType>::g() const
{
    return g_;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicBubbleCloud<CloudType>::pAmbient() const
{
    return pAmbient_;
}


template<class CloudType>
inline Foam::scalar& Foam::KinematicBubbleCloud<CloudType>::pAmbient()
{
    return pAmbient_;
}


template<class CloudType>
//inline const typename CloudType::parcelType::forceType&
inline const typename Foam::KinematicBubbleCloud<CloudType>::forceType&
Foam::KinematicBubbleCloud<CloudType>::forces() const
{
    return forces_;
}


template<class CloudType>
inline typename Foam::KinematicBubbleCloud<CloudType>::forceType&
Foam::KinematicBubbleCloud<CloudType>::forces()
{
    return forces_;
}


template<class CloudType>
inline typename Foam::KinematicBubbleCloud<CloudType>::functionType&
Foam::KinematicBubbleCloud<CloudType>::functions()
{
    return functions_;
}


template<class CloudType>
inline const Foam::InjectionModelList<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::injectors() const
{
    return injectors_;
}


template<class CloudType>
inline Foam::InjectionModelList<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::injectors()
{
    return injectors_;
}


template<class CloudType>
inline const Foam::DispersionModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::dispersion() const
{
    return *dispersionModel_;
}


template<class CloudType>
inline Foam::DispersionModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::dispersion()
{
    return *dispersionModel_;
}


template<class CloudType>
inline const Foam::PatchInteractionModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::patchInteraction() const
{
    return *patchInteractionModel_;
}


template<class CloudType>
inline Foam::PatchInteractionModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::patchInteraction()
{
    return *patchInteractionModel_;
}


template<class CloudType>
inline const Foam::StochasticCollisionModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::stochasticCollision() const
{
    return *stochasticCollisionModel_;
}


template<class CloudType>
inline Foam::StochasticCollisionModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::stochasticCollision()
{
    return *stochasticCollisionModel_;
}


template<class CloudType>
inline const Foam::SurfaceFilmModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::surfaceFilm() const
{
    return *surfaceFilmModel_;
}


template<class CloudType>
inline Foam::SurfaceFilmModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::surfaceFilm()
{
    return *surfaceFilmModel_;
}


template<class CloudType>
inline const Foam::PackingModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::packingModel() const
{
    return *packingModel_;
}


template<class CloudType>
inline Foam::PackingModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::packingModel()
{
    return *packingModel_;
}


template<class CloudType>
inline const Foam::DampingModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::dampingModel() const
{
    return *dampingModel_;
}


template<class CloudType>
inline Foam::DampingModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::dampingModel()
{
    return *dampingModel_;
}


template<class CloudType>
inline const Foam::IsotropyModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::isotropyModel() const
{
    return *isotropyModel_;
}


template<class CloudType>
inline Foam::IsotropyModel<Foam::KinematicBubbleCloud<CloudType>>&
Foam::KinematicBubbleCloud<CloudType>::isotropyModel()
{
    return *isotropyModel_;
}


template<class CloudType>
inline const Foam::integrationScheme&
Foam::KinematicBubbleCloud<CloudType>::UIntegrator() const
{
    return *UIntegrator_;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicBubbleCloud<CloudType>::massInSystem() const
{
    scalar sysMass = 0.0;
    for (const parcelType& p : *this)
    {
        sysMass += p.nParticle()*p.mass();
    }

    return sysMass;
}


template<class CloudType>
inline Foam::vector
Foam::KinematicBubbleCloud<CloudType>::linearMomentumOfSystem() const
{
    vector linearMomentum(Zero);

    for (const parcelType& p : *this)
    {
        linearMomentum += p.nParticle()*p.mass()*p.U();
    }

    return linearMomentum;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicBubbleCloud<CloudType>::totalParticlePerParcel() const
{
    scalar parPerParcel = 0;

    for (const parcelType& p : *this)
    {
        parPerParcel += p.nParticle();
    }

    return parPerParcel;
}


template<class CloudType>
inline Foam::scalar
Foam::KinematicBubbleCloud<CloudType>::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0;

    for (const parcelType& p : *this)
    {
        linearKineticEnergy += p.nParticle()*0.5*p.mass()*(p.U() & p.U());
    }

    return linearKineticEnergy;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicBubbleCloud<CloudType>::Dij
(
    const label i,
    const label j
) const
{
    scalar si = 0.0;
    scalar sj = 0.0;
    for (const parcelType& p : *this)
    {
        si += p.nParticle()*pow(p.d(), i);
        sj += p.nParticle()*pow(p.d(), j);
    }

    reduce(si, sumOp<scalar>());
    reduce(sj, sumOp<scalar>());
    sj = max(sj, VSMALL);

    return si/sj;
}


template<class CloudType>
inline Foam::scalar Foam::KinematicBubbleCloud<CloudType>::Dmax() const
{
    scalar d = -GREAT;
    for (const parcelType& p : *this)
    {
        d = max(d, p.d());
    }

    reduce(d, maxOp<scalar>());

    return max(0.0, d);
}


template<class CloudType>
inline Foam::Random& Foam::KinematicBubbleCloud<CloudType>::rndGen() const
{
    return rndGen_;
}


template<class CloudType>
inline Foam::List<Foam::DynamicList<typename CloudType::particleType*>>&
Foam::KinematicBubbleCloud<CloudType>::cellOccupancy()
{
    if (!cellOccupancyPtr_)
    {
        buildCellOccupancy();
    }

    return *cellOccupancyPtr_;
}


template<class CloudType>
inline const Foam::scalarField&
Foam::KinematicBubbleCloud<CloudType>::cellLengthScale() const
{
    return cellLengthScale_;
}


template<class CloudType>
inline void Foam::KinematicBubbleCloud<CloudType>::transferToCarrier
(
    const parcelType& p,
    const typename parcelType::trackingData& td
)
{
    const scalar m = p.nParticle()*p.mass();

    rhokTrans()[p.cell()] += m;

    UTrans()[p.cell()] += m*p.U();
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::rhokTrans()
{
    return *rhokTrans_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::rhokTrans() const
{
    return *rhokTrans_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::UTrans()
{
    return *UTrans_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::UTrans() const
{
    return *UTrans_;
}


template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::UCoeff()
{
    return *UCoeff_;
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::UCoeff() const
{
    return *UCoeff_;
}

//- ML: Return rate of change of mass
template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::mDotBubble()
{
    return *mDotBubble_;
}

//- ML: Return rate of change of mass
template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::mDotBubble() const
{
    return *mDotBubble_;
}

//- ML: Return source term for Lagrangian Liquid Volume Fraction
template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::Dbetal()
{
    return *Dbetal_;
}

//- ML: Return source term for Lagrangian Liquid Volume Fraction
template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::KinematicBubbleCloud<CloudType>::Dbetal() const
{
    return *Dbetal_;
}

template<class CloudType>
inline Foam::tmp<Foam::volScalarField::Internal>
Foam::KinematicBubbleCloud<CloudType>::Srhok() const
{
    if (debug)
    {
        Pout<< "rhokTrans min/max = " << min(rhokTrans()).value() << ", "
            << max(rhokTrans()).value() << endl;
    }

    if (this->solution().coupled())
    {
        return rhokTrans()/this->db().time().deltaT()/this->mesh().V();
    }

    return tmp<volScalarField::Internal>::New
    (
        this->newIOobject(IOobject::scopedName(this->name(), "rhokTrans")),
        this->mesh(),
        dimensionedScalar
        (
            rhokTrans().dimensions()/dimTime/dimVolume, Zero
        )
    );
}


template<class CloudType>
inline Foam::tmp<Foam::fvVectorMatrix>
Foam::KinematicBubbleCloud<CloudType>::SU(volVectorField& U, bool incompressible)
const
{
    if (debug)
    {
        Pout<< "UTrans min/max = " << min(UTrans()).value() << ", "
            << max(UTrans()).value() << nl
            << "UCoeff min/max = " << min(UCoeff()).value() << ", "
            << max(UCoeff()).value() << endl;
    }

    dimensionSet dim(dimForce);
    if (incompressible)
    {
        dim.reset(dimForce/dimDensity);
    }

    if (solution_.coupled())
    {
        if (solution_.semiImplicit("U"))
        {
            volScalarField::Internal
                Vdt(mesh_.V()*this->db().time().deltaT());

            if (incompressible)
            {
                Vdt.dimensions() *= dimDensity;
            }

            return UTrans()/Vdt - fvm::Sp(UCoeff()/Vdt, U) + UCoeff()/Vdt*U;
        }
        else
        {
            auto tfvm = tmp<fvVectorMatrix>::New(U, dim);
            auto& fvm = tfvm.ref();

            fvm.source() = -UTrans()/(this->db().time().deltaT());

            return tfvm;
        }
    }

    return tmp<fvVectorMatrix>::New(U, dim);
}

//-ML: Return tmp mass rate of bubble
template<class CloudType>
inline Foam::tmp<Foam::volScalarField::Internal>
Foam::KinematicBubbleCloud<CloudType>::SU_mDotBubble() const
{
    
    return mDotBubble()/this->db().time().deltaT()/this->mesh().V();

}

template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicBubbleCloud<CloudType>::vDotSweep() const
{
    auto tvDotSweep = tmp<volScalarField>::New
    (
        this->newIOobject(IOobject::scopedName(this->name(), "vDotSweep")),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& vDotSweep = tvDotSweep.ref();

    for (const parcelType& p : *this)
    {
        const label celli = p.cell();

        vDotSweep[celli] += p.nParticle()*p.areaP()*mag(p.U() - U_[celli]);
    }

    vDotSweep.primitiveFieldRef() /= mesh_.V();
    vDotSweep.correctBoundaryConditions();

    return tvDotSweep;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicBubbleCloud<CloudType>::theta() const
{
    auto ttheta = tmp<volScalarField>::New
    (
        this->newIOobject(IOobject::scopedName(this->name(), "theta")),
        mesh_,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& theta = ttheta.ref();

    for (const parcelType& p : *this)
    {
        const label celli = p.cell();

        theta[celli] += p.nParticle()*p.volume();
    }

    theta.primitiveFieldRef() /= mesh_.V();
    theta.correctBoundaryConditions();

    return ttheta;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicBubbleCloud<CloudType>::alpha() const
{
    auto talpha = tmp<volScalarField>::New
    (
        this->newIOobject(IOobject::scopedName(this->name(), "alpha")),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    scalarField& alpha = talpha.ref().primitiveFieldRef();
    for (const parcelType& p : *this)
    {
        const label celli = p.cell();

        alpha[celli] += p.nParticle()*p.mass();
    }

    alpha /= (mesh_.V()*rho_);

    return talpha;
}


template<class CloudType>
inline const Foam::tmp<Foam::volScalarField>
Foam::KinematicBubbleCloud<CloudType>::rhoEff() const
{
    auto trhoEff = tmp<volScalarField>::New
    (
        this->newIOobject(IOobject::scopedName(this->name(), "rhoEff")),
        mesh_,
        dimensionedScalar(dimDensity, Zero)
    );

    scalarField& rhoEff = trhoEff.ref().primitiveFieldRef();
    for (const parcelType& p : *this)
    {
        const label celli = p.cell();

        rhoEff[celli] += p.nParticle()*p.mass();
    }

    rhoEff /= mesh_.V();

    return trhoEff;
}


// ************************************************************************* //
