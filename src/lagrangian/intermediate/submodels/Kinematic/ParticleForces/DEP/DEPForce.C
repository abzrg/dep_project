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

#include "DEPForce.H"
#include "demandDrivenData.H"
#include "electromagneticConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/**
 * 'coeffs()'
 *     Return the force coefficients dictionary
 *     SIGNATURE: inline const dictionary& coeffs() const;
 *     FILE: $FOAM_SRC/lagrangian/intermediate/submodels/Kinematic/ParticleForces/ParticleForce/ParticleForce.H
 *
 * 'getScalar()'
 *     a macro-generated function; equivalent to dictionary::get<scalar>
 *     SIGNATURE: template<class T> T get (const word& keyword, enum keyType::option matchOpt = keyType::REGEX) const;
 *     FILE: $FOAM_SRC/OpenFOAM/db/dictionary/dictionary.H
 */

template<class CloudType>
Foam::DEPForce<CloudType>::DEPForce
(
    CloudType& owner, // Reference to owner cloud
    const fvMesh& mesh, // Reference to the mesh data base
    const dictionary& dict // Reference to the dictionary that contains particle forces coefficients
)
:
    // Call to the base class constructor
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    // Initializing the private data member
    // Remember to initialize in the same order they have been declared
    EdotGradEName_
    (
        this->coeffs().template getOrDefault<word>("EdotGradE", "EdotGradE")
    ),
    EdotGradEInterpPtr_(nullptr),
    epsilonm_
    (
        this->coeffs().getScalar("epsilonm")
    ),
    corrFactor_
    (
        this->coeffs().getScalar("corrFactor")
    ),
    CMFactor_
    (
        this->coeffs().getScalar("CMFactor")
    )
{}


template<class CloudType>
Foam::DEPForce<CloudType>::DEPForce
(
    const DEPForce& df
)
:
    // Call to the base class copy constructor (possible slicing?)
    ParticleForce<CloudType>(df),
    // Initializing the derived class private members
    EdotGradEName_(df.EdotGradEName_),
    EdotGradEInterpPtr_(df.EdotGradEInterpPtr_),
    epsilonm_(df.epsilonm_),
    corrFactor_(df.corrFactor_),
    CMFactor_(df.CMFactor_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

// Apparently data members can clean themselves according to RAII
template<class CloudType>
Foam::DEPForce<CloudType>::~DEPForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::DEPForce<CloudType>::cacheFields(const bool store)
{
    if (store)
    {
        const volVectorField& EdotGradE =
            this->mesh().template lookupObject<volVectorField>(EdotGradEName_);

        EdotGradEInterpPtr_ = interpolation<vector>::New
        (
            this->owner().solution().interpolationSchemes(),
            EdotGradE
        ).ptr();
    }
    else
    {
        deleteDemandDrivenData(EdotGradEInterpPtr_);
    }
}


template<class CloudType>
Foam::forceSuSp Foam::DEPForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt, // Current time interval
    const scalar mass, // Particle mass
    const scalar Re, // Reynolds number
    const scalar muc // Continuous phase viscosity
) const
{
    /**
     * Construct zero-initialized content
     *
     * from: zero.H (class zero)
     * Zero: Global zero (0)
     *     static constexpr const zero Zero;
     */
    forceSuSp value(Zero);

    const interpolation<vector>& EdotGradEInterp = *EdotGradEInterpPtr_;

    value.Su() =
       corrFactor_
      *constant::mathematical::twoPi                            // twoPi = 2.0 * constant::mathematical::pi
      *epsilonm_*constant::electromagnetic::epsilon0.value()    // epsilon0: dimensionedScalar
      *CMFactor_
      *pow(0.5*p.d(), 3)
      *EdotGradEInterp.interpolate(p.coordinates(), p.currentTetIndices());

    return value;
}


// ************************************************************************* //
