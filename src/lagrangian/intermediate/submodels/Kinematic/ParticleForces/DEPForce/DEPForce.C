/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

template<class CloudType>
Foam::DEPForce<CloudType>::DEPForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
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
    ParticleForce<CloudType>(df),
    EdotGradEName_(df.EdotGradEName_),
    EdotGradEInterpPtr_(df.EdotGradEInterpPtr_),
    epsilonm_(df.epsilonm_),
    corrFactor_(df.corrFactor_),
    CMFactor_(df.CMFactor_)
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
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero);

    const interpolation<vector>& EdotGradEInterp = *EdotGradEInterpPtr_;

    value.Su() =
       corrFactor_
      *constant::mathematical::twoPi
      *epsilonm_*constant::electromagnetic::epsilon0.value()
      *CMFactor_
      *pow(0.5*p.d(), 3)
      *EdotGradEInterp.interpolate(p.coordinates(), p.currentTetIndices());

    return value;
}


// ************************************************************************* //
