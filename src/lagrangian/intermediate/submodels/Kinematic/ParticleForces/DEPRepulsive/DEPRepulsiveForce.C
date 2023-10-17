/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "DEPRepulsiveForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DEPRepulsiveForce<CloudType>::DEPRepulsiveForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    // true: it read coefficients
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    channelheight_
    (
        this->coeffs().getScalar("channelheight")
    ),
    channelwidth_
    (
        this->coeffs().getScalar("channelwidth")
    ),
    channeltype_
    (
        this->coeffs().template getOrDefault<word>("channeltype", "straight")
    ),
    repulsiveForce_
    (
        this->coeffs().getScalar("repulsiveForce")
    )
{}


template<class CloudType>
Foam::DEPRepulsiveForce<CloudType>::DEPRepulsiveForce(const DEPRepulsiveForce& bf)
:
    ParticleForce<CloudType>(bf),
    channelheight_(bf.channelheight_),
    channelwidth_(bf.channelwidth_),
    channeltype_(bf.channeltype_),
    repulsiveForce_(bf.repulsiveForce_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DEPRepulsiveForce<CloudType>::~DEPRepulsiveForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::DEPRepulsiveForce<CloudType>::calcNonCoupled
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

    // Get particle position
    scalar pos_x = p.position().x();
    scalar pos_y = p.position().y();
    scalar pos_z = p.position().z();

    // Get the particle radius
    scalar pRadius = p.d() / 2;

    if (channeltype_ == "straight")
    {
        // down and up patch

        // If particle is too close to the left wall (z = channelheight)
        if (pos_z > (channelheight_ - pRadius))
        {
            value.Su() = vector (0, 0, -repulsiveForce_);
        }
        // If particle is too close the right wall (z = 0)
        else if (pos_z < pRadius)
        {
            value.Su() = vector(0, 0, repulsiveForce_);
        }
        // If particle is too close to left wall (y = channelwdith)
        else if (pos_y > (channelwidth_ - pRadius))
        {
            value.Su() = vector(0, -repulsiveForce_, 0);
        }
        // If particle is too close to right wall (y = channelwdith)
        else if (pos_y < pRadius)
        {
            value.Su() = vector(0, repulsiveForce_, 0);
        }
        // If particle is somewhere in between
        else
        {
            value.Su() = vector(0, 0, 0);
        }
    }
    else if (channeltype_ == "curved") // TODO: make it better + replae magic numbers
    {
        if (p.d() == 0.000005 && pos_x > 0.0000975  && pos_y > 0.000010)
        {
            value.Su() = vector (-repulsiveForce_, 0, 0);
        }
        else if (p.d() == 0.000006 && pos_x > 0.000097  && pos_y > 0.000010)
        {
          value.Su() = vector (-repulsiveForce_,0,0);
        }
        else if (p.d() == 0.000007 && pos_x > 0.0000965  && pos_y > 0.000010)
        {
          value.Su() = vector (-repulsiveForce_,0,0);
        }
        else if (p.d() == 0.000008 && pos_x > 0.000096  && pos_y > 0.000010)
        {
          value.Su() = vector (-repulsiveForce_,0,0);
        }
        else if (p.d() == 0.000009 && pos_x > 0.0000955  && pos_y > 0.000010)
        {
          value.Su() = vector (-repulsiveForce_,0,0);
        }
        else if (p.d() == 0.000010 && pos_x > 0.000095  && pos_y > 0.000010)
        {
          value.Su() = vector (-repulsiveForce_,0,0);
        }
        else
        {
          value.Su() = vector(0,0,0);
        }
    }
    else
    {
        value.Su() = vector(0,0,0);
    }

    return value;
}


// ************************************************************************* //
