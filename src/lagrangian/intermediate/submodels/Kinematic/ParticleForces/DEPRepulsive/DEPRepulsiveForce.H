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

Class
    Foam::DEPRepulsiveForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Calculates a specialized Body Force that prevents particles to hit the
    channel wall when they experience a Dielectrophoretic forces in the channel

    This implementation is a hack. There are things that can be done to make it
    a bit better

SourceFiles
    DEPRepulsiveForceI.H
    DEPRepulsiveForce.C

\*---------------------------------------------------------------------------*/

#ifndef DEPRepulsiveForce_H
#define DEPRepulsiveForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class DEPRepulsiveForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DEPRepulsiveForce
:
    public ParticleForce<CloudType>
{
    // Private data

        //- Channel height
        const scalar channelheight_;

        //- Channel width
        const scalar channelwidth_;

        //- Type of the channel. options: straight (default: straight)
        const word channeltype_;

        //- Repulsive force
        const scalar repulsiveForce_;


public:

    //- Runtime type information
    //  The name that is used in the
    //  kinematicProperties.subModels.particleForces subdictionary
    TypeName("DEPRepulsive");


    // Constructors

        //- Construct from mesh
        DEPRepulsiveForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        DEPRepulsiveForce(const DEPRepulsiveForce& bf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new DEPRepulsiveForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DEPRepulsiveForce();


    // Member Functions

        // Access

            //- Return the channel height
            scalar channelheight() const;

            //- Return the channel width
            scalar channelwidth() const;

            //- Return the type of the channel
            const word& channeltype() const;

            //- Return the repulsive force
            scalar repulsiveForce() const;


        // Evaluation

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DEPRepulsiveForceI.H"

#ifdef NoRepository
    #include "DEPRepulsiveForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
