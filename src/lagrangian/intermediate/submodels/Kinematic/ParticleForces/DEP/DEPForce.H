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
    Foam::DEPForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Calculates particle Dielectrophoretic (DEP) force

SourceFiles
    DEPI.H
    DEP.C

\*---------------------------------------------------------------------------*/

#ifndef DEP_H
#define DEP_H

#include "ParticleForce.H"  // Abstract base class for particle forces
#include "interpolation.H"  // TODO: Used for what ?

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class DEPForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DEPForce
:
    public ParticleForce<CloudType>
{
    // Private data

        //- Name of electric field strength field - default = "EdotGradE"
        const word EdotGradEName_;

        //- EdotGradE interpolator
        const interpolation<vector>* EdotGradEInterpPtr_;

        //- Relative permitivity of the medium/fluid, []
        const scalar epsilonm_;

        //- Correction factor (to account for what?), []
        const scalar corrFactor_;

        //- Clausius-Mossotti factor, []
        const scalar CMFactor_;


public:

    //- Runtime type information
    TypeName("DEP");


    // Constructors

        //- Construct from mesh
        DEPForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy (df: DEP Force)
        DEPForce(const DEPForce& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new DEPForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DEPForce();


    // Member Functions

        // Access

            //- Return the name of electric field strength field
            const word& EdotGradEName() const;

            //- Return Relative permittivity of the medium/fluid
            scalar epsilonm() const;

            //- Return Correction factor
            scalar corrFactor() const;

            //- Return Clausius Mossotti factor
            scalar CMFactor() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

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

#include "DEPForceI.H"

#ifdef NoRepository
    #include "DEPForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
