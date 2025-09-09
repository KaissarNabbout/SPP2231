/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::kEpsilonRelevantTimeBPG

Description
    Lagrangian and Eulerian integral time scalar due to turbulence.

    Reference:
    \verbatim
        Sommerfeld, M. (Ed.). (2008). Best practice guidelines for
        computational fluid dynamics of dispersed multiphase flows. Ercoftac.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonRelevantTimeBPG_H
#define kEpsilonRelevantTimeBPG_H

#include "ParticleRelevantTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class kEpsilonRelevantTimeBPG Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class kEpsilonRelevantTimeBPG
:
    public ParticleRelevantTime<CloudType>
{
public:

    //- Runtime type information
    TypeName("kEpsilonRelevantTimeBPG");


    // Constructors

        //- Construct from mesh
        kEpsilonRelevantTimeBPG
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        kEpsilonRelevantTimeBPG(const kEpsilonRelevantTimeBPG<CloudType>& rt);

        //- Construct and return a clone
        virtual autoPtr<ParticleRelevantTime<CloudType>> clone() const
        {
            return autoPtr<ParticleRelevantTime<CloudType>>
            (
                new kEpsilonRelevantTimeBPG<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~kEpsilonRelevantTimeBPG();


    // Member Functions

        // Evaluation

            //- Calculate the particle relevant times 
            virtual Foam::scalar calcRelevantTime
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar Re
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonRelevantTimeBPG.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
