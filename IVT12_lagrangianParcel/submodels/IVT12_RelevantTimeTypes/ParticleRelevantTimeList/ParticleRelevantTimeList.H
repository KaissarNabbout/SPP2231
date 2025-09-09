/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::ParticleRelevantTimeList

Description
    List of particle relevant times

SourceFiles
    ParticleRelevantTimeListI.H
    ParticleRelevantTimeList.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleRelevantTimeList_H
#define ParticleRelevantTimeList_H

#include "ParticleRelevantTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ParticleRelevantTimeList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleRelevantTimeList
:
    public PtrList<ParticleRelevantTime<CloudType>>
{
    // Private Data

        //- Reference to the owner cloud
        CloudType& owner_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Relevant times dictionary
        const dictionary dict_;

        //- Calculate relevant times flag
        bool calcRelevantTime_;


public:

    // Constructors

        //- Null constructor
        ParticleRelevantTimeList(CloudType& owner, const fvMesh& mesh);

        //- Construct from mesh
        ParticleRelevantTimeList
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool readFields
        );

        //- Construct copy
        ParticleRelevantTimeList(const ParticleRelevantTimeList& prtl);


    //- Destructor
    virtual ~ParticleRelevantTimeList();


    // Member Functions

        // Access

            //- Return const access to the cloud owner
            inline const CloudType& owner() const;

            //- Return references to the cloud owner
            inline CloudType& owner();

            //- Return the mesh database
            inline const fvMesh& mesh() const;

            //- Return the relevant times dictionary
            inline const dictionary& dict() const;

            //- Set the relevant time flag
            inline void setCalcRelevantTime(bool flag);


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the relevant times
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

#include "ParticleRelevantTimeListI.H"

#ifdef NoRepository
    #include "ParticleRelevantTimeList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
