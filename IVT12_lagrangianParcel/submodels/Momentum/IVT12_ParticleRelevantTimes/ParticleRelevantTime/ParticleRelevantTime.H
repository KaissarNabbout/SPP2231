/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::ParticleRelevantTime

Description
    Abstract base class for particle relevant times

SourceFiles
    ParticleRelevantTimeI.H
    ParticleRelevantTime.C
    ParticleRelevantTimeNew.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleRelevantTime_H
#define ParticleRelevantTime_H

#include "dictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ParticleRelevantTime Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleRelevantTime
{
    // Private Data

        //- Reference to the owner cloud
        CloudType& owner_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Relevant times coefficients dictionary
        const dictionary coeffs_;


public:

    //- Runtime type information
    TypeName("ParticleRelevantTime");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ParticleRelevantTime,
        dictionary,
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (owner, mesh, dict)
    );

    // Constructors

        //- Construct from mesh
        ParticleRelevantTime
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& relevantTimeType,
            const bool readCoeffs
        );

        //- Construct copy
        ParticleRelevantTime(const ParticleRelevantTime& prt);

        //- Construct and return a clone
        virtual autoPtr<ParticleRelevantTime<CloudType>> clone() const
        {
            return autoPtr<ParticleRelevantTime<CloudType>>
            (
                new ParticleRelevantTime<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleRelevantTime();


    //- Selector
    static autoPtr<ParticleRelevantTime<CloudType>> New
    (
        CloudType& owner,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& relevantTimeType
    );


    // Member Functions

        // Access

            //- Return const access to the cloud owner
            inline const CloudType& owner() const;

            //- Return references to the cloud owner
            inline CloudType& owner();

            //- Return the mesh database
            inline const fvMesh& mesh() const;

            //- Return the relevant times coefficients dictionary
            inline const dictionary& coeffs() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the relevant time
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

#include "ParticleRelevantTimeI.H"

#ifdef NoRepository
    #include "ParticleRelevantTime.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParticleRelevantTimeModel(CloudType)                                      \
                                                                               \
    typedef Foam::CloudType::momentumCloudType momentumCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::ParticleRelevantTime<momentumCloudType>, 0);                          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            ParticleRelevantTime<momentumCloudType>,                                 \
            dictionary                                                         \
        );                                                                     \
    }


#define makeParticleRelevantTimeModelType(SS, CloudType)                              \
                                                                               \
    typedef Foam::CloudType::momentumCloudType momentumCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<momentumCloudType>, 0);      \
                                                                               \
    Foam::ParticleRelevantTime<momentumCloudType>::                                  \
        adddictionaryConstructorToTable<Foam::SS<momentumCloudType>>          \
            add##SS##CloudType##momentumCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
