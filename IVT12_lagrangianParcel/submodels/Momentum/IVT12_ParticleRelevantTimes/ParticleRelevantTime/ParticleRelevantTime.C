/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "ParticleRelevantTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleRelevantTime<CloudType>::ParticleRelevantTime
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& relevantTimeType,
    const bool readCoeffs
)
:
    owner_(owner),
    mesh_(mesh),
    coeffs_
    (
        readCoeffs
      ? dict.optionalSubDict(relevantTimeType + "Coeffs")
      : dictionary::null
    )
{
    if (readCoeffs && coeffs_.isNull())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Relevant Time " << relevantTimeType << " must be specified as a dictionary"
            << exit(FatalIOError);
    }
}


template<class CloudType>
Foam::ParticleRelevantTime<CloudType>::ParticleRelevantTime(const ParticleRelevantTime& prt)
:
    owner_(prt.owner_),
    mesh_(prt.mesh_),
    coeffs_(prt.coeffs_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleRelevantTime<CloudType>::~ParticleRelevantTime()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleRelevantTime<CloudType>::cacheFields(const bool store)
{}


template<class CloudType>
Foam::scalar Foam::ParticleRelevantTime<CloudType>::calcRelevantTime
(
    const typename CloudType::parcelType&,
    const typename CloudType::parcelType::trackingData& td,
    const scalar Re
) const
{
    scalar dtNew = 0.0;
    return dtNew;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleRelevantTimeNew.C"

// ************************************************************************* //
