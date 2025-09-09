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

\*---------------------------------------------------------------------------*/

#include "ParticleRelevantTimeList.H"
#include "entry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleRelevantTimeList<CloudType>::ParticleRelevantTimeList
(
    CloudType& owner,
    const fvMesh& mesh
)
:
    PtrList<ParticleRelevantTime<CloudType>>(),
    owner_(owner),
    mesh_(mesh),
    dict_(dictionary::null),
    calcRelevantTime_(true)
{}


template<class CloudType>
Foam::ParticleRelevantTimeList<CloudType>::ParticleRelevantTimeList
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool readFields
)
:
    PtrList<ParticleRelevantTime<CloudType>>(),
    owner_(owner),
    mesh_(mesh),
    dict_(dict),
    calcRelevantTime_(true)
{
    if (readFields)
    {
        wordList modelNames(dict.toc());

        Info<< "Constructing relevant times list" << endl;

        if (modelNames.size() > 0)
        {
            this->setSize(modelNames.size());

            label i = 0;
            forAllConstIter(IDLList<entry>, dict, iter)
            {
                const word& model = iter().keyword();
                if (iter().isDict())
                {
                    this->set
                    (
                        i++,
                        ParticleRelevantTime<CloudType>::New
                        (
                            owner,
                            mesh,
                            iter().dict(),
                            model
                        )
                    );
                }
                else
                {
                    this->set
                    (
                        i++,
                        ParticleRelevantTime<CloudType>::New
                        (
                            owner,
                            mesh,
                            dictionary::null,
                            model
                        )
                    );
                }
            }
        }
        else
        {
            Info<< "    none" << endl;
        }
    }
}


template<class CloudType>
Foam::ParticleRelevantTimeList<CloudType>::ParticleRelevantTimeList
(
    const ParticleRelevantTimeList& prtl
)
:
    PtrList<ParticleRelevantTime<CloudType>>(prtl),
    owner_(prtl.owner_),
    mesh_(prtl.mesh_),
    dict_(prtl.dict_)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleRelevantTimeList<CloudType>::~ParticleRelevantTimeList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleRelevantTimeList<CloudType>::cacheFields(const bool store)
{
    forAll(*this, i)
    {
        this->operator[](i).cacheFields(store);
    }
}


template<class CloudType>
Foam::scalar Foam::ParticleRelevantTimeList<CloudType>::calcRelevantTime
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar Re
) const
{
    scalar dtNew(great);
    scalar dtOld;
    if (calcRelevantTime_)
    {
        forAll(*this, i)
        {
            dtOld = dtNew;
            dtNew = this->operator[](i).calcRelevantTime(p, td, Re);
            dtNew = min(dtNew,dtOld);
        }
    }

    return dtNew/(td.trackTime() + rootVSmall);
}


// ************************************************************************* //
