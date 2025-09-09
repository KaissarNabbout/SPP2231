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

\*---------------------------------------------------------------------------*/

#include "kEpsilonRelevantTimeBPG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::kEpsilonRelevantTimeBPG<CloudType>::kEpsilonRelevantTimeBPG
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleRelevantTime<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::kEpsilonRelevantTimeBPG<CloudType>::kEpsilonRelevantTimeBPG
(
    const kEpsilonRelevantTimeBPG<CloudType>& rt
)
:
    ParticleRelevantTime<CloudType>(rt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::kEpsilonRelevantTimeBPG<CloudType>::~kEpsilonRelevantTimeBPG()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::kEpsilonRelevantTimeBPG<CloudType>::calcRelevantTime
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar Re
) const
{
    const auto& dispersionModel = this->owner().dispersion();
    
    const auto* dispersionRASModelPtr = 
        dynamic_cast<const DispersionRASModel<CloudType>*>(&dispersionModel);
        
    if (!dispersionRASModelPtr)
    {
        FatalErrorInFunction
            << "Failed to cast dispersionModel to DispersionRASModel<CloudType>!"
            << " Possible reason: no dispersion model is being used."
            << exit(FatalError);
    }
        
    const scalar k = dispersionRASModelPtr->kcInterp().interpolate
    (
        p.coordinates(),
        p.currentTetIndices(td.mesh)
    );
        
    const scalar epsilon = dispersionRASModelPtr->epsiloncInterp().interpolate
    (
        p.coordinates(),
        p.currentTetIndices(td.mesh)
    );
            
    return 0.48*k/(sqr(epsilon) + rootVSmall);
}


// ************************************************************************* //
