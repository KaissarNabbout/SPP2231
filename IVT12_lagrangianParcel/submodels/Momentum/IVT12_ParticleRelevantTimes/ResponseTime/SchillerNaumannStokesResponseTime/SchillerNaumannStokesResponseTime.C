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

#include "SchillerNaumannStokesResponseTime.H"

// * * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::SchillerNaumannStokesResponseTime<CloudType>::CdRe(const scalar Re)
{
    if (Re < 0.5)
    {
        return 24.0;
    }
    else if (Re > 1000.0)
    {
        return 0.44*Re;
    }
    else
    {
        return 24.0*(1.0 + 0.15*pow(Re, 0.687));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SchillerNaumannStokesResponseTime<CloudType>::SchillerNaumannStokesResponseTime
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleRelevantTime<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::SchillerNaumannStokesResponseTime<CloudType>::SchillerNaumannStokesResponseTime
(
    const SchillerNaumannStokesResponseTime<CloudType>& rt
)
:
    ParticleRelevantTime<CloudType>(rt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SchillerNaumannStokesResponseTime<CloudType>::~SchillerNaumannStokesResponseTime()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::SchillerNaumannStokesResponseTime<CloudType>::calcRelevantTime
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar Re
) const
{
    return (4.0*p.rho()*sqr(p.d())) / (3.0*td.muc()*CdRe(Re) + small);
}


// ************************************************************************* //
