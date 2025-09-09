/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "SchillerNaumannStokesDragForce.H"

// * * * * * * * * * * * * *  Static Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::SchillerNaumannStokesDragForce<CloudType>::CdRe(const scalar Re)
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
Foam::SchillerNaumannStokesDragForce<CloudType>::SchillerNaumannStokesDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::SchillerNaumannStokesDragForce<CloudType>::SchillerNaumannStokesDragForce
(
    const SchillerNaumannStokesDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SchillerNaumannStokesDragForce<CloudType>::~SchillerNaumannStokesDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::SchillerNaumannStokesDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    // Remove the "const" protection of "p" to be able to change the value of Cd.
    // This is a lazy implementation.
    // Not a good practice, but it should work fine here.
    const_cast<typename CloudType::parcelType&>(p).Cd() = CdRe(Re) / Re;
    
    return forceSuSp(Zero, mass*0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d())));
}


// ************************************************************************* //
