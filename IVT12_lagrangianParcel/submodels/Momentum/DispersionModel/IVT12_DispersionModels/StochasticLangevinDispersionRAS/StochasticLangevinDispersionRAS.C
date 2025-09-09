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

#include "StochasticLangevinDispersionRAS.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticLangevinDispersionRAS<CloudType>::StochasticLangevinDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner)
{}


template<class CloudType>
Foam::StochasticLangevinDispersionRAS<CloudType>::StochasticLangevinDispersionRAS
(
    const StochasticLangevinDispersionRAS<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticLangevinDispersionRAS<CloudType>::~StochasticLangevinDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::StochasticLangevinDispersionRAS<CloudType>::update
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const label celli,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb
)
{
    distributions::standardNormal& stdNormal = this->owner().stdNormal();
    
    const scalar k =
        this->kcInterp().interpolate
        (
            p.coordinates(),
            p.currentTetIndices(td.mesh)
        );
    
    const scalar epsilon =
        this->epsiloncInterp().interpolate
        (
            p.coordinates(),
            p.currentTetIndices(td.mesh)
        );
   
    // Coefficients of the method
    const scalar cT = 0.24;
    const scalar cL = 3.0;
    const vector unitVector(1.0, 1.0, 1.0);
    
    // Sigma is the RMS value of the fluid velocity (here simplified
    // with the assumption of isotropic turbulence)
    const scalar sigma = sqrt(2.0*k/3.0);
    
    // Relative velocity between fluid and parcel
    const vector Urel = Uc - U;
    
    // Relative displacement between fluid parcel
    const vector deltaR = Urel*dt;
    const scalar magDeltaR = mag(deltaR);
    
    // Lagrangian integral time scale
    const scalar TL = cT*sqr(sigma)/(epsilon + rootVSmall);
    
    // Eulerian integral length scale (turbulence)
    const scalar LE = cL*TL*sigma;
    
    // Correlation coefficients
    const scalar fDeltaR = exp(-magDeltaR/(LE + rootVSmall));
    const scalar gDeltaR = fDeltaR*(1.0 - magDeltaR/(2.0*LE + rootVSmall));
        
    // Lagrangian component of the correction function RP
    const scalar RL = exp(-dt/TL);
    
    // Eulerian component of the correlation function RP
    const vector sqrDeltaR = cmptMultiply(deltaR,deltaR);
    const vector RE = (fDeltaR - gDeltaR)*sqrDeltaR/sqr(magDeltaR) + 
                    gDeltaR*unitVector;
    
    // Correction function
    const vector RP = RL*RE;
    
    // Independent Wiener processes with zero mean
    const vector dW(stdNormal.sample(),stdNormal.sample(),stdNormal.sample());
    
    // Update velocity fluctuation
    const vector sqrRP = unitVector - cmptMultiply(RP,RP);
    const vector sqrtRP(sqrt(sqrRP.x()), sqrt(sqrRP.y()), sqrt(sqrRP.z()));   
    UTurb = cmptMultiply(RP,UTurb) + sigma*cmptMultiply(sqrtRP,dW);

    return Uc + UTurb;
}


// ************************************************************************* //
