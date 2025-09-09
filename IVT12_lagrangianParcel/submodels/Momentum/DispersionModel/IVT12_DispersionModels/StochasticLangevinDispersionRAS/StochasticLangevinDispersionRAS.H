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
    Foam::StochasticLangevinDispersionRAS

Description

\*---------------------------------------------------------------------------*/

#ifndef StochasticLangevinDispersionRAS_H
#define StochasticLangevinDispersionRAS_H

#include "DispersionRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class StochasticLangevinDispersionRAS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class StochasticLangevinDispersionRAS
:
    public DispersionRASModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("StochasticLangevinDispersionRAS");
    // Constructors

        //- Construct from components
        StochasticLangevinDispersionRAS(const dictionary& dict, CloudType& owner);

        //- Construct copy
        StochasticLangevinDispersionRAS(const StochasticLangevinDispersionRAS<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const
        {
            return autoPtr<DispersionModel<CloudType>>
            (
                new StochasticLangevinDispersionRAS<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~StochasticLangevinDispersionRAS();


    // Member Functions

        //- Update (disperse particles)
        virtual vector update
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const scalar dt,
            const label celli,
            const vector& U,
            const vector& Uc,
            vector& UTurb,
            scalar& tTurb
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "StochasticLangevinDispersionRAS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
