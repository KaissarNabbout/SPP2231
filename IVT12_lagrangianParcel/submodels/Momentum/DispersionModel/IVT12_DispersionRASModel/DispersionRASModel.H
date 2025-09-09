/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::DispersionRASModel

Description
    Base class for particle dispersion models based on RAS turbulence.

\*---------------------------------------------------------------------------*/

#ifndef DispersionRASModel_H
#define DispersionRASModel_H

#include "DispersionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class DispersionRASModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DispersionRASModel
:
    public DispersionModel<CloudType>
{
protected:

    // Protected data

        // Locally cached turbulence fields

            //- Turbulence k
            const volScalarField* kPtr_;

            //- Take ownership of the k field
            mutable bool ownK_;

            //- Turbulence epsilon
            const volScalarField* epsilonPtr_;

            //- Take ownership of the epsilon field
            mutable bool ownEpsilon_;

        // Interpolators for continuous phase fields

            //- Turbulence k interpolator
            autoPtr<interpolation<scalar>> kcInterp_;

            //- Turbulence epsilon interpolator
            autoPtr<interpolation<scalar>> epsiloncInterp_;

    // Protected Functions

        //- Return the k field from the turbulence model
        tmp<volScalarField> kModel() const;

        //- Return the epsilon field from the turbulence model
        tmp<volScalarField> epsilonModel() const;


public:

    //- Runtime type information
    TypeName("dispersionRASModel");


    // Constructors

        //- Construct from components
        DispersionRASModel(const dictionary& dict, CloudType& owner);

        //- Construct copy
        DispersionRASModel(const DispersionRASModel<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~DispersionRASModel();


    // Member Functions

        // Access

            //- Return const access to the interpolator for continuous
            //  phase turbulent kinematic energy field
            inline const interpolation<scalar>& kcInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase turbulence dissipation energy field
            inline const interpolation<scalar>& epsiloncInterp() const;

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
        ) = 0;

        //- Cache carrier fields
        virtual void cacheFields(const bool store);


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DispersionRASModelI.H"

#ifdef NoRepository
    #include "DispersionRASModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
