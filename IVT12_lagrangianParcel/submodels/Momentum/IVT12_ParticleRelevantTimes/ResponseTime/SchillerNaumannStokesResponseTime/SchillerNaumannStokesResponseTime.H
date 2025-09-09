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

Class
    Foam::SchillerNaumannStokesResponseTime

Description
    Particle response time based on Schiller-Naumann 
    drag model for spheres.

    Reference:
    \verbatim
        Naumann, Z., & Schiller, L. (1935).
        A drag coefficient correlation.
        Z Ver Deutsch Ing, 77, 318-323.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef SchillerNaumannStokesResponseTime_H
#define SchillerNaumannStokesResponseTime_H

#include "ParticleRelevantTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SchillerNaumannStokesResponseTime Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SchillerNaumannStokesResponseTime
:
    public ParticleRelevantTime<CloudType>
{
public:

    // Static Member Functions

        //- Drag coefficient multiplied by Reynolds number
        static scalar CdRe(const scalar Re);


    //- Runtime type information
    TypeName("SchillerNaumannStokesResponseTime");


    // Constructors

        //- Construct from mesh
        SchillerNaumannStokesResponseTime
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        SchillerNaumannStokesResponseTime(const SchillerNaumannStokesResponseTime<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleRelevantTime<CloudType>> clone() const
        {
            return autoPtr<ParticleRelevantTime<CloudType>>
            (
                new SchillerNaumannStokesResponseTime<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SchillerNaumannStokesResponseTime();


    // Member Functions

        // Evaluation

            //- Calculate the characteristic time
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

#ifdef NoRepository
    #include "SchillerNaumannStokesResponseTime.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
