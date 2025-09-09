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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DispersionRASModel<CloudType>::kcInterp() const
{
    if (!kcInterp_.valid())
    {
        FatalErrorInFunction
            << "Carrier phase kc interpolation object not set. "
               "Possible cause: use of dispersion with laminar simulation."
            << abort(FatalError);
    }
      
    return kcInterp_();
}

template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::DispersionRASModel<CloudType>::epsiloncInterp() const
{
    if (!epsiloncInterp_.valid())
    {
        FatalErrorInFunction
            << "Carrier phase epsilonc interpolation object not set. "
               "Possible cause: use of dispersion with laminar simulation."
            << abort(FatalError);
    }
      
    return epsiloncInterp_();
}

// ************************************************************************* //
