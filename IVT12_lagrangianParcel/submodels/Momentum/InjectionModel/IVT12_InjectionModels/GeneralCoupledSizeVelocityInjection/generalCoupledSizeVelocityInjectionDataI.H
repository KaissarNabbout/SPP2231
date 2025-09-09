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

#include "generalCoupledSizeVelocityInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Access

inline Foam::scalar Foam::generalCoupledSizeVelocityInjectionData::dMain() const
{
    return dMain_;
}

inline Foam::scalar Foam::generalCoupledSizeVelocityInjectionData::dAux() const
{
    return dAux_;
}

inline Foam::scalarList Foam::generalCoupledSizeVelocityInjectionData::VxMain() const
{
    return VxMain_;
}

inline Foam::scalarList Foam::generalCoupledSizeVelocityInjectionData::VxAux() const
{
    return VxAux_;
}

inline Foam::scalarList Foam::generalCoupledSizeVelocityInjectionData::VyMain() const
{
    return VyMain_;
}

inline Foam::scalarList Foam::generalCoupledSizeVelocityInjectionData::VyAux() const
{
    return VyAux_;
}

inline Foam::scalarList Foam::generalCoupledSizeVelocityInjectionData::VzMain() const
{
    return VzMain_;
}

inline Foam::scalarList Foam::generalCoupledSizeVelocityInjectionData::VzAux() const
{
    return VzAux_;
}

//Edit

inline Foam::scalar& Foam::generalCoupledSizeVelocityInjectionData::dMain()
{
    return dMain_;
}

inline Foam::scalar& Foam::generalCoupledSizeVelocityInjectionData::dAux()
{
    return dAux_;
}

inline Foam::scalarList& Foam::generalCoupledSizeVelocityInjectionData::VxMain()
{
    return VxMain_;
}

inline Foam::scalarList& Foam::generalCoupledSizeVelocityInjectionData::VxAux()
{
    return VxAux_;
}

inline Foam::scalarList& Foam::generalCoupledSizeVelocityInjectionData::VyMain()
{
    return VyMain_;
}

inline Foam::scalarList& Foam::generalCoupledSizeVelocityInjectionData::VyAux()
{
    return VyAux_;
}

inline Foam::scalarList& Foam::generalCoupledSizeVelocityInjectionData::VzMain()
{
    return VzMain_;
}

inline Foam::scalarList& Foam::generalCoupledSizeVelocityInjectionData::VzAux()
{
    return VzAux_;
}


// ************************************************************************* //
