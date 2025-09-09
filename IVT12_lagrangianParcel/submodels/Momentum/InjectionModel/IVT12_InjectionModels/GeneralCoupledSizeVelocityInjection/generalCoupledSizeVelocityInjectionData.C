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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(generalCoupledSizeVelocityInjectionData, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::generalCoupledSizeVelocityInjectionData::generalCoupledSizeVelocityInjectionData()
:
    dMain_(),
    dAux_(),
    VxMain_(),
    VxAux_(),
    VyMain_(),
    VyAux_(),
    VzMain_(),
    VzAux_()  
{}


Foam::generalCoupledSizeVelocityInjectionData::generalCoupledSizeVelocityInjectionData
(
    const dictionary& dict
)
:
    dMain_(readScalar(dict.lookup("dClass"))),
    dAux_(readScalar(dict.lookup("dAux"))),
    VxMain_(dict.lookup("VxClass")),
    VxAux_(dict.lookup("VxProb")),
    VyMain_(dict.lookup("VyClass")),
    VyAux_(dict.lookup("VyProb")),
    VzMain_(dict.lookup("VzClass")),
    VzAux_(dict.lookup("VzProb"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::generalCoupledSizeVelocityInjectionData::~generalCoupledSizeVelocityInjectionData()
{}


// ************************************************************************* //
