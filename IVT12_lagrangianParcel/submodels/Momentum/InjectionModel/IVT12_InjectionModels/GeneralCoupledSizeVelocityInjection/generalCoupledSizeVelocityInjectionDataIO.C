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

Foam::generalCoupledSizeVelocityInjectionData::generalCoupledSizeVelocityInjectionData(Istream& is)
{
    is.check("reading dClass");
    is >> dMain_;

    is.check("reading dProb");
    is >> dAux_;
    
    is.check("reading VxClass");
    is >> VxMain_;

    is.check("reading VxProb");
    is >> VxAux_;
    
    is.check("reading VyClass");
    is >> VyMain_;

    is.check("reading VyProb");
    is >> VyAux_;
    
    is.check("reading VzClass");
    is >> VzMain_;

    is.check("reading VzProb");
    is >> VzAux_;

    is.check("generalCoupledSizeVelocityInjectionData(Istream& is)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const generalCoupledSizeVelocityInjectionData& data
)
{
    os << data.dMain_ << data.dAux_
       << data.VxMain_ << data.VxAux_
       << data.VyMain_ << data.VyAux_
       << data.VzMain_ << data.VzAux_;

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, generalCoupledSizeVelocityInjectionData& data)
{
    is.check("reading dClass");
    is >> data.dMain_;

    is.check("reading dProb");
    is >> data.dAux_;
    
    is.check("reading VxClass");
    is >> data.VxMain_;

    is.check("reading VxProb");
    is >> data.VxAux_;
    
    is.check("reading VyClass");
    is >> data.VyMain_;

    is.check("reading VyProb");
    is >> data.VyAux_;
    
    is.check("reading VzClass");
    is >> data.VzMain_;

    is.check("reading VzProb");
    is >> data.VzAux_;

    is.check("operator(Istream&, generalCoupledSizeVelocityInjectionData&)");

    return is;
}


// ************************************************************************* //
