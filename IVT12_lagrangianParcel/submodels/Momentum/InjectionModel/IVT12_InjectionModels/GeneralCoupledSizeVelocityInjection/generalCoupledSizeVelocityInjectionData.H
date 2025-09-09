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
    Foam::generalCoupledSizeVelocityInjectionData

Description
    Container class to provide injection data for momentum parcels

SourceFiles
    generalCoupledSizeVelocityInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef generalCoupledSizeVelocityInjectionData_H
#define generalCoupledSizeVelocityInjectionData_H

#include "dictionary.H"
#include "vector.H"
#include "point.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class generalCoupledSizeVelocityInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const generalCoupledSizeVelocityInjectionData&
);

Istream& operator>>
(
    Istream&,
    generalCoupledSizeVelocityInjectionData&
);

/*---------------------------------------------------------------------------*\
               Class generalCoupledSizeVelocityInjectionData Declaration
\*---------------------------------------------------------------------------*/

class generalCoupledSizeVelocityInjectionData
{
protected:

    // Parcel properties

    //- Class center for diameter [m]
    scalar dMain_;

    //- Cumulative probability for diameter []
    scalar dAux_;
    
    //- Class center for Vx [m/s]
    scalarList VxMain_;
    
    //- Cumulative probability for Vx []
    scalarList VxAux_;
    
    //- Class center for Vy [m/s]
    scalarList VyMain_;
    
    //- Cumulative probability for Vy []
    scalarList VyAux_;
    
    //- Class center for Vz [m/s]
    scalarList VzMain_;
    
    //- Cumulative probability for Vz []
    scalarList VzAux_;	
		
public:

    //- Runtime type information
    TypeName("generalCoupledSizeVelocityInjectionData");

    // Constructors

    //- Null constructor
    generalCoupledSizeVelocityInjectionData();

    //- Construct from dictionary
    generalCoupledSizeVelocityInjectionData(const dictionary& dict);

    //- Construct from Istream
    generalCoupledSizeVelocityInjectionData(Istream& is);


    //-Destructor
    virtual ~generalCoupledSizeVelocityInjectionData();


    // Access

    //- Return const access to the class center for diameter
    inline scalar dMain() const;

    //- Return const access to the cumulative probability for diameter
    inline scalar dAux() const;

    //- Return const access to the class center for Vx
    inline scalarList VxMain() const;

    //- Return const access to the cumulative probability for Vx
    inline scalarList VxAux() const;

    //- Return const access to the class center for Vy
    inline scalarList VyMain() const;

    //- Return const access to the cumulative probability for Vy
    inline scalarList VyAux() const;

    //- Return const access to the class center for Vz
    inline scalarList VzMain() const;

    //- Return const access to the cumulative probability for Vz
    inline scalarList VzAux() const;
    
    // Edit

    //- Return const access to the class center for diameter
    inline scalar& dMain();

    //- Return const access to the cumulative probability for diameter
    inline scalar& dAux();

    //- Return const access to the class center for Vx
    inline scalarList& VxMain();

    //- Return const access to the cumulative probability for Vx
    inline scalarList& VxAux();

    //- Return const access to the class center for Vy
    inline scalarList& VyMain();

    //- Return const access to the cumulative probability for Vy
    inline scalarList& VyAux();

    //- Return const access to the class center for Vz
    inline scalarList& VzMain();

    //- Return const access to the cumulative probability for Vz
    inline scalarList& VzAux();

    // I-O

    // Friend Operators

        friend bool operator==
        (
            const generalCoupledSizeVelocityInjectionData& a,
            const generalCoupledSizeVelocityInjectionData& b
        )
        {
            notImplemented
            (
                "operator=="
                "("
                    "const generalCoupledSizeVelocityInjectionData&, "
                    "const generalCoupledSizeVelocityInjectionData&"
                ")"
            );
            return false;
        }

        friend bool operator!=
        (
            const generalCoupledSizeVelocityInjectionData& a,
            const generalCoupledSizeVelocityInjectionData& b
        )
        {
            notImplemented
            (
                "operator=="
                "("
                    "const generalCoupledSizeVelocityInjectionData&, "
                    "const generalCoupledSizeVelocityInjectionData&"
                ")"
            );
            return false;
        }

        //- Ostream operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const generalCoupledSizeVelocityInjectionData& data
        );

        //- Istream operator
        friend Istream& operator>>
        (
            Istream& is,
            generalCoupledSizeVelocityInjectionData& data
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "generalCoupledSizeVelocityInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
