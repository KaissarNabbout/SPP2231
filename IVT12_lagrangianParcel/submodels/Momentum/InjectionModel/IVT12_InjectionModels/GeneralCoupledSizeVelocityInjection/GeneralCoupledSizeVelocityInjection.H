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
    Foam::GeneralCoupledSizeVelocityInjection

Description
Generalist injection model to introduce parcels into the domain based on:
- normal distribution for size and velocity (not correlated)
- CDF for size and normal distribution for velocity (correlated)
- CDF for size and velocity (correlated)

Correlated means that the provided velocity distribution is correlated to the
the data from the randomly chosen class for the droplet size.

SourceFiles
    GeneralCoupledSizeVelocityInjection.C

\*---------------------------------------------------------------------------*/

#ifndef GeneralCoupledSizeVelocityInjection_H
#define GeneralCoupledSizeVelocityInjection_H

#include "InjectionModel.H"
#include "generalCoupledSizeVelocityInjectionDataIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class GeneralCoupledSizeVelocityInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class GeneralCoupledSizeVelocityInjection
:
    public InjectionModel<CloudType>
{
    // Private data

        //- Name of file containing positions data
        const word positionsFile_;

        //- Field of parcel positions
/*        GlobalIOField<vector> positions_;*/
        List<point> positions_;

        //- Name of file containing samples data
        const word samplesFileName_;

        //- List of samples for the injector
        generalCoupledSizeVelocityInjectionDataIOList samples_;

        //- List of coordinates corresponding to injector positions
        List<barycentric> injectorCoordinates_;

        //- List of cell labels corresponding to injector positions
        labelList injectorCells_;

        //- List of tetFace labels corresponding to injector positions
        labelList injectorTetFaces_;

        //- List of tetPt labels corresponding to injector positions
        labelList injectorTetPts_;

        //- Injection duration [s]
        scalar duration_;

        //- Mass flow rate relative to SOI []
        const autoPtr<Function1<scalar>> massFlowRate_;

        //- Number of parcels to introduce per second
        const autoPtr<Function1<scalar>> parcelsPerSecond_;
        
        //- Flag to indicate if the PPS is to be automatically adjusted based on nP
        bool autoAdjustPPS_;
        
        //- Flag to indicate to randomise injection positions
        bool rndInjection_;
        
        //- Limits (dx, dy, dz) where the parcels can be randomly injected
        //  relative to the injection position (x, y, z) [m]
        const vector rndInjectionBox_;
        
        //- Flag to suppress errors if particle injection site is out-of-bounds
        Switch ignoreOutOfBounds_;
        
        //- Random number generator - used by some injection routines
        randomGenerator& rndGen_;

        //- Parcel size distribution model (to be properly replaced in the future)
        //  true = CDF ; false = MeanStd
        bool sizeDistributionCDF_;

        //- Parcel velocity distribution model (to be properly replaced in the future)
        //  true = CDF ; false = MeanStd
        bool velocityDistributionCDF_;
        
public:

    //- Runtime type information
    TypeName("GeneralCoupledSizeVelocityInjection");


    // Constructors

    //- Construct from dictionary
    GeneralCoupledSizeVelocityInjection
    (
        const dictionary& dict,
        CloudType& owner,
        const word& modelName
    );

    //- Construct copy
    GeneralCoupledSizeVelocityInjection(const GeneralCoupledSizeVelocityInjection<CloudType>& im);

    //- Construct and return a clone
    virtual autoPtr<InjectionModel<CloudType> > clone() const
    {
        return autoPtr<InjectionModel<CloudType> >
        (
            new GeneralCoupledSizeVelocityInjection<CloudType>(*this)
        );
    }

    //- Destructor
    virtual ~GeneralCoupledSizeVelocityInjection();

    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void topoChange();
        
        //- Check if the samples are correctly set
        virtual void checkSamples();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label nParcelsToInject(const scalar time0, const scalar time1);

        //- Parcel mass to introduce relative to SOI
        virtual scalar massToInject(const scalar time0, const scalar time1);

    // Injection geometry

        //- Set the injection position and owner cell, tetFace and tetPt
        virtual void setPositionAndCell
        (
            const label parcelI,
            const label nParcels,
            const scalar time,
            barycentric& coordinates,
            label& celli,
            label& tetFacei,
            label& tetPti,
            label& facei
        );

        //- Set parcel properties (size, velocity, etc.)
        virtual void setProperties
        (
            const label parcelI,
            const label nParcels,
            const scalar time,
            typename CloudType::parcelType::trackingData& td,
            typename CloudType::parcelType& parcel
        );
        
        //- Get ID of the size class to inject parcel
        virtual label getClassID(label classID);
        
        //- Select size to inject based on classID
        virtual scalar sampleSize(const label classID);
        
        //- Select velocity of the parcel based on classID
        virtual vector sampleVelocity(const label classID);

        //- Flag to identify whether model fully describes the parcel
        virtual bool fullyDescribed() const;

        //- Return flag to identify whether or not injection of parcelI is
        //  permitted
        virtual bool validInjection(const label parcelI);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GeneralCoupledSizeVelocityInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
