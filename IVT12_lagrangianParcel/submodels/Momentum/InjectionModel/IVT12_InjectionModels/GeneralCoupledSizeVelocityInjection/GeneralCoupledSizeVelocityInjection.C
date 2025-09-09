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

#include "GeneralCoupledSizeVelocityInjection.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::GeneralCoupledSizeVelocityInjection<CloudType>::GeneralCoupledSizeVelocityInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    positions_(this->coeffDict().lookup("position")),
    samplesFileName_(this->coeffDict().lookup("samplesFile")),
    samples_
    (
        IOobject
        (
            samplesFileName_,
            owner.db().time().constant() + "/Samples",
            owner.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
	),
    injectorCoordinates_(barycentric::uniform(NaN)),
    injectorCells_(positions_.size(), -1),
    injectorTetFaces_(positions_.size(), -1),
    injectorTetPts_(positions_.size(), -1),
    duration_(this->readDuration(dict, owner)),
    massFlowRate_(this->readMassFlowRate(dict, owner, duration_)),
    parcelsPerSecond_(this->readParcelsPerSecond(dict, owner)),
    autoAdjustPPS_
    (
        this->coeffDict().lookupOrDefault("autoAdjustPPS", false)
    ),
    rndInjection_
    (
        this->coeffDict().lookupOrDefault("randomInjection", false)
    ),
    rndInjectionBox_
    (
        this->coeffDict().lookupOrDefault("rndInjectionBox", vector::zero)
    ),
    ignoreOutOfBounds_
    (
        this->coeffDict().lookupOrDefault("ignoreOutOfBounds", false)
    ),
    rndGen_(this->owner().rndGen()),
    sizeDistributionCDF_
    (
    this->coeffDict().lookupOrDefault("dSampling", true)
    ),
    velocityDistributionCDF_
    (
    this->coeffDict().lookupOrDefault("vSampling", true)
    )
{    
    topoChange();
    
    checkSamples();
    
    if (!sizeDistributionCDF_ && velocityDistributionCDF_)
    {
    FatalErrorInFunction
        << type() << nl << "Combination not possible: "
        << "sizeDistributionCDF_ = " << sizeDistributionCDF_
        << " and velocityDistributionCDF_ = " << velocityDistributionCDF_
        << abort(FatalError);
    }
}

template<class CloudType>
Foam::GeneralCoupledSizeVelocityInjection<CloudType>::GeneralCoupledSizeVelocityInjection
(
    const GeneralCoupledSizeVelocityInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    //--------------------------
    positionsFile_(im.positionsFile_),
    positions_(im.positions_),
    samplesFileName_(im.samplesFileName_),
    samples_(im.samples_),
    injectorCoordinates_(im.injectorCoordinates_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    duration_(im.duration_),
    massFlowRate_(im.massFlowRate_, false),
    parcelsPerSecond_(im.parcelsPerSecond_),
    rndInjection_(im.rndInjection_),
    rndInjectionBox_(im.rndInjectionBox_),
    ignoreOutOfBounds_(im.ignoreOutOfBounds_),
    rndGen_(im.rndGen_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::GeneralCoupledSizeVelocityInjection<CloudType>::~GeneralCoupledSizeVelocityInjection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::GeneralCoupledSizeVelocityInjection<CloudType>::topoChange()
{
    forAll(positions_, i)
    {
//        injectorCells_[i] = this->owner().mesh().findNearestCell(positions_[i]);
        
        this->findCellAtPosition
        (
            positions_[i],
            injectorCoordinates_[i],
            injectorCells_[i],
            injectorTetFaces_[i],
            injectorTetPts_[i],
            !ignoreOutOfBounds_
        );
    }
}

template<class CloudType>
void Foam::GeneralCoupledSizeVelocityInjection<CloudType>::checkSamples()
{
    label nEntries = samples_.size();
    
    // Check if first value is negative. In case it is not, it is not needed to 
    // check the other values because the next checks make sure it is ascending
    if (samples_[0].dMain() < 0 || samples_[0].dAux() < 0 ||
        samples_[0].VxAux()[0] < 0 || samples_[0].VyAux()[0] < 0 ||
        samples_[0].VzAux()[0] < 0 )
    {
        FatalErrorInFunction
            << type() << nl << "Samples check: "
            << "Sizes and probabilities can't be negative."
            << abort(FatalError);
    }
        
    for (label i=1; i<nEntries; i++)
    {
        if ((samples_[i].dMain() <= samples_[i-1].dMain()) ||
            (samples_[i].dAux() <= samples_[i-1].dAux()))
        {
            FatalErrorInFunction
                << type() << nl <<  "Samples check (d): "
                << "Classes or probabilities are not in ascending order."
                << abort(FatalError);
        }
        
        if (samples_[i].VxMain().size() == samples_[i].VxAux().size())
        {
            for (label j=0; j<samples_[i].VxMain().size(); j++)
            {
                if ((j > 0 && samples_[i].VxMain()[j] < samples_[i].VxMain()[j-1]) ||
                    (j > 0 && samples_[i].VxAux()[j] < samples_[i].VxAux()[j-1]))
                {
                    FatalErrorInFunction
                        << type() << nl <<  "Samples check (Vx): "
                        << "Classes or probabilities are not in ascending order."
                        << abort(FatalError);
                }
            }
         }
         else
         {
            FatalErrorInFunction
                << type() << nl << "Samples check (Vx): "
                << "Lists do not have the same size."
                << abort(FatalError);
         }
        
        if (samples_[i].VyMain().size() == samples_[i].VyAux().size())
        {
            for (label j=0; j<samples_[i].VyMain().size(); j++)
            {
                if ((j > 0 && samples_[i].VyMain()[j] < samples_[i].VyMain()[j-1]) ||
                    (j > 0 && samples_[i].VyAux()[j] < samples_[i].VyAux()[j-1]))
                {
                    FatalErrorInFunction
                        << type() << nl << "Samples check (Vy): "
                        << "Classes or probabilities are not in ascending order."
                        << abort(FatalError);
                }
            }
         }
         else
         {
            FatalErrorInFunction
                << type() << nl << "Samples check (Vy): "
                << "Lists do not have the same size."
                << abort(FatalError);
         }
        
        if (samples_[i].VzMain().size() == samples_[i].VzAux().size())
        {
            for (label j=0; j<samples_[i].VzMain().size(); j++)
            {
                if ((j > 0 && samples_[i].VzMain()[j] < samples_[i].VzMain()[j-1]) ||
                    (j > 0 && samples_[i].VzAux()[j] < samples_[i].VzAux()[j-1]))
                {
                    FatalErrorInFunction
                        << type() << nl << "Samples check (Vz): "
                        << "Classes or probabilities are not in ascending order."
                        << abort(FatalError);
                }
            }
         }
         else
         {
            FatalErrorInFunction
                << type() << nl << "Samples check (Vz): "
                << "Lists do not have the same size."
                << abort(FatalError);
         }
     }
}

template<class CloudType>
Foam::scalar Foam::GeneralCoupledSizeVelocityInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::GeneralCoupledSizeVelocityInjection<CloudType>::nParcelsToInject
(
    const scalar time0,
    const scalar time1
)
{    
    if (time0 >= 0 && time0 < duration_)
    {
        if (this->uniformParcelSize_ == Foam::injectionModel::uniformParcelSize::nParticle 
            && autoAdjustPPS_)
        {
            // Add the influence of nP in the number of parcels to be injected.
            // For instance, if the PDA measured 1000 droplets per second, it
            // means that parcelsPerSecond_ = 1000 always (read in the dict). 
            // Therefore, it is necessary to divide it by nP to have the real
            // number of parcels to be injected.
            
            return 
                floor 
                (
                    parcelsPerSecond_->integral(0, time1) / this->nParticleFixed_ 
                    - this->parcelsAddedTotal()
                );
        }
        else
        {
            return 
                floor 
                (
                    parcelsPerSecond_->integral(0, time1) - this->parcelsAddedTotal()
                );
        }
    }
    else
    {
        return 0;
    }
}

template<class CloudType>
Foam::scalar Foam::GeneralCoupledSizeVelocityInjection<CloudType>::massToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (time0 >= 0 && time0 < duration_)
    {
        return massFlowRate_->integral(time0, time1);
    }
    else
    {
        return 0;
    }
}

template<class CloudType>
void Foam::GeneralCoupledSizeVelocityInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    barycentric& coordinates,
    label& celli,
    label& tetFacei,
    label& tetPti,
    label& facei
)
{
    // Randomly choose a position to inject within the limits set by rndInjectionBox
    if (rndInjection_)
    {
        vector positionDisplacement_ = Zero;
        
        if (magSqr(rndInjectionBox_) > 0) // check if the vector is not nul
        {  
            // Generate random displacement on master processor
            if (Pstream::master())
            {
                positionDisplacement_.x() = rndInjectionBox_.x()*rndGen_.scalarAB(-1,1);
                positionDisplacement_.y() = rndInjectionBox_.y()*rndGen_.scalarAB(-1,1);
                positionDisplacement_.z() = rndInjectionBox_.z()*rndGen_.scalarAB(-1,1);
            }
        }
        
        // Broadcast the random displacement to all processors to ensure that they
        // use the same value later to search for the cell
        Pstream::scatter(positionDisplacement_.x());
        Pstream::scatter(positionDisplacement_.y());
        Pstream::scatter(positionDisplacement_.z());
        
        // Because there is only one injection position, positions_[0] works here.
        // In the case of more than one injection position, the randomise_ concept from
        // MomentumLookupTableInjection can be applied.
        vector position = positions_[0] + positionDisplacement_;
        
        this->findCellAtPosition
        (
            position,
            coordinates,
            celli,
            tetFacei,
            tetPti,
            !ignoreOutOfBounds_
        );
    }
    else
    {
        coordinates = injectorCoordinates_[0];
        celli = injectorCells_[0];
        tetFacei = injectorTetFaces_[0];
        tetPti = injectorTetPts_[0];
    }
}

template<class CloudType>
void Foam::GeneralCoupledSizeVelocityInjection<CloudType>::setProperties
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    typename CloudType::parcelType::trackingData& td,
    typename CloudType::parcelType& parcel
)
{
    //- Get the data from injector samples
    label classID = 0;
    
    if (sizeDistributionCDF_)
    {
        classID = getClassID(classID);
    }
	parcel.d() = sampleSize(classID);
	parcel.U() = sampleVelocity(classID);
}

template<class CloudType>
Foam::label Foam::GeneralCoupledSizeVelocityInjection<CloudType>::getClassID
(
    label classID
)
{
    // Find which diameter's class contains the probability rndProb    
    scalar rndProb = rndGen_.scalar01();
    
    label classIDmax = samples_.size() - 1;

    while (rndProb > samples_[classID].dAux() && classID < classIDmax) 
    {
      classID++;
    }
    
    return classID;
}

template<class CloudType>
Foam::scalar Foam::GeneralCoupledSizeVelocityInjection<CloudType>::sampleSize
(
    const label classID
)
{
    distributions::standardNormal& stdNormal = this->owner().stdNormal();
    
    scalar dInject;
    
    // Calculate randomly the diameter to inject within the class
    if (sizeDistributionCDF_)
    {
        scalar dMain = samples_[classID].dMain();
        
        scalar dMainSize = samples_[classID].dMain() - samples_[classID-1].dMain();
        
        scalar dRandom = dMainSize*rndGen_.scalar01();
        
        dInject = dMain - dRandom;
    }
    // Calculate the diameter to inject based on a normal distribution
    else
    {
        dInject = -1.0;
        while (dInject < 0.0) // avoid negative diameters
        {
            scalar dMain = samples_[classID].dMain();
            
            scalar dRandom = samples_[classID].dAux()*stdNormal.sample();
            
            dInject = dMain + dRandom;
        }
    }

    return dInject;
}

template<class CloudType>
Foam::vector Foam::GeneralCoupledSizeVelocityInjection<CloudType>::sampleVelocity
(
    const label classID
)
{
    distributions::standardNormal& stdNormal = this->owner().stdNormal();
    
    vector V;
    label VxID = 0;
    label VyID = 0;
    label VzID = 0;
    
    if (sizeDistributionCDF_ && velocityDistributionCDF_)
    {
        // Find which Vx's class contains the probability rndProb
        scalar rndProb = rndGen_.scalar01();
        
        label VxIDmax = samples_[classID].VxAux().size() - 1;

        while (rndProb > samples_[classID].VxAux()[VxID] && VxID < VxIDmax) 
        {
          VxID++;
        }
        
        // Find which Vy's class contains the probability rndProb
        rndProb = rndGen_.scalar01();
        
        label VyIDmax = samples_[classID].VyAux().size() - 1;

        while (rndProb > samples_[classID].VyAux()[VyID] && VyID < VyIDmax) 
        {
          VyID++;
        }
        
        // Find which Vz's class contains the probability rndProb
        rndProb = rndGen_.scalar01();
        
        label VzIDmax = samples_[classID].VzAux().size() - 1;

        while (rndProb > samples_[classID].VzAux()[VzID] && VzID < VzIDmax) 
        {
          VzID++;
        }
        
        // Calculate randomly the Vx to inject within the class
        scalar VxClass = samples_[classID].VxMain()[VxID];
        scalar VxClassSize = samples_[classID].VxMain()[VxID] - samples_[classID].VxMain()[VxID-1];
        scalar VxRandom = VxClassSize*rndGen_.scalar01();
        V.x() = VxClass - VxRandom;
        
        // Calculate randomly the Vy to inject within the class
        scalar VyClass = samples_[classID].VyMain()[VyID];
        scalar VyClassSize = samples_[classID].VyMain()[VyID] - samples_[classID].VyMain()[VyID-1];
        scalar VyRandom = VyClassSize*rndGen_.scalar01();
        V.y() = VyClass - VyRandom;
        
        // Calculate randomly the Vz to inject within the class
        scalar VzClass = samples_[classID].VzMain()[VzID];
        scalar VzClassSize = samples_[classID].VzMain()[VzID] - samples_[classID].VzMain()[VzID-1];
        scalar VzRandom = VzClassSize*rndGen_.scalar01();
        V.z() = VzClass - VzRandom;
    }
    else
    {
        // Calculate Vx to inject based on a normal distribution
        scalar VxMain = samples_[classID].VxMain()[VxID];
        scalar VxRandom = samples_[classID].VxAux()[VxID]*stdNormal.sample();
        V.x() = VxMain + VxRandom;
        
        // Calculate Vy to inject based on a normal distribution
        scalar VyMain = samples_[classID].VyMain()[VyID];
        scalar VyRandom = samples_[classID].VyAux()[VyID]*stdNormal.sample();
        V.y() = VyMain + VyRandom;
        
        // Calculate Vz to inject based on a normal distribution
        scalar VzMain = samples_[classID].VzMain()[VzID];
        scalar VzRandom = samples_[classID].VzAux()[VzID]*stdNormal.sample();
        V.z() = VzMain + VzRandom;
    }
    
    return V;
}

template<class CloudType>
bool Foam::GeneralCoupledSizeVelocityInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::GeneralCoupledSizeVelocityInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
