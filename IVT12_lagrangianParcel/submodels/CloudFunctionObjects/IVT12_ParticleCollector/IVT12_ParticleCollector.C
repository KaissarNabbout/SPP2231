/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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

#include "IVT12_ParticleCollector.H"
#include "Pstream.H"
#include "surfaceWriter.H"
#include "randomGenerator.H"
#include "triangle.H"
#include "cloud.H"
#include "axesRotation.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::IVT12_ParticleCollector<CloudType>::makeLogFile
(
    const faceList& faces,
    const Field<point>& points,
    const Field<scalar>& area
)
{
    // Create the output file if not already created
    if (log_)
    {
        if (debug)
        {
            Info<< "Creating output file" << endl;
        }

        if (Pstream::master())
        {
            // Create directory if does not exist
            mkDir(this->writeTimeDir());

            // Open new file at start up
            outputFilePtr_.reset
            (
                new OFstream(this->writeTimeDir()/(type() + ".dat"))
            );

            outputFilePtr_()
                << "Face" << tab
                << "Time [s]" << tab
                << "Position [m]" << tab
                << "Size [m]" << tab
                << "Velocity [m/s]" << endl;
        }
    }
}


template<class CloudType>
void Foam::IVT12_ParticleCollector<CloudType>::initPolygons
(
    const List<Field<point>>& polygons
)
{
    mode_ = mtPolygon;

    label nPoints = 0;
    forAll(polygons, polyI)
    {
        label np = polygons[polyI].size();
        if (np < 3)
        {
            FatalIOErrorInFunction(this->coeffDict())
                << "polygons must consist of at least 3 points"
                << exit(FatalIOError);
        }

        nPoints += np;
    }

    label pointOffset = 0;
    points_.setSize(nPoints);
    faces_.setSize(polygons.size());
    area_.setSize(polygons.size());
    forAll(faces_, facei)
    {
        const Field<point>& polyPoints = polygons[facei];
        face f(identityMap(polyPoints.size()) + pointOffset);
        UIndirectList<point>(points_, f) = polyPoints;
        area_[facei] = f.mag(points_);
        faces_[facei].transfer(f);

        pointOffset += polyPoints.size();
    }
}


template<class CloudType>
void Foam::IVT12_ParticleCollector<CloudType>::initConcentricCircles()
{
    mode_ = mtConcentricCircle;

    vector origin(this->coeffDict().lookup("origin"));

    this->coeffDict().lookup("radius") >> radius_;
    nSector_ = this->coeffDict().template lookup<label>("nSector");

    label nS = nSector_;

    vector refDir;
    if (nSector_ > 1)
    {
        refDir = this->coeffDict().lookup("refDir");
        refDir -= normal_[0]*(normal_[0] & refDir);
        refDir /= mag(refDir);
    }
    else
    {
        // Set 4 quadrants for single sector cases
        nS = 4;
        refDir = normalised(perpendicular(normal_[0]));
    }

    scalar dTheta = degToRad(5.0);
    scalar dThetaSector = degToRad(360.0)/scalar(nS);
    label intervalPerSector = max(1, ceil(dThetaSector/dTheta));
    dTheta = dThetaSector/scalar(intervalPerSector);

    label nPointPerSector = intervalPerSector + 1;

    label nPointPerRadius = nS*(nPointPerSector - 1);
    label nPoint = radius_.size()*nPointPerRadius;
    label nFace = radius_.size()*nS;

    // Add origin
    nPoint++;

    points_.setSize(nPoint);
    faces_.setSize(nFace);
    area_.setSize(nFace);

    coordSys_ =
        coordinateSystems::cylindrical("coordSys", origin, normal_[0], refDir);

    List<label> ptIDs(identityMap(nPointPerRadius));

    points_[0] = origin;

    // Points
    forAll(radius_, radI)
    {
        label pointOffset = radI*nPointPerRadius + 1;

        for (label i = 0; i < nPointPerRadius; i++)
        {
            label pI = i + pointOffset;
            point pCyl(radius_[radI], i*dTheta, 0.0);
            points_[pI] = coordSys_.globalPosition(pCyl);
        }
    }

    // Faces
    DynamicList<label> facePts(2*nPointPerSector);
    forAll(radius_, radI)
    {
        if (radI == 0)
        {
            for (label secI = 0; secI < nS; secI++)
            {
                facePts.clear();

                // Append origin point
                facePts.append(0);

                for (label ptI = 0; ptI < nPointPerSector; ptI++)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = ptIDs.fcIndex(i - 1) + 1;
                    facePts.append(id);
                }

                label facei = secI + radI*nS;

                faces_[facei] = face(facePts);
                area_[facei] = faces_[facei].mag(points_);
            }
        }
        else
        {
            for (label secI = 0; secI < nS; secI++)
            {
                facePts.clear();

                label offset = (radI - 1)*nPointPerRadius + 1;

                for (label ptI = 0; ptI < nPointPerSector; ptI++)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = offset + ptIDs.fcIndex(i - 1);
                    facePts.append(id);
                }
                for (label ptI = nPointPerSector-1; ptI >= 0; ptI--)
                {
                    label i = ptI + secI*(nPointPerSector - 1);
                    label id = offset + nPointPerRadius + ptIDs.fcIndex(i - 1);
                    facePts.append(id);
                }

                label facei = secI + radI*nS;

                faces_[facei] = face(facePts);
                area_[facei] = faces_[facei].mag(points_);
            }
        }
    }
}


template<class CloudType>
void Foam::IVT12_ParticleCollector<CloudType>::collectParcelPolygon
(
    const point& p1,
    const point& p2
) const
{
    forAll(faces_, facei)
    {
        const label facePoint0 = faces_[facei][0];

        const point& pf = points_[facePoint0];

        const scalar d1 = normal_[facei] & (p1 - pf);
        const scalar d2 = normal_[facei] & (p2 - pf);

        if (sign(d1) == sign(d2))
        {
            // Did not cross polygon plane
            continue;
        }

        // Intersection point
        const point pIntersect = p1 + (d1/(d1 - d2))*(p2 - p1);

        // Identify if point is within the bounds of the face. Create triangles
        // between the intersection point and each edge of the face. If all the
        // triangle normals point in the same direction as the face normal, then
        // the particle is within the face. Note that testing for pointHits on
        // the face's decomposed triangles does not work due to ambiguity along
        // the diagonals.
        const face& f = faces_[facei];
        const vector a = f.area(points_);
        bool inside = true;
        for (label i = 0; i < f.size(); ++ i)
        {
            const label j = f.fcIndex(i);
            const triPointRef t(pIntersect, points_[f[i]], points_[f[j]]);
            if ((a & t.area()) < 0)
            {
                inside = false;
                break;
            }
        }

        // Add to the list of hits
        if (inside)
        {
            const label proci = Pstream::myProcNo();
            hitFaceIndices_.append(facei);
            hitFacePositions_[proci].append(pIntersect);
        }
    }
}


template<class CloudType>
void Foam::IVT12_ParticleCollector<CloudType>::collectParcelConcentricCircles
(
    const point& p1,
    const point& p2
) const
{
    label secI = -1;

    const scalar d1 = normal_[0] & (p1 - coordSys_.origin());
    const scalar d2 = normal_[0] & (p2 - coordSys_.origin());

    if (sign(d1) == sign(d2))
    {
        // Did not cross plane
        return;
    }
    
    // Intersection point
    const point pIntersect = p1 + (d1/(d1 - d2))*(p2 - p1);

    // Intersection point in cylindrical co-ordinate system
    const point pCyl = coordSys_.localPosition(pIntersect);

    scalar r = pCyl[0];

    if (r < radius_.last())
    {
        label radI = 0;
        while (r > radius_[radI])
        {
            radI++;
        }

        if (nSector_ == 1)
        {
            secI = 4*radI;
        }
        else
        {
            scalar theta = pCyl[1] + constant::mathematical::pi;

            secI =
                nSector_*radI
              + floor
                (
                    scalar(nSector_)*theta/constant::mathematical::twoPi
                );
        }
    }

    if (secI != -1)
    {
        hitFaceIndices_.append(secI);
        const label proci = Pstream::myProcNo();
        hitFacePositions_[proci].append(pIntersect);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::IVT12_ParticleCollector<CloudType>::write()
{
    Pstream::gatherList(hitFaceIDsList_);
    
    Pstream::gatherList(hitFaceTimes_);
    
    Pstream::gatherList(hitFacePositions_);
    
    Pstream::gatherList(hitFaceSizes_);
    
    Pstream::gatherList(hitFaceVelocities_);
    
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        labelList hitFaceIDs;
        scalarList hitParcelTimes;
        List<point> hitParcelPosition;
        scalarList hitParcelSizes;
        List<vector> hitParcelVelocity;
        
        hitFaceIDs = hitFaceIDsList_[proci];
        hitParcelTimes = hitFaceTimes_[proci];
        hitParcelPosition = hitFacePositions_[proci];
        hitParcelSizes = hitFaceSizes_[proci];
        hitParcelVelocity = hitFaceVelocities_[proci];
        
        forAll(hitFaceIDs, i)
        {
            if (outputFilePtr_.valid() && hitFaceIDs.size() > 0)
            {
                outputFilePtr_()
                    << hitFaceIDs[i] << tab 
                    << hitParcelTimes[i] << tab
                    << hitParcelPosition[i] << tab
                    << hitParcelSizes[i] << tab  
                    << hitParcelVelocity[i] << endl;
            }
        }
    }
    
    hitFaceIDsList_.clear();
    hitFaceTimes_.clear();
    hitFacePositions_.clear();
    hitFaceSizes_.clear();
    hitFaceVelocities_.clear(); 
               
    hitFaceIDsList_.setSize(Pstream::nProcs());
    hitFaceTimes_.setSize(Pstream::nProcs());
    hitFacePositions_.setSize(Pstream::nProcs());
    hitFaceSizes_.setSize(Pstream::nProcs());
    hitFaceVelocities_.setSize(Pstream::nProcs());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IVT12_ParticleCollector<CloudType>::IVT12_ParticleCollector
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    mode_(mtUnknown),
    parcelType_(this->coeffDict().lookupOrDefault("parcelType", -1)),
    removeCollected_(this->coeffDict().lookup("removeCollected")),
    points_(),
    faces_(),
    nSector_(0),
    radius_(),
    coordSys_("coordSys", vector::zero, axesRotation(sphericalTensor::I)),
    normal_(),
    surfaceFormat_(this->coeffDict().lookup("surfaceFormat")),
    resetOnWrite_(this->coeffDict().lookup("resetOnWrite")),
    totalTime_(0.0),
    mass_(),
    massTotal_(),
    massFlowRate_(),
    log_(this->coeffDict().lookup("log")),
    outputFilePtr_(),
    timeOld_(owner.mesh().time().value()),
    hitFaceIndices_(),
    hitFaceIDsList_(),
    hitFaceTimes_(),
    hitFacePositions_(),
    hitFaceSizes_(),
    hitFaceVelocities_()
{
    normal_ /= mag(normal_);

    word mode(this->coeffDict().lookup("mode"));
    if (mode == "polygon")
    {
        List<Field<point>> polygons(this->coeffDict().lookup("polygons"));

        initPolygons(polygons);

        vector n0(this->coeffDict().lookup("normal"));
        normal_ = vectorField(faces_.size(), n0);
    }
    else if (mode == "polygonWithNormal")
    {
        List<Tuple2<Field<point>, vector>> polygonAndNormal
        (
            this->coeffDict().lookup("polygons")
        );

        List<Field<point>> polygons(polygonAndNormal.size());
        normal_.setSize(polygonAndNormal.size());

        forAll(polygons, polyI)
        {
            polygons[polyI] = polygonAndNormal[polyI].first();
            normal_[polyI] = polygonAndNormal[polyI].second();
            normal_[polyI] /= mag(normal_[polyI]) + rootVSmall;
        }

        initPolygons(polygons);
    }
    else if (mode == "concentricCircle")
    {
        vector n0(this->coeffDict().lookup("normal"));
        normal_ = vectorField(1, n0);

        initConcentricCircles();
    }
    else
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "Unknown mode " << mode << ".  Available options are "
            << "polygon, polygonWithNormal and concentricCircle"
            << exit(FatalIOError);
    }

    mass_.setSize(faces_.size(), 0.0);
    massTotal_.setSize(faces_.size(), 0.0);
    massFlowRate_.setSize(faces_.size(), 0.0);
    
    hitFaceIDsList_.setSize(Pstream::nProcs());
    hitFaceTimes_.setSize(Pstream::nProcs());
    hitFacePositions_.setSize(Pstream::nProcs());
    hitFaceSizes_.setSize(Pstream::nProcs());
    hitFaceVelocities_.setSize(Pstream::nProcs());

    makeLogFile(faces_, points_, area_);
}


template<class CloudType>
Foam::IVT12_ParticleCollector<CloudType>::IVT12_ParticleCollector
(
    const IVT12_ParticleCollector<CloudType>& pc
)
:
    CloudFunctionObject<CloudType>(pc),
    mode_(pc.mode_),
    parcelType_(pc.parcelType_),
    removeCollected_(pc.removeCollected_),
    points_(pc.points_),
    faces_(pc.faces_),
    nSector_(pc.nSector_),
    radius_(pc.radius_),
    coordSys_(pc.coordSys_),
    normal_(pc.normal_),
    surfaceFormat_(pc.surfaceFormat_),
    resetOnWrite_(pc.resetOnWrite_),
    totalTime_(pc.totalTime_),
    mass_(pc.mass_),
    massTotal_(pc.massTotal_),
    massFlowRate_(pc.massFlowRate_),
    log_(pc.log_),
    outputFilePtr_(),
    timeOld_(0.0),
    hitFaceIndices_(),
    hitFaceIDsList_(),
    hitFaceTimes_(),
    hitFacePositions_(),
    hitFaceSizes_(),
    hitFaceVelocities_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::IVT12_ParticleCollector<CloudType>::~IVT12_ParticleCollector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::IVT12_ParticleCollector<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point& position0,
    bool& keepParticle
)
{
    if ((parcelType_ != -1) && (parcelType_ != p.typeId()))
    {
        return;
    }

    hitFaceIndices_.clear();

    switch (mode_)
    {
        case mtPolygon:
        {
            collectParcelPolygon
            (
                position0,
                p.position(this->owner().mesh())
            );
            break;
        }
        case mtConcentricCircle:
        {
            collectParcelConcentricCircles
            (
                position0,
                p.position(this->owner().mesh())
            );
            break;
        }
        default:
        {}
    }
    
    const label proci = Pstream::myProcNo();

    forAll(hitFaceIndices_, i)
    {
        label facei = hitFaceIndices_[i];
        
        hitFaceIDsList_[proci].append(facei);
        hitFaceTimes_[proci].append(this->owner().mesh().time().value());
        hitFaceSizes_[proci].append(p.d());
        hitFaceVelocities_[proci].append(p.U());

        if (removeCollected_)
        {
            keepParticle = false;
        }
    }
}


// ************************************************************************* //
