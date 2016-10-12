#include "triSurfaceAndRay.H"

#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

namespace Foam
{
    defineTypeNameAndDebug(triSurfaceAndRay, 0);
    addToRunTimeSelectionTable(holeSource,triSurfaceAndRay, holeSource);
}


Foam::triSurfaceAndRay::triSurfaceAndRay
(
    const fvMesh& mesh,
    const dictionary& io
)
:
    holeSource(mesh, io),
    surfaceName_(dict_.lookup("surfaceName")),
    surface_(surfaceName_),
    ray_(dict_.lookup("ray")),
    flip_(dict_.lookupOrDefault("flip", false))
{}

Foam::label
Foam::triSurfaceAndRay::markCells(labelList& map, const label mark) const
{
    const List<vector>& ccs = mesh_.cellCentres();

    if(map.size() != ccs.size())
    {
        FatalErrorIn
        (
            "triSurfaceAndRay::markCells"
            "(labelList& map, const label mark)"
        )
        << "Size of map and mesh doesn't match."
        << exit(FatalError);
    }

    //- Points of the surface
    const Field<point>& points = surface_.localPoints();

    //- Faces of the surface
    const List<labelledTri>& faces = surface_.localFaces();

    const List<point>& fc = surface_.faceCentres();

    //- Faces of the surface
    const List<vector>& normals = surface_.faceNormals();

    //- Bounding box
    boundBox bound(points);

    Info<<"Bounding box of triSurface: "<<bound<<endl;
    Info<<"Flip: "<<flip_<<endl;

    //- Nummber of marked cells
    label nMarked = 0;

    forAll(ccs, cellI)
    {
        if(bound.contains(ccs[cellI]))
        {
            const point cc = ccs[cellI];
//             Info<<bound<<" contains "<<cc<<endl;
            label intersections = 0;
            forAll(fc,faceI)
            {
                pointHit hit = faces[faceI].ray(cc,ray_,points, intersection::HALF_RAY);
                if(hit.hit())
                {
                    intersections++;
                }
            }
            if( (intersections%2 == 0) == flip_)
            {
                map[cellI] = mark;
                nMarked++;
            }
//             Info<<"There were "<<intersections
//                 <<" faces hit by the ray starting at centre of cell "<<cellI
//                 <<endl;
        }
        else if (flip_)
        {
            map[cellI] = mark;
            nMarked++;
        }
    }

    return nMarked;
}
