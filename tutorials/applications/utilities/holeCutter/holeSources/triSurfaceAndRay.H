#ifndef TRISURFACEANDRAY_H
#define TRISURFACEANDRAY_H

#include "holeSource.H"
#include "triSurface.H"
#include "Switch.H"

namespace Foam
{

class triSurfaceAndRay : public holeSource
{
private:
    //- Name of the tri surface
    const fileName surfaceName_;

    //- Tri surface
    triSurface surface_;

    //- Ray
    const vector ray_;

    //- Flip inside/outside
    const Switch flip_;

public:
    //- Construct from mesh and dictionary
    triSurfaceAndRay(const fvMesh&, const dictionary&);

    TypeName("triSurfaceAndRay");
    //- Mark cells
    virtual label markCells(labelList&, const label mark) const;
};

}

#endif // TRISURFACEANDRAY_H