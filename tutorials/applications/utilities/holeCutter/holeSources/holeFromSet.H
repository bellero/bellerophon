#ifndef HOLEFROMSET_H
#define HOLEFROMSET_H

#include "holeSource.H"
#include "triSurface.H"
#include "Switch.H"

namespace Foam
{

class holeFromSet : public holeSource
{
private:
    //- Name of the tri surface
    const word setName_;

public:
    //- Construct from mesh and dictionary
    holeFromSet(const fvMesh&, const dictionary&);

    TypeName("holeFromSet");
    //- Mark cells
    virtual label markCells(labelList&, const label mark) const;
};

}

#endif // HOLEFROMSET_H