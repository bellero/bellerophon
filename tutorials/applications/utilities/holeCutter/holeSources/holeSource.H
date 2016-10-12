#ifndef HOLESOURCE_H
#define HOLESOURCE_H

#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "word.H"

namespace Foam
{

class holeSource
{
protected:
    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Reference to the dictionary
    const dictionary& dict_;

public:
    TypeName("holeSource");

    declareRunTimeSelectionTable
    (
        autoPtr,
        holeSource,
        holeSource,
        (const fvMesh& mesh, const dictionary& dict),
        (mesh, dict)
    );

    //- Construct from mesh and dict
    holeSource(const fvMesh&, const dictionary&);

    //- Destructor
    virtual ~holeSource();

    //- Return selected turbulence model
    static autoPtr<holeSource> New(const fvMesh& mesh, const dictionary& dict);

    //- Return name of hole hole
    static word sourceType(const dictionary);

    //- Return hole map
    virtual autoPtr<labelList> holeMap() const;

    //- Return hole map
    virtual label markCells(labelList& map, const label mark = 1) const = 0;
};

}

#endif // HOLESOURCE_H