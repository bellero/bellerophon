#include "holeFromSet.H"

#include "addToRunTimeSelectionTable.H"
#include "cellSet.H"

namespace Foam
{
    defineTypeNameAndDebug(holeFromSet, 0);
    addToRunTimeSelectionTable(holeSource,holeFromSet, holeSource);
}


Foam::holeFromSet::holeFromSet
(
    const fvMesh& mesh,
    const dictionary& io
)
:
    holeSource(mesh, io),
    setName_(dict_.lookup("setName"))
{}

Foam::label
Foam::holeFromSet::markCells(labelList& map, const label mark) const
{
    cellSet cSet(mesh_, setName_,IOobject::MUST_READ,IOobject::NO_WRITE);

    Info<<"Read cellSet "<<setName_<<" with "<<cSet.size()<<" cells."<<endl;

    forAllConstIter(labelHashSet, cSet, iter)
    {
        map[iter.key()] = mark;
    }

    return cSet.size();
}
