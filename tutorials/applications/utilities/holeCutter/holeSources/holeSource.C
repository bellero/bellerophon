#include "holeSource.H"

namespace Foam
{
    defineTypeNameAndDebug(holeSource, 0);
    defineRunTimeSelectionTable(holeSource, holeSource);
}

Foam::holeSource::holeSource
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
mesh_(mesh),
dict_(dict.subDict(sourceType(dict)+"Coeffs"))
{}

Foam::holeSource::~holeSource()
{}

Foam::autoPtr<Foam::holeSource> Foam::holeSource::New
(
    const fvMesh& mesh,
    const dictionary& dict
){
    const word sourceT = sourceType(dict);

    Info<<"looking for type "<<sourceT<<endl;

    holeSourceConstructorTable::iterator constructorIter =
        holeSourceConstructorTablePtr_->find(sourceT);

    Info<<"was looking for type "<<sourceT<<endl;

    if (constructorIter == holeSourceConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "holeSource::New", dict
        )   << "Unknown hole source "
            << sourceType(dict) << nl << nl
            << "Valid hole sources are :" << endl
            << holeSourceConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<holeSource>
    (
        constructorIter()
        (
            mesh, dict
        )
    );
}

Foam::word Foam::holeSource::sourceType(const dictionary dict)
{
    return dict.lookup("holeSource");
}

Foam::autoPtr< Foam::labelList > Foam::holeSource::holeMap() const
{
    autoPtr<labelList> hmPtr(new labelList(mesh_.nCells()));
    labelList& hm = hmPtr();
    this->markCells(hm);
    return hmPtr;
}
