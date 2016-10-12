/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef donorCellInterpolationMethod_H
#define donorCellInterpolationMethod_H

#include "bellerophonInterpolationMethod.H"

namespace Foam
{
class donorCellInterpolationMethod
:
public bellerophonInterpolationMethod
{
public:
    TypeName("donorCell");

    //- Construct from mesh, primary donor cells and deltas
    donorCellInterpolationMethod
    (
        const Foam::fvMesh& mesh,
        const Foam::labelList& primaryDonorCells,
        const Foam::labelList& acceptorCells,
        const Foam::vectorField& deltas,
        const Foam::autoPtr<Foam::List< Foam::searchItem> >& donorItemsPtr,
        const dictionary& dict
    );

    //- Destructor
    virtual ~donorCellInterpolationMethod();

    //- update interpolation coeffs
    virtual void update();
};

}
#endif // donorCellInterpolationMethod_H