/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef donorCellTetInterpolationMethod_H
#define donorCellTetInterpolationMethod_H

#include "bellerophonInterpolationMethod.H"

namespace Foam
{
class donorCellTetInterpolationMethod
:
public bellerophonInterpolationMethod
{
public:
    TypeName("donorCellTet");

    //- Construct from mesh, primary donor cells and deltas
    donorCellTetInterpolationMethod
    (
        const Foam::fvMesh& mesh,
        const Foam::labelList& primaryDonorCells,
        const Foam::labelList& acceptorCells,
        const Foam::vectorField& deltas,
        const Foam::autoPtr<Foam::List< Foam::searchItem> >& donorItemsPtr,
        const dictionary& dict
    );

    //- Destructor
    virtual ~donorCellTetInterpolationMethod();

    //- update interpolation coeffs
    virtual void update();
};

}
#endif // donorCellTetInterpolationMethod_H