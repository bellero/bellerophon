/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef donorCellGradientInterpolationMethod_H
#define donorCellGradientInterpolationMethod_H

#include "bellerophonInterpolationMethod.H"

namespace Foam
{
class donorCellGradientInterpolationMethod
:
public bellerophonInterpolationMethod
{
public:
    TypeName("donorCellGradient");

    //- Construct from mesh, primary donor cells and deltas
    donorCellGradientInterpolationMethod
    (
        const Foam::fvMesh& mesh,
        const Foam::labelList& primaryDonorCells,
        const Foam::labelList& acceptorCells,
        const Foam::vectorField& deltas,
        const Foam::autoPtr<Foam::List< Foam::searchItem> >& donorItemsPtr,
        const dictionary& dict
    );

    //- Destructor
    virtual ~donorCellGradientInterpolationMethod();

    //- update interpolation coeffs
    virtual void update();
};

}
#endif // donorCellGradientInterpolationMethod_H