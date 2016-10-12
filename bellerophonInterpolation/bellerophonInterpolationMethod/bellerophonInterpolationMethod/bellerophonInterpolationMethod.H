/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef bellerophonInterpolationMethod_H
#define bellerophonInterpolationMethod_H

#include "List.H"
#include "autoPtr.H"
#include "volFields.H"
#include "searchItem.H"
#include "interpolationItem.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
class bellerophonInterpolationMethod
{
protected:

        //- Mesh
        const fvMesh& mesh_;

        //- Primary donor cells
        const labelList& primaryDonorCells_;

        //- Acceptor cells
        const labelList& acceptorCells_;

        //- Vectors from primary donor cell centres to acceptor cell centres
        const vectorField& deltas_;

        //- Search items on this proc
        const autoPtr<List <searchItem> >& donorItemsPtr_;

        //- Dictionary
        const dictionary dict_;

        //- Interpolation weights for primary donor cells
        autoPtr<scalarField> primaryDonorWeightsPtr_;

        //- Interpolation items on own processor
        autoPtr< List<interpolationItem> > ownInterpolationItemsPtr_;

        //- Interpolation items for other processors
        autoPtr< List< List<interpolationItem> > >
            neighbourInterpolationItemsPtr_;

        // Maps interpolation contribution of interpolation items on
        // neighbour proc to interpolated field
        autoPtr< labelListList > neighbourValueToFieldMapPtr_;

        //- Index of the donor on the local domain stored separatly for
        //  each domain the acceptor cell is in
        //    Amul: Holds the index of interpolated Psi which needs to be send
        //          to the corresponding proc
        //          (transfer raw value, needs to be multiplied by diagnoal
        //          element for that row)
        //    Tmul: Holds the index of the row the received value from the
        //          corresponding proc needs to be added to in interpolated Psi
        //          (transfer value mutilplied with diagonal coeff,
        //          directly added to that row )
        autoPtr<labelListList> donorColPtr_;

        //- Rows of the acceptor cells on the local domain stored separatly for
        //  each domain the donor cell is in
        //    Amul: Holds the index of the row the received value from the
        //          corresponding proc needs to be added to
        //          (raw value transfered, needs to be multiplied by diagnoal
        //          element for that row)
        //    Tmul: Holds the index of Psi which needs to be send to the
        //          corresponding proc
        //          (transfer value mutilplied with diagonal coeff,
        //          directly added to that row )
        autoPtr<labelListList> acceptorRowPtr_;

    // Distribute interpolation information to neighbours and update matrix
    // coeffs
    void distribute(const labelList&);

public:
    TypeName("bellerophonInterpolationMethod");

    declareRunTimeSelectionTable
    (
        autoPtr,
        bellerophonInterpolationMethod,
        bellerophonInterpolationMethod,
        (
            const fvMesh& mesh,
            const labelList& primaryDonorCells,
            const labelList& acceptorCells,
            const vectorField& deltas,
            const autoPtr<List <searchItem> >& donorItemsPtr,
            const dictionary& dict
        ),
        (mesh, primaryDonorCells, acceptorCells, deltas, donorItemsPtr, dict)
    );

    //- Construct from mesh, primary donor cells and deltas_
    bellerophonInterpolationMethod
    (
        const fvMesh& mesh,
        const labelList& primaryDonorCells,
        const labelList& acceptorCells,
        const vectorField& deltas,
        const autoPtr<List <searchItem> >& donorItemsPtr,
        const dictionary& dict
    );

    //- Destructor
    virtual ~bellerophonInterpolationMethod();

    //- Return selected turbulence model
    static autoPtr<bellerophonInterpolationMethod> New
    (
        const fvMesh& mesh,
        const labelList& primaryDonorCells,
        const labelList& acceptorCells,
        const vectorField& deltas,
        const autoPtr<List <searchItem> >& donorItemsPtr,
        const dictionary& dict
    );

    //- Return name of interpolation method
    static word methodType(const dictionary&);

    //- Update interpolation coefficients
    virtual void update() = 0;

    // Access

        //- Access to primary donor weights
        inline const scalarField& primaryDonorWeights() const;

        //- Interpolation items in own domain
        inline const List<interpolationItem>& ownInterpolationItems() const;

        //- Interpolation items for neighbour domains
        inline const List< List<interpolationItem> >&
            neighbourInterpolationItems() const;

        //- Maps interpolation contribution of interpolation items on
        //  neighbour proc to interpolated field
        inline const labelListList& neighbourValueToFieldMap() const;

        //- Index of the donor on the local domain stored separatly for
        //  each domain the acceptor cell is in
        inline const labelListList& donorCols() const;

        //- Rows of the acceptor cells on the local domain stored separatly for
        //  each domain the donor cell is in
        inline const labelListList& acceptorRows() const;

};

} // End namespace Foam

#include "bellerophonInterpolationMethodI.H"

#endif // bellerophonInterpolationMethod_H