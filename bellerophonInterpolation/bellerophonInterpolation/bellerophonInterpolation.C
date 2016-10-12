/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bellerophonInterpolation.H"
#include "bellerophonInterface.H"
#include "bellerophonInterfaceField.H"
#include "bellerophonLduMatrix.H"
#include "FieldField.H"
#include "gradientSearch.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "faceSet.H"
#include "cellSet.H"
#include "OFstream.H"
#include "PstreamBuffers.H"
#include "Time.H"

#include "vector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellerophon, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bellerophon::bellerophon()
:
dictPtr_(),
meshPtr_(NULL),
patchPtrs_(0),
cifs_(0),
nAcceptors_(0),
acceptorCells_(0),
nDonors_(0),
primaryDonorCells_(0),
holeCells_(0),
newLiveCells_(0),
acceptorMapPtr_(),
deltas_(0),
interpolationMethod_(),
zoneSeeds_(),
zoneIDs_(0),
interfaceFaces_(0),
interfaceFlipMap_(0),
topologyUpToDate_(false),
interpolationUpToDate_(false),
donorItemsPtr_(),
continuityFields_(List<word>(0)),
forceZeroFields_(List<word>(0))
{
    if(debug)
    {
        Info<<"Creating bellerophon interpolation instance."<<nl<<endl;
    }
}

// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::bellerophon::~bellerophon()
{
    if(debug)
    {
        Info<<"Destructing Bellerophon interpolation instance."<<nl<<endl;
    }

    clearTopology();
}

// * * * * * * * * * * * * Private Member functions  * * * * * * * * * * * * //

void Foam::bellerophon::clearTopology() const
{
    // Clearing topology based information
    // If this is used for meshes with changes, maybe this can be
    // improved by only updating/resetting patches with change.


    if(debug)
    {
        Info<<"Clearing topological information of Bellerophon "
            <<"interpolation."<<nl<<endl;
    }

    clearInterpolation();

    topologyUpToDate_ = false;

    nAcceptors_= 0;
    acceptorCells_.setSize(0);
    primaryDonorCells_.setSize(0);
    holeCells_.setSize(0);
    deltas_.setSize(0);
    donorItemsPtr_.clear();
    interpolationMethod_.clear();

    zoneIDs_.setSize(0);
    zoneSeeds_.setSize(0);

}

void Foam::bellerophon::clearInterpolation() const
{
    if(debug)
    {
        Info<<"Clearing interpolation information of Bellerophon "
            <<"interpolation."<<nl<<endl;
    }

    interpolationUpToDate_ = false;

    nDonors_=0;
}

// * * * * * * * * * * * * Public Member functions   * * * * * * * * * * * * //

Foam::bellerophon& Foam::bellerophon::Interpolation()
{
    //- Create unique instance of bellerophon Interpolation
    static bellerophon bellerophonInterpolation;

    return bellerophonInterpolation;
}


Foam::label Foam::bellerophon::addInterface(Foam::bellerophonInterface& i)
{
    // If this is used for meshes with topological changes, maybe this can be
    // improved by only updating/resetting patches with change.

    if(debug>1)
    {
        Info<<"Adding interface to Bellerophon interpolation."<<nl
            <<"New size is "<< patchPtrs_.size()+1 << endl;
    }

    const bellerophonInterface* iPtr = &i;

    label index=-1;
    forAll(patchPtrs_,ptrI)
    {
        if(iPtr == patchPtrs_[ptrI])
        {
            WarningIn
            (
                "void Foam::bellerophon::addInterface"
                "("
                "    Foam::bellerophonInterface& i"
                ")"
            )<<"Trying to add patch which is already present in list."<<endl;
            index=ptrI;
        }
    }

    if(index<0)
    {
        index = patchPtrs_.size();
        patchPtrs_.setSize(index+1,const_cast<bellerophonInterface*>(iPtr));
        clearTopology();
    }

    cifs_.setSize(0);

    clearTopology();

    return index;
}

void Foam::bellerophon::removeInterface(Foam::bellerophonInterface& i)
{
    // If this is used for meshes with topological changes, maybe this can be
    // improved by only updating/resetting patches with change.

    const bellerophonInterface* iPtr = &i;

    bool found = false;
    forAll(patchPtrs_,ptrI)
    {
        if(found)
        {
            // Update position in list
            patchPtrs_[ptrI]->decreaseIndex();

            // Close hole from removing item
            patchPtrs_[ptrI-1]=patchPtrs_[ptrI];
        }
        else if(iPtr == patchPtrs_[ptrI])
        {
            found = true;
            if(debug>1)
            {
                Info<<"Removing interface "<<ptrI<<" from Bellerophon"
                    <<"interpolation."<<nl<<"New size is "
                    << patchPtrs_.size()-1 << endl;
            }
        }
    }

    if(found)
    {
        patchPtrs_.setSize(patchPtrs_.size()-1);
        clearTopology();
    }
    else
    {
        WarningIn
        (
            "void Foam::bellerophon::removeInterface"
            "("
            "    Foam::bellerophonInterface& i"
            ")"
        )<<"Trying to remove patch which is not present in list."<<endl;
    }

    clearTopology();

    cifs_.setSize(0);
}

const Foam::fvMesh& Foam::bellerophon::mesh() const
{
    if(!meshPtr_)
    {
        if(patchPtrs_.size() > 0)
        {
            meshPtr_ = &patchPtrs_[0]->mesh();
        }
        else
        {
            FatalError <<"No patches set"<<exit(FatalError);
        }
    }

    return *meshPtr_;
}

const Foam::dictionary& Foam::bellerophon::dict() const
{
    if(!dictPtr_.valid())
    {
        dictPtr_.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    "bellerophonDict",
                    mesh().time().constant(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
        dictPtr_().lookup("continuityFields") >> continuityFields_;
        dictPtr_().lookup("forceZeroFields") >> forceZeroFields_;
    }
    return dictPtr_();
}

void Foam::bellerophon::markHole() const
{
    // Resize and set to false
    newLiveCellsMarkup_.setSize(mesh().nCells());
    forAll(newLiveCellsMarkup_,i) newLiveCellsMarkup_[i] = false;

    forAll(acceptorCells_, i) newLiveCellsMarkup_[acceptorCells_[i]] = true;
    forAll(holeCells_, i) newLiveCellsMarkup_[holeCells_[i]] = true;
}

bool Foam::bellerophon::enforceContinuity(const word& name) const
{
    if(!dictPtr_.valid())
    {
        dict();
    }

    forAll(continuityFields_,fieldI)
    {
        if(continuityFields_[fieldI] == name)
        {
            if(debug)
            {
                Info<<"Continuity will be enforced on field "<<name<<endl;
            }
            return true;
        }
    }
    if(debug)
    {
        Info<<"Continuity will NOT be enforced on field "<<name<<endl;
    }
    return false;
}

bool Foam::bellerophon::forceZero(const word& name) const
{
    if(!dictPtr_.valid())
    {
        dict();
    }

    forAll(forceZeroFields_,fieldI)
    {
        if(forceZeroFields_[fieldI] == name)
        {
            if(debug)
            {
                Info<<"Zero will be enforced on field "<<name<<endl;
            }
            return true;
        }
    }
    if(debug)
    {
        Info<<"Zero will NOT be enforced on field "<<name<<endl;
    }
    return false;
}

template<class Type>
void Foam::bellerophon::interpolate
(
    const word& fieldName,
    UPstream::commsTypes commsType,
    bool bound
) const
{
    updateInterpolation();

    if
    (
        !mesh().thisDb().foundObject
                <GeometricField< Type, fvPatchField, volMesh > >(fieldName)
    )
    {
        if(debug) Info << "Field "<<fieldName<<" not found, not interpolating." <<endl;
        return;
    }

    if(debug)
    {
        Info<<"Interpolating "<<fieldName<<endl;
    }

    // Don't try this at home! De-consting the internal field to enforce
    // interpolation to acceptor cells. This might be doing stuff
    // The Wrong Way (TM), but that's how get access to the field values

    GeometricField<Type, fvPatchField, volMesh>& iField =
        const_cast< GeometricField<Type, fvPatchField, volMesh>& >
        (
            mesh().thisDb().lookupObject
                <GeometricField< Type, fvPatchField, volMesh > >(fieldName)
        );

    Field<Type> interpolatedPsi(nDonors_);

    const scalarList& pDW = primaryDonorWeights();

    forAll(primaryDonorCells_, donorI)
    {
        interpolatedPsi[donorI] =
            iField[primaryDonorCells_[donorI]] * pDW[donorI];
    }

    const List<interpolationItem>& oII = ownInterpolationItems();
    forAll(oII, itemI)
    {
        const interpolationItem& curItem = oII[itemI];
        interpolatedPsi[curItem.donorID()] +=
            iField[curItem.cellID()] * curItem.weight();
    }

    if(Pstream::parRun())
    {
        const scalar myProcNo = Pstream::myProcNo();

        autoPtr<PstreamBuffers> pBufsPtr
        (
            new PstreamBuffers(Pstream::nonBlocking)
        );

        forAll(neighbourInterpolationItems(), procI)
        {
            if(procI != myProcNo)
            {
                const List<interpolationItem>& neighbourItems
                    = neighbourInterpolationItems()[procI];

                UOPstream toDomain(procI, pBufsPtr());

                forAll(neighbourItems, itemI)
                {
                    const interpolationItem& curItem = neighbourItems[itemI];
                    toDomain << iField[curItem.cellID()] * curItem.weight();
                }
            }
        }

        pBufsPtr().finishedSends();

        forAll(neighbourValueToFieldMap(), procI)
        {
            if(procI != myProcNo)
            {
                const labelList& neighbourFieldMap =
                    neighbourValueToFieldMap()[procI];

                UIPstream str(procI, pBufsPtr());

                forAll(neighbourFieldMap, valueI)
                {
                    Type value;
                    str >> value;

                    interpolatedPsi[neighbourFieldMap[valueI]] += value;
                }
            }
        }

        if(bound)
        {
            forAll(interpolatedPsi, psiI)
            {
                Type zero = SMALL * pTraits<Type>::one;
                Type& curPsi = interpolatedPsi[psiI];
                curPsi = max(curPsi, zero);
            }
        }

        const labelListList& donorCols_ = donorCols();
        const labelListList& acceptorRows_ = acceptorRows();

        pBufsPtr.reset(new PstreamBuffers(Pstream::nonBlocking));

        forAll(donorCols_,procI)
        {
            if(procI != myProcNo)
            {
                UOPstream toDomain(procI, pBufsPtr());
                const labelList& curDonorCols = donorCols_[procI];
                const label* const __restrict__ curDonorColsPtr =
                    curDonorCols.begin();
                register const label nDonors = curDonorCols.size();
                for(register label donorI = 0; donorI < nDonors; donorI++)
                {
                    toDomain << interpolatedPsi[curDonorColsPtr[donorI]];
                }
            }
        }

        pBufsPtr().finishedSends();

        forAll(acceptorRows_,procI)
        {
            if(procI != myProcNo)
            {
                UIPstream str(procI, pBufsPtr());
                const labelList& curAcceptorRows = acceptorRows_[procI];
                const label* const __restrict__ curAcceptorRowPtr =
                    curAcceptorRows.begin();
                register const label nAcceptors = curAcceptorRows.size();
                for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
                {
                    register label row = curAcceptorRowPtr[acceptorI];
                    Type received;
                    str >> received;
                    iField[row] = received;
                }
            }
        }

        const labelList& curDonorCols = donorCols_[myProcNo];
        const labelList& curAcceptorRows = acceptorRows_[myProcNo];

        const label* const __restrict__ curDonorColsPtr =
            curDonorCols.begin();
        const label* const __restrict__ curAcceptorRowPtr =
            curAcceptorRows.begin();

        register const label nAcceptors = curAcceptorRows.size();

        for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
        {
            register const label row = curAcceptorRowPtr[acceptorI];
            iField[row] = interpolatedPsi[curDonorColsPtr[acceptorI]];

        }

    }
    else
    {
        if(bound)
        {
            forAll(acceptorCells_, acceptorI)
            {
                Type zero = SMALL * pTraits<Type>::one;
                iField[acceptorCells_[acceptorI]] =
                    max(zero, interpolatedPsi[acceptorI]);
            }
        }
        else
        {
            forAll(acceptorCells_, acceptorI)
            {
                iField[acceptorCells_[acceptorI]] = interpolatedPsi[acceptorI];
            }
        }
    }
}


// void Foam::bellerophon::updateMatrix
// (
//     Foam::bellerophonLduMatrix& matrix,
//     const Foam::bellerophonInterfaceField* interface
// ) const
// {
//     // Modification of the bellerophonLduMatrix matrix...
//     updateInterpolation();
//
//     //TODO:
//     // - store and reuse if mesh is static or at least for one time step!
//     // - Special treatment for serial runs - YES
//     // - move to update function as this independed from the matrix so solve for
//
//     // Get reference to lower and upper coeffs
//     scalar* __restrict__ lowerPtr = matrix.modLower().begin();
//     scalar* __restrict__ upperPtr = matrix.modUpper().begin();
//     const scalarField diag = matrix.diag();
//
//     const cellList& cells = mesh().cells();
//     const labelList& owner = mesh().faceOwner();
//     const label nIntFaces = mesh().nInternalFaces();
//
//     for ( label acceptorI = 0 ; acceptorI < nAcceptors_ ; acceptorI++ )
//     {
//         // Get faces of cell
//         const label curCell = acceptorCells_[acceptorI];
//
//     // First of all, we need to delete all off-diagonal coeffs in rows
//     // of the faceCells, since the value of these cells will be found
//     // from interpolation
//         labelList faces = cells[curCell];
//         const label nFaces=faces.size();
//         for ( label faceI = 0 ; faceI < nFaces ; faceI++ )
//         {
//             // check of cell is owner of face
//             const label curFace=faces[faceI];
//             if(curFace<nIntFaces)
//             {
//                 if(curCell == owner[curFace])
//                 {
//                     upperPtr[curFace]=0.0;
//                 }
//                 else
//                 {
//                     lowerPtr[curFace]=0.0;
//                 }
//             }
//         }
//     }
//
// }


void Foam::bellerophon::updateMatrix
(
    Foam::bellerophonLduMatrix& matrix,
    const Foam::bellerophonInterfaceField* interface
) const
{
    // Modification of the bellerophonLduMatrix matrix...
    updateInterpolation();

    // Get reference to lower and upper coeffs
    scalar* __restrict__ lowerPtr = matrix.modLower().begin();
    scalar* __restrict__ upperPtr = matrix.modUpper().begin();

    labelList map(mesh().nCells(), REGULAR_CELL);

    forAll(holeCells_, cellI) map[holeCells_[cellI]] = HOLE_CELL;
    forAll(acceptorCells_, cellI) map[acceptorCells_[cellI]] = ACCEPTOR_CELL;

    const labelList& own  = mesh().owner();
    const labelList& nei  = mesh().neighbour();

    forAll(own, faceI)
    {
//         if ( map[own[faceI]] == HOLE_CELL || map[own[faceI]] == ACCEPTOR_CELL )
        if ( map[own[faceI]] == ACCEPTOR_CELL )
        {
            upperPtr[faceI] = 0.0;
        }
//         if ( map[nei[faceI]] == HOLE_CELL || map[nei[faceI]] == ACCEPTOR_CELL )
        if ( map[nei[faceI]] == ACCEPTOR_CELL )
        {
            lowerPtr[faceI] = 0.0;
        }
    }
}

Foam::tmp<Foam::scalarField> Foam::bellerophon::correctSource
(
    const Foam::scalarField& oldSource,
    const Foam::scalarField& diag,
    const direction cmpt
) const
{
    tmp<scalarField> tSource(oldSource);
    scalarField& newSource = tSource();

    for(label acceptorI=0; acceptorI<nAcceptors_; acceptorI++)
    {
            newSource[acceptorCells_[acceptorI]] = 0.0;
    }

    return tSource;
}

template void Foam::bellerophon::interpolate<Foam::scalar>
(
    const Foam::word&,
    Foam::UPstream::commsTypes,
    bool
) const;

template void Foam::bellerophon::interpolate<Foam::vector>
(
    const Foam::word&,
    Foam::UPstream::commsTypes,
    bool
) const;

template void Foam::bellerophon::interpolate<Foam::tensor>
(
    const Foam::word&,
    Foam::UPstream::commsTypes,
    bool
) const;

template void Foam::bellerophon::interpolate<Foam::sphericalTensor>
(
    const Foam::word&,
    Foam::UPstream::commsTypes,
    bool
) const;

template void Foam::bellerophon::interpolate<Foam::symmTensor>
(
    const Foam::word&,
    Foam::UPstream::commsTypes,
    bool
) const;

