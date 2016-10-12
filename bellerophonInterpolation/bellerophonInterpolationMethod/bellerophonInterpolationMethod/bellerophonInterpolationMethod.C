/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bellerophonInterpolationMethod.H"

namespace Foam
{
    defineTypeNameAndDebug(bellerophonInterpolationMethod, 0);
    defineRunTimeSelectionTable(bellerophonInterpolationMethod, bellerophonInterpolationMethod);
}

Foam::bellerophonInterpolationMethod::bellerophonInterpolationMethod
(
    const Foam::fvMesh& mesh,
    const Foam::labelList& primaryDonorCells,
    const Foam::labelList& acceptorCells,
    const Foam::vectorField& deltas,
    const Foam::autoPtr<Foam::List< Foam::searchItem> >& donorItemsPtr,
    const dictionary& dict
)
:
mesh_(mesh),
primaryDonorCells_(primaryDonorCells),
acceptorCells_(acceptorCells),
deltas_(deltas),
donorItemsPtr_(donorItemsPtr),
dict_(dict.subDict(this->methodType(dict)+"Coeffs"))
{}

Foam::bellerophonInterpolationMethod::~bellerophonInterpolationMethod()
{}

Foam::autoPtr<Foam::bellerophonInterpolationMethod>
Foam::bellerophonInterpolationMethod::New
(
    const Foam::fvMesh& mesh,
    const Foam::labelList& primaryDonorCells,
    const Foam::labelList& acceptorCells,
    const Foam::vectorField& deltas,
    const Foam::autoPtr<Foam::List< Foam::searchItem> >& donorItemsPtr,
    const dictionary& dict
)
{
    const word sourceT = methodType(dict);

    bellerophonInterpolationMethodConstructorTable::iterator constructorIter =
    bellerophonInterpolationMethodConstructorTablePtr_->find(sourceT);

    if
    (
        constructorIter
        ==
        bellerophonInterpolationMethodConstructorTablePtr_->end()
    )
    {
        FatalIOErrorIn
        (
            "bellerophonInterpolationMethod::New", dict
        )   << "Unknown interpolation method "
            << methodType(dict) << nl << nl
            << "Valid methods are :" << endl
            << bellerophonInterpolationMethodConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<bellerophonInterpolationMethod>
    (
        constructorIter()
        (
            mesh, primaryDonorCells, acceptorCells, deltas, donorItemsPtr, dict
        )
    );
}

Foam::word Foam::bellerophonInterpolationMethod::methodType
(
    const dictionary& dict
)
{
    return dict.lookup("interpolationMethod");
}

void Foam::bellerophonInterpolationMethod::distribute
(
    const labelList& procToPatchMap
)
{
    List< List<interpolationItem> >& neighbourInterpolationItems =
        neighbourInterpolationItemsPtr_();

    const label nDonors = primaryDonorWeights().size();

    // Communicate interpolation items to neighbours
    if(Pstream::parRun())
    {
        const scalar myProcNo = Pstream::myProcNo();
        const scalar nProcs = Pstream::nProcs();

        PstreamBuffers pBufs(Pstream::nonBlocking);

        // Communicate values of non-local donors
        for (label procI = 0; procI < nProcs; procI++)
        {
            if (procI != myProcNo)
            {
                UOPstream toDomain(procI, pBufs);
                toDomain << neighbourInterpolationItems[procI];
            }
        }

        pBufs.finishedSends();

        // Collect the received values for acceptors
        for (label procI = 0; procI < nProcs; procI++)
        {
            if (procI != myProcNo)
            {
                UIPstream str(procI, pBufs);
                str >> neighbourInterpolationItems[procI];
            }
        }

        // Get face cell for the neighbourInterpolationItems
        forAll(neighbourInterpolationItems, procI)
        {
            if(neighbourInterpolationItems[procI].size() > 0)
            {
                List<interpolationItem>& neighbourItems =
                    neighbourInterpolationItems[procI];

                const label patchI = procToPatchMap[procI];

                if(patchI < 0)
                {
                    FatalErrorIn
                    (
                        "void Foam::bellerophon::updateInterpolation() "
                        "const"
                    )
                    << "Received interpolation items from proc " << procI
                    << "but could not find processor patch to this proc."
                    << abort(FatalError);
                }

                const labelList& faceCells =
                    mesh_.boundary()[patchI].faceCells();

                forAll(neighbourItems, itemI)
                {
                    // "Translate" face on the processor patch to cell next
                    // to the face
                    neighbourItems[itemI].cellID() =
                        faceCells[neighbourItems[itemI].cellID()];
                }
            }
        }
    }

    // Update matrix transfer data
    if(Pstream::parRun())
    {
        // Number of procs
        const label nProcs = Pstream::nProcs();
        const label myProcNo = Pstream::myProcNo();

        // Local search items
        const List<searchItem>& donorItems = donorItemsPtr_();

        // Number of local donors per proc
        labelList donorsPerProc(nProcs,0);

        forAll(donorItems, itemI)
        {
            donorsPerProc[donorItems[itemI].procID()]++;
        }

        {
            //- Local phi to proc addressing for Amul
            donorColPtr_.reset(new labelListList(nProcs));
            labelListList& donorCol = donorColPtr_();

            //- Value from proc to local phi addressing for Amul
            acceptorRowPtr_.reset(new labelListList(nProcs));
            labelListList& acceptorRow = acceptorRowPtr_();

            // Size lists for donor data
            forAll(donorsPerProc, procI)
            {
                const label nItems = donorsPerProc[procI];
                donorsPerProc[procI] = 0;
                donorCol[procI].setSize(nItems);
                acceptorRow[procI].setSize(nItems);
            }

            // Collect donor data
            forAll(donorItems, itemI)
            {
                const searchItem& curItem = donorItems[itemI];
                const label procI = curItem.procID();
                const label donorI = donorsPerProc[procI]++;
                donorCol[procI][donorI] = itemI ; // curItem.seed();
                acceptorRow[procI][donorI] = curItem.cellLabel();
            }

            PstreamBuffers pBufs(Pstream::nonBlocking);

            // Communicate rows to procs
            for (label procI = 0; procI < nProcs; procI++)
            {
                if (procI != myProcNo)
                {
                    UOPstream toDomain(procI, pBufs);
                    toDomain << acceptorRow[procI];
                }
            }

            pBufs.finishedSends();

            // Collect the received values for rows
            for (label procI = 0; procI < nProcs; procI++)
            {
                if (procI != myProcNo)
                {
                    UIPstream str(procI, pBufs);
                    str >> acceptorRow[procI];
                }
            }
        }
    }
    else
    {
        donorColPtr_.reset(new labelListList(1));
        //- Local phi to proc addressing for Amul
        labelList& donorCol = donorColPtr_()[0];

        acceptorRowPtr_.reset(new labelListList(1));
        //- Value from proc to local phi addressing for Amul
        labelList& acceptorRow = acceptorRowPtr_()[0];

        // Size lists for donor data
        donorCol.setSize(nDonors);
        acceptorRow.setSize(nDonors);

        // Collect donor data
        forAll(acceptorCells_, cellI)
        {
            donorCol[cellI] = cellI;
            acceptorRow[cellI] =
                acceptorCells_[cellI];
        }
    }
}
