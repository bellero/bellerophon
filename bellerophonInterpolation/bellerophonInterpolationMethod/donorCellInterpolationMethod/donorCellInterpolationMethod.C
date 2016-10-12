/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "donorCellInterpolationMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "searchItem.H"
#include "PstreamBuffers.H"

namespace Foam
{
    defineTypeNameAndDebug(donorCellInterpolationMethod, 0);
    addToRunTimeSelectionTable(bellerophonInterpolationMethod,donorCellInterpolationMethod, bellerophonInterpolationMethod);
}

// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

Foam::donorCellInterpolationMethod::donorCellInterpolationMethod
(
    const Foam::fvMesh& mesh,
    const Foam::labelList& primaryDonorCells,
    const Foam::labelList& acceptorCells,
    const Foam::vectorField& deltas,
    const Foam::autoPtr<Foam::List< Foam::searchItem> >& donorItemsPtr,
    const dictionary& dict
)
:
bellerophonInterpolationMethod
(
    mesh,
    primaryDonorCells,
    acceptorCells,
    deltas,
    donorItemsPtr,
    dict
)
{}

// * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * //

Foam::donorCellInterpolationMethod::
~donorCellInterpolationMethod()
{}

// * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * * //

void Foam::donorCellInterpolationMethod::update()
{
    // Number of local donor cells
    const label nDonors = primaryDonorCells_.size();

    // Weights of primary donors cells are 1.0
    primaryDonorWeightsPtr_.reset(new scalarField(nDonors,1.0));

    // Number of procs
    const label nProcs = Pstream::nProcs();

    // No other interpolation donors, set sizes to zero
    ownInterpolationItemsPtr_.reset
    (
        new List<interpolationItem>(0)
    );

    neighbourInterpolationItemsPtr_.reset
    (
        new List< List<interpolationItem> >(nProcs, List<interpolationItem>(0))
    );

    neighbourValueToFieldMapPtr_.reset
    (
        new labelListList(nProcs,labelList(0))
    );

//     // Update matrix transfer data
//     if(Pstream::parRun())
//     {
//         // Number of procs
// 
//         // Local search items
//         const List<searchItem>& donorItems = donorItemsPtr_();
// 
//         // Number of local donors per proc
//         labelList donorsPerProc(nProcs,0);
// 
//         forAll(donorItems, itemI)
//         {
//             donorsPerProc[donorItems[itemI].procID()]++;
//         }
// 
//         {
//             //- Local phi to proc addressing for Amul
//             donorColPtr_.reset(new labelListList(nProcs));
//             labelListList& donorCol = donorColPtr_();
// 
//             //- Value from proc to local phi addressing for Amul
//             acceptorRowPtr_.reset(new labelListList(nProcs));
//             labelListList& acceptorRow = acceptorRowPtr_();
// 
//             // Size lists for donor data
//             forAll(donorsPerProc, procI)
//             {
//                 const label nItems = donorsPerProc[procI];
//                 donorsPerProc[procI] = 0;
//                 donorCol[procI].setSize(nItems);
//                 acceptorRow[procI].setSize(nItems);
//             }
// 
//             // Collect donor data
//             forAll(donorItems, itemI)
//             {
//                 const searchItem& curItem = donorItems[itemI];
//                 const label procI = curItem.procID();
//                 const label donorI = donorsPerProc[procI]++;
//                 donorCol[procI][donorI] = itemI ; // curItem.seed();
//                 acceptorRow[procI][donorI] = curItem.cellLabel();
//             }
// 
//             PstreamBuffers pBufs(Pstream::nonBlocking);
// 
//             // Communicate rows to procs
//             for (label procI = 0; procI < nProcs; procI++)
//             {
//                 if (procI != myProcNo)
//                 {
//                     UOPstream toDomain(procI, pBufs);
//                     toDomain << acceptorRow[procI];
//                 }
//             }
// 
//             pBufs.finishedSends();
// 
//             // Collect the received values for rows
//             for (label procI = 0; procI < nProcs; procI++)
//             {
//                 if (procI != myProcNo)
//                 {
//                     UIPstream str(procI, pBufs);
//                     str >> acceptorRow[procI];
//                 }
//             }
//         }
//     }
//     else
//     {
//         donorColPtr_.reset(new labelListList(1));
//         //- Local phi to proc addressing for Amul
//         labelList& donorCol = donorColPtr_()[0];
// 
//         acceptorRowPtr_.reset(new labelListList(1));
//         //- Value from proc to local phi addressing for Amul
//         labelList& acceptorRow = acceptorRowPtr_()[0];
// 
//         // Size lists for donor data
//         donorCol.setSize(nDonors);
//         acceptorRow.setSize(nDonors);
// 
//         // Collect donor data
//         forAll(acceptorCells_, cellI)
//         {
//             donorCol[cellI] = cellI;
//             acceptorRow[cellI] =
//                 acceptorCells_[cellI];
//         }
//     }

    distribute(labelList(Pstream::nProcs(),-1));
}

