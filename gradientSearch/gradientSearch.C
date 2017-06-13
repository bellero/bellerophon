/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "gradientSearch.H"
#include "GeometricField.H"
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "ops.H"
#include "SLList.H"
#include "PstreamBuffers.H"
#include "polyMesh.H"
#include "meshSearch.H"

// #include <fstream>
// #include <unistd.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gradientSearch, 0);
}

// * * * * * * * * * * * * * * * Contructors * * * * * * * * * * * * * * * * //

Foam::gradientSearch::gradientSearch(const Foam::fvMesh& mesh)
:
mesh_(mesh),
owner_(mesh.owner()),
neighbour_(mesh.neighbour()),
neighbourFaceCells_(0),
visited_(labelList(10,-1)),
cellTreePtrs_(List< autoPtr<indexedOctree<treeDataCell> > >(mesh_.cellZones().size()+1)),
zoneBoundPtrs_(List< autoPtr<boundBox> >(mesh_.cellZones().size()+1))
{
    if(Pstream::parRun())
    {
        neighbourFaceCells_.setSize(Pstream::nProcs(),labelList(0));
        const polyBoundaryMesh& bm = mesh_.boundaryMesh();
        forAll(bm, patchI)
        {
            if
            (
                isA<processorPolyPatch>(bm[patchI])
                &&
                !isA<processorCyclicPolyPatch>(bm[patchI])
            )
            {
                const processorPolyPatch& np =
                    refCast<const processorPolyPatch>(bm[patchI]);
                neighbourFaceCells_[np.neighbProcNo()] = np.faceCells();
            }
        }

        // scatter proc hit data
        PstreamBuffers pBufs(Pstream::nonBlocking);

        const label nProcs = Pstream::nProcs();
        const label myProcNo = Pstream::myProcNo();

        for (label procI = 0; procI < nProcs; procI++)
        {
            if (procI != myProcNo)
            {
                UOPstream toDomain(procI, pBufs);
                toDomain << neighbourFaceCells_[procI];
            }
        }

        pBufs.finishedSends();

        for (label procI = 0; procI < nProcs; procI++)
        {
            if (procI != myProcNo)
            {
                UIPstream str(procI, pBufs);
                str >> neighbourFaceCells_[procI];
            }
        }
    }
}

// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * //

Foam::gradientSearch::~gradientSearch() {}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::label Foam::gradientSearch::searchLocal
(
    Foam::searchItem& item,
    label fromProc = -1
) const
{
//     bool log = false;
//     autoPtr<std::ofstream> osPtr;
//     const vectorField& ccs = mesh_.cellCentres();
//
//     if(item.procID()==1 && item.cellLabel()==22414)
//     {
//         log = true;
//         Pout<<"Searching for item..."<<endl;
//         fileName path = mesh_.time().rootPath()+"/"+mesh_.time().globalCaseName(); // mesh_.time().systemPath();
//         osPtr.set(new std::ofstream((path+"/searchPath.obj").c_str(),std::ofstream::app | std::ofstream::out));
//         osPtr()
//             <<"#starting search for cell at "
//             <<item.position().x()<<" "<<item.position().y()<<" "
//             <<item.position().z()<<" "<<" with seed "<<item.seed()
//             <<" on proc "<<Pstream::myProcNo()
//             <<std::endl;
//     }

    const Foam::point& curPoint(item.position());
    // Get seed cell
    label& seedCell = item.seed();

    clearVisited();

    if(seedCell >= 0)
    {
        // Collect mesh data
        const label nIntFaces = mesh_.nInternalFaces();

        const cellList& cells(mesh_.cells());
        const vectorField& faceCentres_(mesh_.faceCentres());
        const vectorField& faceNormals_(mesh_.faceAreas());

        const faceList& faces = mesh_.faces();
        const vectorField& cc = mesh_.cellCentres();
        const vectorField& points = mesh_.points();

        // Set value for proc hit checking
        label neighbProcNo = -1;
        label procHit = -2;

        // We do not have the containing cell yet...
        bool found = false;

        // ...but there is a cell to look at.
        bool nextCell = true;

        while(!found && nextCell)
        {
            // Asume this is the right cell until opposit proven
            bool inside = true;
            bool keepGoing = false;

            // But we do not know any next cell to look at if we fail
            nextCell = false;

            const labelList& cellFaces = cells[seedCell];
            if(cellFaces.size() == 0)
            {
                inside = false;
            }

            const vector& curCentre = cc[seedCell];
            const vector CP = curPoint - curCentre;
            for
            (
                label faceI = 0;
                !nextCell && !found && (keepGoing || inside)
             && ( faceI < cellFaces.size() );
                faceI++
            )
            {
                const label curFace = cellFaces[faceI];

                // Asume current cell is owner of the face. This was isOwner
                // is true for boundary faces
                bool isOwner = true;
                if(curFace<nIntFaces)
                {
                    // Check ownership for internal faces
                    isOwner = (seedCell == owner_[curFace]);
                }

                bool isInPyramid = true;
                const face& f = faces[curFace];
                for(label pointI = 0; pointI < f.size() && isInPyramid; pointI++)
                {
                    const vector CS = points[f[pointI]] - curCentre;
                    const vector CE = points[f.nextLabel(pointI)] - curCentre;
                    if( ( ( ( CE^CS ) & CP ) > 0.0 ) == isOwner )
                    {
                        isInPyramid = false;
                    }
                }

                if(isInPyramid)
                {
                    // if projection is zero down to machine precision,
                    // the point is on owner side of the face
                    const scalar projection =
                        (curPoint - faceCentres_[curFace]) & faceNormals_[curFace];
                    if(isOwner)
                    {
                        if
                        (
                            projection > SMALL
                        )
                        {
                            // In owner cell and projection is greater than zero
                            // Point is on other side of the face, so the point is
                            // not inside the current cell
                            inside = false;
                            if(curFace < nIntFaces && !visited(neighbour_[curFace]))
                            {
    //                             if(log)
    //                             {
    //                                 const vector& old = ccs[seedCell];
    //                                 const vector& neW = ccs[neighbour_[curFace]];
    //                                 osPtr()
    //                                     <<"# from cell "<<seedCell<<" on proc "
    //                                     <<Pstream::myProcNo()<<":\n"
    //                                     <<"v "<<old.x()<<" "<<old.y()<<" "<<old.z()
    //                                     <<"\n"
    //                                     <<"# to cell "<<neighbour_[curFace]
    //                                     <<" on proc "<<Pstream::myProcNo()<<":\n"
    //                                     <<"v "<<neW.x()<<" "<<neW.y()<<" "<<neW.z()
    //                                     <<"\nl -2 -1"<<std::endl;
    //                             }
                                seedCell = neighbour_[curFace];
                                setVisited(seedCell);
                                nextCell = true;
                            }
                        }
                        else
                        {
                            found = true;
                        }
                    }
                    else
                    {
                        if
                        (
                            projection <= 0
                        )
                        {
                            // In neighbour cell and projection is not greater than
                            // zero. Point is on other side of the face
                            inside = false;
                            if(curFace < nIntFaces && !visited(owner_[curFace]))
                            {
    //                             if(log)
    //                             {
    //                                 const vector& old = ccs[seedCell];
    //                                 const vector& neW = ccs[owner_[curFace]];
    //                                 osPtr()
    //                                     <<"# from cell "<<seedCell<<" on proc "
    //                                     <<Pstream::myProcNo()<<":\n"
    //                                     <<"v "<<old.x()<<" "<<old.y()<<" "<<old.z()
    //                                     <<"\n"
    //                                     <<"# to cell "<<owner_[curFace]
    //                                     <<" on proc "<<Pstream::myProcNo()<<":\n"
    //                                     <<"v "<<neW.x()<<" "<<neW.y()<<" "<<neW.z()
    //                                     <<"\nl -2 -1"<<std::endl;
    //                             }
                                seedCell = owner_[curFace];
                                setVisited(seedCell);
                                nextCell = true;
                            }
                        }
                        else
                        {
                            found = true;
                        }
                    }

                    if(!inside && (curFace >= nIntFaces))
                    {
                        if(Pstream::parRun())
                        {
                            // Point is behind the current face which is a boundary
                            // face. If this is a processor face, we might want to
                            // continue on the other processor

                            const label patchID =
                                mesh_.boundaryMesh().whichPatch(curFace);
                            if
                            (
                                isA<processorPolyPatch>
                                    (mesh_.boundaryMesh()[patchID])
                                &&
                                !isA<processorCyclicPolyPatch>
                                    (mesh_.boundaryMesh()[patchID])
                            )
                            {
                                const processorPolyPatch& interface =
                                    dynamic_cast<const processorPolyPatch&>
                                    (
                                        mesh_.boundaryMesh()[patchID]
                                    );

                                // Don't go back to the proc you came from
                                if(interface.neighbProcNo() != fromProc)
                                {
    //                                 if(log)
    //                                 {
    //                                     const vector& old = ccs[seedCell];
    //                                     const vector& neW = interface.neighbFaceCellCentres()[curFace-interface.start()];
    //                                     osPtr()
    //                                         <<"# proc jump\n"
    //                                         <<"# from cell "<<seedCell<<" on proc "
    //                                         <<Pstream::myProcNo()<<":\n"
    //                                         <<"v "<<old.x()<<" "<<old.y()<<" "<<old.z()
    //                                         <<"\n"
    //                                         <<"# accross face "<<curFace-interface.start()
    //                                         <<" to proc "<<interface.neighbProcNo()<<":\n"
    //                                         <<"v "<<neW.x()<<" "<<neW.y()<<" "<<neW.z()
    //                                         <<"\nl -2 -1"<<std::endl;
    //                                 }
                                    neighbProcNo = interface.neighbProcNo();

                                    //- Cell label in neighbour domain behind the
                                    //  processor patch face
                                    procHit =
                                        neighbourFaceCells_[neighbProcNo]
                                            [curFace-interface.start()];
                                }
                            }
                        }

                        // Hit a boundary face, so we want to continue and hope to
                        // "creep" around the obstacle...
                        // even if it is a processor boundary. Go as far on this
                        // proc as possible.
                        keepGoing = true;
                    }
                }
                else if(faceI == cellFaces.size()-1)
                {
                    inside = false;
                }
            } // End of loop over all faces
            found = inside;

            // First cell completed, next cell will be a different one, so we
            // could jump back to the proc we came from
            fromProc = -1;
        }

        if(found)
        {
            // Maybe one wants to insert zoneChecking here, if dealing with special
            // cases like periodic overset grids or other stuff without clearly
            // separated cell zones
            return Pstream::myProcNo();
        }
        if( !found && (procHit > -1) )
        {
            seedCell = procHit;
            return neighbProcNo;
        }
    }

    return -1;
}

bool Foam::gradientSearch::visited(const Foam::label cellI) const
{
    const label* const __restrict__ vPtr = visited_.begin();
    for(label visitedI = visited_.size()-1; visitedI > 0; visitedI--)
    {
        if( vPtr[visitedI] == cellI) return true;
    }
    return false;
}

void Foam::gradientSearch::setVisited(const Foam::label cellI) const
{
    label* __restrict__ vPtr = visited_.begin();
    const label nVisited = visited_.size();
    for(label visitedI = nVisited - 1 ; visitedI > 0; visitedI--)
    {
        vPtr[visitedI] = vPtr[visitedI-1];
    }
    vPtr[0] = cellI;
}

void Foam::gradientSearch::clearVisited() const
{
    label* __restrict__ vPtr = visited_.begin();
    const label nVisited = visited_.size();
    for(label visitedI = nVisited - 1 ; visitedI > 0; visitedI--)
    {
        vPtr[visitedI] = -1;
    }
}

Foam::labelList Foam::gradientSearch::zoneSeeds(const label zone = -1) const
{
    labelList seedSeeds(3,-1);
    if(zone != -1)
    {
        const labelList& zoneCells =
            mesh_.cellZones()[zone];
        const label nZoneCells = zoneCells.size();
        if(nZoneCells > 2)
        {
            seedSeeds[2]=zoneCells[(nZoneCells/2)];
        }
        else
        {
            seedSeeds.setSize(2);
        }

        if(nZoneCells > 1)
        {
            seedSeeds[1]=zoneCells[nZoneCells-1];
        }
        else
        {
            seedSeeds.setSize(1);
        }

        if(nZoneCells > 0)
        {
            seedSeeds[0]=zoneCells[0];
        }
        else
        {
            seedSeeds.setSize(0);
        }
    }
    else
    {
        const label nCells = mesh_.nCells();
        if(nCells > 2)
        {
            seedSeeds[2] = nCells/2;
        }
        else
        {
            seedSeeds.setSize(2);
        }

        if(nCells > 1)
        {
            seedSeeds[1] = nCells-1;
        }
        else
        {
            seedSeeds.setSize(1);
        }

        if(nCells > 0)
        {
            seedSeeds[0] = 0;
        }
        else
        {
            seedSeeds.setSize(0);
        }

    }
    return seedSeeds;

}

// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

void Foam::gradientSearch::search
(
    Foam::autoPtr< Foam::List< Foam::searchItem > >& searchableItemsPtr
) const
{
    // Search for items
    // 1. Perform gradient search from seed cell
    // 2. If this fails, send items to processors with point in bounding box
    // 3. Perform octree search with these items
    // 4. Check which items were found
    // 5. Send still failed items back to their proc

    const label nLocalItems = searchableItemsPtr().size();
    label nGlobalItems = nLocalItems;

    Pstream::gather(nGlobalItems,sumOp<scalar>());
    Pstream::scatter(nGlobalItems);

    const label myProcNo = Pstream::myProcNo();
    const label nProcs = Pstream::nProcs();

    if(debug)
    {
        Pout<<"Starting gradient search for "<<searchableItemsPtr().size()
            <<" items."<<nl<<endl;
    }

    label nBruteForced = 0;

    if(failItemsPtr_.valid() && failItemsPtr_().size() > 0)
    {
        WarningIn
        (
            "void Foam::gradientSearch::search("
            "Foam::autoPtr< Foam::List< Foam::searchItem > >& "
            "searchableItemsPtr) const"
        )
            <<"There were "<<failItemsPtr_().size()<<" failed items in the "
            <<"previous search, which have not been handled and will be "
            <<"deleted now."<<endl;
        failItemsPtr_.clear();
    }

    failItemsPtr_.reset( new List<searchItem>(0) );

    List<searchItem>& failItems = failItemsPtr_();

    // Reference to search items in question
    List<searchItem>& searchableItems = searchableItemsPtr();

    // Global success of search algorithm
    labelListList finallyFoundPerProc(nProcs, labelList(0));
    finallyFoundPerProc[myProcNo].setSize(nLocalItems, 0);

    PstreamBuffers* pBufs = new PstreamBuffers(Pstream::nonBlocking);

    for (label procI = 0; procI < nProcs; procI++)
    {
        if (procI != myProcNo)
        {
            UOPstream toDomain(procI, *pBufs);
            toDomain << nLocalItems;
        }
    }

    pBufs->finishedSends();

    for (label procI = 0; procI < nProcs; procI++)
    {
        if (procI != myProcNo)
        {
            UIPstream str(procI, *pBufs);
            label nItems;
            str >> nItems;
            finallyFoundPerProc[procI].setSize(nItems,0);
        }
    }

    delete pBufs;
    pBufs = NULL;

    // Note group and running index
    forAll(searchableItems, itemI)
    {
        searchableItems[itemI].groupID() = myProcNo;
        searchableItems[itemI].groupIndex() = itemI;
    }

    // Store initial items to communicate failed items back to original proc
    const List<searchItem> initialItems(searchableItems);

    // Linked list of items found on this proc
    SLList< List<searchItem>* > localItemList;

    // Linked list of failed items (buffer)
    SLList< searchItem > failedItemsBuffer;


    bool anyProcHit = false;

    // Store the label, the item came from
    labelList hitFromProc(searchableItems.size(), -1);

    if(searchableItems.size() > 0)
    {
        // There are some points which hit us from this proc...
        anyProcHit = true;
    }

    // Perform gradient search with proc passing
    // Items which have no cell are pushed to failedItemsBuffer
    for
    (
        Pstream::gather<bool>(anyProcHit, orOp<bool>()),
            Pstream::scatter<bool>(anyProcHit);
        anyProcHit;
        Pstream::gather<bool>(anyProcHit, orOp<bool>()),
            Pstream::scatter<bool>(anyProcHit)
    )
    {
        label nSearchableItems = searchableItems.size();

        // Counting differnt result types
        label nLocalItems = 0;
        label nFailItems = 0;

        // Proc hits per other proc
        labelList hitsPerProc(nProcs, 0);

        // Store returned proc ID for each search
        labelList procNos(nSearchableItems);

        forAll(searchableItems,itemI)
        {
            // Perform local search for this point and store returned proc ID
            label& procNo = procNos[itemI];

            if(debug>1)
            {
                Pout<<"Searching for item "<<itemI<<" / "<<searchableItems.size()<<endl;
            }
            procNo = searchLocal(searchableItems[itemI], hitFromProc[itemI]);

            // Count the different results
            if(procNo == myProcNo)
            {
                nLocalItems++;
            }
            else if(procNo > -1)
            {
                hitsPerProc[procNo]++;
            }
            else
            {
                nFailItems++;
            }
        }

    // Allocate lists for different search results
    // For each proc to be able to scatter and gatter data

        // Append local results, will be filled from proc hits at interfaces to this
        // proc on other procs
        localItemList.append(new List<searchItem>(nLocalItems));
        List<searchItem>& localItems = *localItemList.last();

        // Allocate memory for hits per proc
        List< List< searchItem > > allProcHitItems(nProcs);
        forAll(hitsPerProc, hitsPerProcI)
        {
            label& hpp = hitsPerProc[hitsPerProcI];
            allProcHitItems[hitsPerProcI].setSize( hpp );
            hpp = 0;
        }


        nLocalItems = 0;
        nFailItems = 0;

        anyProcHit = false;

        // Copy items to appropriate list
        forAll(searchableItems, itemI)
        {
            const label procNo = procNos[itemI];
            if(procNo == myProcNo)
            {
                // Copy local item to local item list
                const searchItem& curItem = searchableItems[itemI];

                localItems[nLocalItems] = curItem;
                nLocalItems++;
                label& success =
                    finallyFoundPerProc[curItem.groupID()][curItem.groupIndex()];

                if(success == 0)
                {
                    success = 1;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void Foam::gradientSearch::search"
                        "("
                            "Foam::autoPtr< Foam::List< Foam::searchItem > >& "
                            "searchableItemsPtr"
                        ") const"
                    )
                    << "Multiple success for item "<<curItem
                    << exit(FatalError);
                }
            }
            else if(procNo != -1)
            {
                // Copy item with proc hit to corresponding proc list
                label& hpp = hitsPerProc[procNo];
                allProcHitItems[procNo][hpp]=searchableItems[itemI];
                hpp++;
                anyProcHit = true;
            }
            else
            {
                failedItemsBuffer.append(searchableItems[itemI]);
            }
        }

        Pstream::gather<bool>(anyProcHit,orOp<bool>());
        Pstream::scatter<bool>(anyProcHit);

        if(anyProcHit)
        {
            // scatter proc hit data
            PstreamBuffers pBufs(Pstream::nonBlocking);

            for (label procI = 0; procI < nProcs; procI++)
            {
                if (procI != myProcNo)
                {
                    UOPstream toDomain(procI, pBufs);
                    toDomain << allProcHitItems[procI];
                }
            }

            pBufs.finishedSends();

            // gather all the proc hits
            label procHitsItemCount = 0;

            for (label procI = 0; procI < nProcs; procI++)
            {
                if (procI != myProcNo)
                {
                    UIPstream str(procI, pBufs);
                    str >> allProcHitItems[procI];
                    procHitsItemCount += allProcHitItems[procI].size();
                }
            }

            // Combine the received proc hits
            searchableItems.setSize(procHitsItemCount);
            hitFromProc.setSize(procHitsItemCount);
            procHitsItemCount = 0;
            for (label procI = 0; procI < nProcs; procI++)
            {
                if (procI != myProcNo)
                {
                    const List<searchItem>& curProcHits = allProcHitItems[procI];
                    forAll(curProcHits, hitI)
                    {
                        searchableItems[procHitsItemCount] = curProcHits[hitI];

                        hitFromProc[procHitsItemCount] = procI;
                        procHitsItemCount++;
                    }
                }
            }
        }

    }

    // Mark successful items
    collectGlobalSuccess(finallyFoundPerProc);

    // Gradient search performed
    // Handle failed items
    failItems = failedItemsBuffer;
    failedItemsBuffer.clear();
    label nFailItems = failItems.size();
    bool anyFailItems = nFailItems > 0;

    Pstream::gather<bool>(anyFailItems,orOp<bool>());
    Pstream::scatter<bool>(anyFailItems);

    if(anyFailItems)
    {
        if(debug)
        {
            Pout<<"There were "<<nFailItems<<" failed items"<<endl;
        }

        // Force tet decomposition on all procs simultaniously
        mesh_.tetBasePtIs();

        // Create bouding boxes
        if(bBoxesPtr_.empty())
        {
            bBoxesPtr_.set(new List<boundBox>(nProcs));
            List<boundBox>& bBoxes = bBoxesPtr_();
            const boundBox& myBounds = mesh_.bounds();
            bBoxes[myProcNo] = myBounds;

    // communicate boundBoxes

            pBufs = new PstreamBuffers(Pstream::nonBlocking);

            for (label procI = 0; procI < nProcs; procI++)
            {
                if (procI != myProcNo)
                {
                    UOPstream toDomain(procI, *pBufs);
                    toDomain << myBounds;
                }
            }

            pBufs->finishedSends();

            for (label procI = 0; procI < nProcs; procI++)
            {
                if (procI != myProcNo)
                {
                    UIPstream str(procI, *pBufs);
                    str >> bBoxes[procI];
                }
            }

            delete pBufs;
            pBufs = NULL;
        }

        List<boundBox>& bBoxes = bBoxesPtr_();

        // Pointer to PstreamBuffers
        pBufs = new PstreamBuffers(Pstream::nonBlocking);

        // Collect all failed items in mesh bounding box
        List< List< searchItem > > allFailItems(nProcs);
        {
            // Local failed items found in the boundBox of each domain
            labelListList itemsInDomain(nProcs);

            for (label procI = 0; procI < nProcs; procI++)
            {
                const boundBox& curBox = bBoxes[procI];
                label nItems = 0;
                boolList inBox(failItems.size(), false);

                forAll(failItems,itemI)
                {
                    if( curBox.contains(failItems[itemI].position()) )
                    {
                        inBox[itemI] = true;
                        nItems++;
                    }
                }

                labelList& curList = itemsInDomain[procI];
                curList.setSize(nItems);
                nItems = 0;
                forAll(failItems,itemI)
                {
                    if(inBox[itemI])
                    {
                        curList[nItems] = itemI;
                        nItems++;
                    }
                }

                if(procI == myProcNo)
                {
                    List<searchItem>& curFailItems = allFailItems[procI];
                    curFailItems.setSize(nItems);
                    forAll(curFailItems,itemI)
                    {
                        curFailItems[itemI] = failItems[curList[itemI]];
                    }
                }
            }

            for(label procI = 0; procI < nProcs; procI++)
            {
                if (procI != myProcNo)
                {
                    UOPstream toDomain(procI, *pBufs);
                    const labelList& curList = itemsInDomain[procI];
                    toDomain << curList.size();
                    forAll(curList,itemI)
                    {
                        toDomain << failItems[curList[itemI]];
                    }
                }
            }

            pBufs->finishedSends();

    // Collect global failed items

            // Receive the failed items from foreign procs
            for (label procI = 0; procI < nProcs; procI++)
            {
                if (procI != myProcNo)
                {
                    UIPstream str(procI, *pBufs);
                    label nItems = -1;
                    str >> nItems;
                    List<searchItem>& curFailItems = allFailItems[procI];
                    curFailItems.setSize(nItems);
                    for(label itemI = 0; itemI < nItems; itemI++)
                    {
                        str >> curFailItems[itemI];
                    }
                }
            }
        }

        delete pBufs;
        pBufs = NULL;

        // Clear failItems as they are in allFailItems of their new procs
        failItems.setSize(0);

        // Check all the failed items if there are cells in the zones
        for(label procI = 0; procI < nProcs; procI++)
        {
            // Failed items from procI
            List<searchItem>& curFailItems = allFailItems[procI];

            label nSuccessfulItems = 0;

            boolList successfulItems(curFailItems.size(),false);

            if(debug)
            {
                Pout<<"Starting octree search for "<<curFailItems.size()
                    <<" items."<<endl;
            }

            forAll(curFailItems, itemI)
            {
                searchItem& curItem = curFailItems[itemI];
                label& curFound =
                    finallyFoundPerProc[curItem.groupID()][curItem.groupIndex()];

                const label zoneID = curItem.zoneID();

                if
                (
                    curFound == 0
                    &&
                    (
                        (
                            zoneID == -1
                            &&
                            mesh_.bounds().contains(curItem.position())
                        )
                        ||
                        zoneBounds(zoneID).contains(curItem.position())
                    )
                )
                {
                    const indexedOctree<treeDataCell>& tree =
                        cellTree(curItem.zoneID());

                    const label res = tree.findInside(curItem.position());
                    if(res >= 0)
                    {
                        curItem.seed() = zoneID < 0
                                        ?
                                        res
                                        :
                                        mesh_.cellZones()[zoneID][res];
                        curFound = 1;
                        successfulItems[itemI]=true;
                        nSuccessfulItems++;
                    }
                }
                else if (curFound > 0)
                {
                    FatalErrorIn
                    (
                        "void Foam::gradientSearch::search"
                        "("
                            "Foam::autoPtr< Foam::List< Foam::searchItem > >& "
                            "searchableItemsPtr"
                        ") const"
                    )
                    << "Item already found and still a failed item "<<curItem
                    << exit(FatalError);
                }
            }

            if(debug)
            {
                Pout<<"Found "<<nSuccessfulItems<<" items using octree "
                    <<"search from proc "<<procI<<"."<<endl;
            }

            nBruteForced += nSuccessfulItems;

            localItemList.append(new List<searchItem>(nSuccessfulItems));
            List<searchItem>& localItems = *localItemList.last();

            nSuccessfulItems = 0;
            forAll(curFailItems, itemI)
            {
                if(successfulItems[itemI])
                {
                    localItems[nSuccessfulItems++] = curFailItems[itemI];
                }
            }
        }

        // Mark items found during octree search
        collectGlobalSuccess(finallyFoundPerProc);

        // Count items without success per original processor
        labelList nFailedItemsPerProc(nProcs, 0);
        const labelList& finallyFound = finallyFoundPerProc[myProcNo];
        forAll(initialItems, itemI)
        {
            if(finallyFound[itemI] == 0)
            {
                nFailedItemsPerProc[initialItems[itemI].procID()]++;
            }
        }

        labelListList failedItemsPerProc(nProcs);
        for(label procI = 0; procI < nProcs; procI++)
        {
            failedItemsPerProc[procI].setSize(nFailedItemsPerProc[procI]);
            nFailedItemsPerProc[procI] = 0;
        }

        forAll(initialItems, itemI)
        {
            if(finallyFound[itemI] == 0)
            {
                const label procID = initialItems[itemI].procID();
                failedItemsPerProc[procID][nFailedItemsPerProc[procID]++] =
                    itemI;
            }
        }

        nFailItems = nFailedItemsPerProc[myProcNo];

        pBufs = new PstreamBuffers(Pstream::nonBlocking);

        for (label procI = 0; procI < nProcs; procI++)
        {
            if (procI != myProcNo)
            {
                UOPstream toDomain(procI, *pBufs);
                toDomain << nFailedItemsPerProc[procI];
            }
        }

        pBufs->finishedSends();

        for (label procI = 0; procI < nProcs; procI++)
        {
            if (procI != myProcNo)
            {
                UIPstream str(procI, *pBufs);
                label nItems;
                str >> nItems;
                nFailItems += nItems;
            }
        }

        delete pBufs;
        pBufs = NULL;

        failItems.setSize(nFailItems);
        nFailItems = 0;

        pBufs = new PstreamBuffers(Pstream::nonBlocking);

        for (label procI = 0; procI < nProcs; procI++)
        {
            const labelList& curFailedItems = failedItemsPerProc[procI];
            if (procI != myProcNo)
            {
                UOPstream toDomain(procI, *pBufs);
                toDomain << curFailedItems.size();
                forAll(curFailedItems, itemI)
                {
                    toDomain << initialItems[curFailedItems[itemI]];
                }
            }
            else
            {
                forAll(curFailedItems, itemI)
                {
                    failItems[nFailItems++] = initialItems[curFailedItems[itemI]];
                }
            }
        }

        pBufs->finishedSends();

        for (label procI = 0; procI < nProcs; procI++)
        {
            if (procI != myProcNo)
            {
                UIPstream str(procI, *pBufs);
                label nItems;
                str >> nItems;
                for(label itemI = 0; itemI < nItems; itemI++)
                {
                    str >> failItems[nFailItems++];
                }
            }
        }

        delete pBufs;
        pBufs = NULL;

        if(nFailItems > 0 && debug)
        {
            WarningIn("void Foam::gradientSearch::search("
                "Foam::autoPtr< Foam::List< Foam::searchItem > >& "
                "searchableItemsPtr) const")
            <<"Got unhandled failed items."<<endl;
            Pout<<"Number of failed items: "<<nFailItems<<endl;
        }
    }

    label nFoundItems = 0;
    forAllConstIter(SLList< List<searchItem>* >, localItemList, iter)
    {
        nFoundItems += (**iter).size();
    }

    searchableItems.setSize(nFoundItems);

    nFoundItems = 0;
    forAllIter( SLList< List<searchItem>* >, localItemList, iter)
    {
        const List<searchItem>& curList = **iter;
        forAll(curList,itemI)
        {
            searchableItems[nFoundItems++] = curList[itemI];
        }
        delete *iter;
        *iter = NULL;
    }

    if(debug)
    {
        label nFound = nFoundItems;

        Pstream::gather(nFound, sumOp<scalar>());
        Pstream::scatter(nFound);
        Pstream::gather(nBruteForced, sumOp<scalar>());
        Pstream::scatter(nBruteForced);

        if(nGlobalItems > 0)
        {
            Info<<"gradientSearch finished, found "<<nFound<<" out of "<<nGlobalItems
                <<". "<<(100.0*nBruteForced/nGlobalItems)
                <<"% of them using brute force"<<endl;
        }
        if(debug>1)
        {
            Pout<<"Failed items: "<<failItems<<endl;
        }
    }
}

Foam::autoPtr< Foam::List< Foam::searchItem > >
Foam::gradientSearch::search
(
    const Foam::List<Foam::point>& points,
    const Foam::List<Foam::label>& seeds,
    const Foam::List<Foam::label>& zones = Foam::List<Foam::label>(0)
) const
{
    bool noZones = false;

    if(zones.size() == 0)
    {
       noZones = true;
    }
    if
    (
        ( points.size() != seeds.size() )
        ||
        ( ( points.size() != zones.size() ) && ( !noZones ) )
    )
    {
        FatalErrorIn
        (
            "const Foam::autoPtr<Foam::List<Foam::searchItem>> "
            "Foam::gradientSearch::search "
            "("
                "const Foam::List<Foam::point>& points, "
                "const Foam::List<Foam::label>& seeds"
            ") const"
        )
        << "Size of points, seed cell and zone labels do not match."
        << exit(FatalError);
    }

    // Allocate memory for result
    autoPtr< List< searchItem > > itemsPtr
    (
        new List<searchItem>(points.size())
    );
    List<searchItem>& items(itemsPtr());

    const label myProcNo = Pstream::myProcNo();

    forAll(items, itemI)
    {
        searchItem& curItem = items[itemI];
        curItem.procID() = myProcNo;
        curItem.cellLabel() = itemI;
        curItem.position() = points[itemI];
        curItem.seed() = seeds[itemI];
        if(noZones)
        {
            curItem.zoneID() = -1;
        }
        else
        {
            curItem.zoneID() = zones[itemI];
        }
    }

    this->search(itemsPtr);

    return itemsPtr;
}

Foam::autoPtr< Foam::List< Foam::searchItem > >
Foam::gradientSearch::search
(
    const Foam::List<Foam::point>& points
) const
{
    const labelList seeds(points.size(), 0);
    const labelList zones(points.size(), -1);

    autoPtr< List< searchItem > > tResult = search(points, seeds, zones);

//     List<searchItem>& result = tResult();

    return tResult;
}

Foam::autoPtr< Foam::searchItem >
Foam::gradientSearch::search
(
    const Foam::point& p,
    const Foam::label s = 0,
    const Foam::label z = -1
) const
{
    const List<point> points(1, p);
    const List<label> seeds(1, s);
    const List<label> zones(1, z);

    autoPtr< List< searchItem > > tResultList = search(points, seeds, zones);

    autoPtr<searchItem> tResult;
    searchItem& result = tResult();

    result = tResultList()[0];

    return tResult;
}

void Foam::gradientSearch::collectGlobalSuccess
(
    Foam::labelListList& finallyFoundPerProc
) const
{
    if(!Pstream::parRun())
    {
        return;
    }

    // Collect all the found items and check if any item handed to this proc
    // was not found on any proc
    PstreamBuffers *pBufs = new PstreamBuffers(Pstream::nonBlocking);

    const label myProcNo = Pstream::myProcNo();
    const label nProcs = Pstream::nProcs();

    for (label procI = 0; procI < nProcs; procI++)
    {
        if (procI != myProcNo)
        {
            UOPstream toDomain(procI, *pBufs);
            toDomain << finallyFoundPerProc[procI];
        }
    }

    pBufs->finishedSends();

    labelList& finallyFound = finallyFoundPerProc[myProcNo];

    for (label procI = 0; procI < nProcs; procI++)
    {
        if (procI != myProcNo)
        {
            labelList received;
            UIPstream str(procI, *pBufs);
            str >> received;
            if(received.size() != finallyFound.size())
            {
                FatalErrorIn
                (
                    "void Foam::gradientSearch::collectGlobalSuccess "
                    "( "
                        "const Foam::boolListList& finallyFoundPerProc "
                    ") const"
                )
                << "Size of received list (" << received.size()<< ") from proc "
                << procI << " doesn't match local size (" <<finallyFound.size()
                << ")."
                << exit(FatalError);
            }
            forAll(received, itemI)
            {
                if
                (
                    received[itemI] == 1
                )
                {
                    if(finallyFound[itemI] > 0 && !ignoreMultiple_)
                    {
                        FatalErrorIn
                        (
                            "void Foam::gradientSearch::collectGlobalSuccess "
                            "( "
                                "const Foam::boolListList& finallyFoundPerProc "
                            ") const"
                        )
                        << "Multiple results for item " << itemI << " proc "
                        << procI << "." << nl
                        << "Received: "<< received[itemI] << ", stored "
                        << finallyFound[itemI]
                        << exit(FatalError);
                    }
                    else
                    {
                        finallyFound[itemI] = 2;
                    }
                }
            }
        }
    }

    delete pBufs;
    pBufs = NULL;

    forAll(finallyFound,i)
    {
        if(finallyFound[i]>0)
        {
            finallyFound[i]=2;
        }
    }

    pBufs = new PstreamBuffers(Pstream::nonBlocking);

    for (label procI = 0; procI < nProcs; procI++)
    {
        if (procI != myProcNo)
        {
            UOPstream toDomain(procI, *pBufs);
            toDomain << finallyFound;
        }
    }

    pBufs->finishedSends();

    for (label procI = 0; procI < nProcs; procI++)
    {
        if (procI != myProcNo)
        {
            UIPstream str(procI, *pBufs);
            str >> finallyFoundPerProc[procI];
        }
    }

    delete pBufs;
    pBufs = NULL;
}

const Foam::indexedOctree< Foam::treeDataCell >&
Foam::gradientSearch::cellTree(const label zoneID) const
{
    const label idx = (zoneID == -1) ? mesh_.cellZones().size() : zoneID;

    if(cellTreePtrs_.size() < idx + 1)
    {
        FatalErrorIn
        (
            "const Foam::indexedOctree< Foam::treeDataCell >& "
            "Foam::gradientSearch::cellTree(const label zoneID)"
        )
        << "Requested cellTree with " << idx << ", but only "
        << cellTreePtrs_.size() << " trees defined."
        << exit(FatalError);
    }

    autoPtr< indexedOctree<Foam::treeDataCell> >& treePtr = cellTreePtrs_[idx];

    if(treePtr.empty())
    {
        if(zoneID<0)
        {
            treePtr.reset(
                new indexedOctree<treeDataCell>
                    (
                        treeDataCell
                        (
                            false,
                            mesh_,
                            //polyMesh::FACE_CENTRE_TRIS
                            polyMesh::CELL_TETS
                        ),
                        treeBoundBox(zoneBounds(zoneID)),
                        8,              // maxLevel
                        10,             // leafsize
                        6.0             // duplicity
                    )
            );
        }
        else
        {
            treePtr.reset(
                new indexedOctree<treeDataCell>
                    (
                        treeDataCell
                        (
                            false,
                            mesh_,
                            mesh_.cellZones()[zoneID],
                            //polyMesh::FACE_CENTRE_TRIS
                            polyMesh::CELL_TETS
                        ),
                        treeBoundBox(zoneBounds(zoneID)),
                        8,              // maxLevel
                        10,             // leafsize
                        6.0             // duplicity
                    )
            );
        }
    }

    return treePtr();
}

const Foam::boundBox& Foam::gradientSearch::zoneBounds(const label zoneID) const
{
    const label idx = (zoneID == -1) ? mesh_.cellZones().size() : zoneID;

    if(zoneBoundPtrs_.size() < idx + 1)
    {
        FatalErrorIn
        (
            "const Foam::boundBox& "
            "Foam::gradientSearch::zoneBounds(const label zoneID) const"
        )
        << "Requested bounding box for zone with index with " << idx
        << ", but only " << zoneBoundPtrs_.size() << " boxes defined."
        << exit(FatalError);
    }

    autoPtr<boundBox>& zoneBoundPtr = zoneBoundPtrs_[idx];

    if(zoneBoundPtr.empty())
    {
        if(zoneID < 0)
        {
            zoneBoundPtr.reset(new boundBox(mesh_.bounds()));
        }
        else
        {
            if(mesh_.cellZones()[zoneID].size() > 0)
            {
                const labelList& zoneCells =
                    mesh_.cellZones()[zoneID];
                const labelListList& cellPnts = mesh_.cellPoints();
                const pointField& pts = mesh_.points();

                // Initial bounds of the cell zone
                vector min_(pts[cellPnts[zoneCells[0]][0]]);
                vector max_(min_);
                forAll(zoneCells,cellI)
                {
                    const labelList& curPnts
                        = cellPnts[zoneCells[cellI]];
                    forAll(curPnts,pntI)
                    {
                        const point& curPnt = pts[curPnts[pntI]];
                        min_=min(min_,curPnt);
                        max_=max(max_,curPnt);
                    }
                }
                zoneBoundPtr.reset(new boundBox(min_,max_));
            }
            else
            {
                // Empty box
                zoneBoundPtr.reset(new boundBox(vector::one, -vector::one));
            }
        }
    }

    return zoneBoundPtr();
}
