/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bellerophonInterpolation.H"
#include "bellerophonInterface.H"
#include "bellerophonInterfaceField.H"
#include "bellerophonLduMatrix.H"
#include "gradientSearch.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "faceSet.H"
#include "cellSet.H"
#include "OFstream.H"

#include "Time.H"
#include "interpolationItem.H"
#include "processorFvPatch.H"

void Foam::bellerophon::updateTopology() const
{
    // If this is used for meshes with topological changes, maybe this can be
    // improved by only updating/resetting patches with change.

    if(topologyUpToDate_)
    {
        return;
    }

    if(debug)
    {
        Info<<"Updating topology of Bellerophon interpolation."<<nl<<endl;
    }

    // Clear interpolation data
    clearInterpolation();

// Get the distinct list of Acceptor cells

    const fvMesh& msh = mesh();

    // Number of internal faces
    const label nIntFaces = msh.nInternalFaces();

    const labelList& own = msh.owner();
    const labelList& nei = msh.neighbour();

    // Is the cell an acceptor cell?
    labelList cellMarked(msh.nCells(),-1);

    // Maps internal face to patch if is interface
    labelList faceToInterface(nIntFaces,-1);

    // Flip map for interface faces
    boolList interfaceFlipMap(nIntFaces);

    // Number of interface faces per patch
    labelList nInterfaceFacesPerPatch(patchPtrs_.size(),0);

    forAll(patchPtrs_,ptrI)
    {
        // Reference to the patch
        const bellerophonInterface& p = *(patchPtrs_[ptrI]);

        // Holds the index of the interpolation cell in the list of
        // interpolation cell labels
        labelList patchInterpolationCells(msh.nCells(),0);

        // Get the cell labels of acceptor cells
        const labelList& acceptorFaceCells = p.faceCells();

        forAll(acceptorFaceCells,faceI)
        {
            label curAcceptor=acceptorFaceCells[faceI];
            if(cellMarked[curAcceptor] == -1)
            {
                cellMarked[curAcceptor] = ptrI;
                nAcceptors_++;
            }
        }

        for(label faceI = 0; faceI < nIntFaces; faceI++)
        {
            const label o = own[faceI];
            const label n = nei[faceI];
            if
            (
                cellMarked[o] == ptrI
                &&
                cellMarked[n] == -1
            )
            {
                faceToInterface[faceI]=ptrI;
                interfaceFlipMap[faceI]=true;
                nInterfaceFacesPerPatch[ptrI]++;
            }
            else if
            (
                cellMarked[o] == -1
                &&
                cellMarked[n] == ptrI
            )
            {
                faceToInterface[faceI]=ptrI;
                interfaceFlipMap[faceI]=false;
                nInterfaceFacesPerPatch[ptrI]++;
            }
        }
    }

    interfaceFaces_.setSize(patchPtrs_.size());
    interfaceFlipMap_.setSize(patchPtrs_.size());

    forAll(patchPtrs_,ptrI)
    {
        interfaceFaces_[ptrI].setSize(nInterfaceFacesPerPatch[ptrI]);
        interfaceFlipMap_[ptrI].setSize(nInterfaceFacesPerPatch[ptrI]);
        nInterfaceFacesPerPatch[ptrI] = 0;
    }

    forAll(faceToInterface, faceI)
    {
        const label ptrI = faceToInterface[faceI];
        if(ptrI != -1)
        {
            const label nInterfaceFace = nInterfaceFacesPerPatch[ptrI]++;
            interfaceFaces_[ptrI][nInterfaceFace] = faceI;
            interfaceFlipMap_[ptrI][nInterfaceFace] = interfaceFlipMap[faceI];
        }
    }

    nInterfaceFacesPerPatch.setSize(0);
    interfaceFlipMap.setSize(0);

    // Changing meaning of cell cellMarked:
    // After this, cellMarked holds the index of the acceptor/phi
    // interpolation cell in the corresponding list
    acceptorCells_.setSize(nAcceptors_);
    zoneIDs_.setSize(nAcceptors_);
    nAcceptors_ = 0;

    forAll(cellMarked,cellI)
    {
        if(cellMarked[cellI] != -1)
        {
            // Cell is acceptor cell
            acceptorCells_[nAcceptors_] = cellI;
            zoneIDs_[nAcceptors_] = patchPtrs_[cellMarked[cellI]]->donorZoneID();

            // Remember donor number for face to acceptor mapping
            cellMarked[cellI]=nAcceptors_;
            nAcceptors_++;
        }
    }

    const cellZoneMesh& cellZones = msh.cellZones();
    zoneSeeds_.setSize(cellZones.size());
    forAll(zoneSeeds_, zoneI)
    {
        if(cellZones[zoneI].size() > 0)
        {
            zoneSeeds_[zoneI]=cellZones[zoneI][0];
        }
        else
        {
            zoneSeeds_[zoneI]=-1;
        }
    }

    interpolationMethod_.clear();
    interpolationMethod_ = bellerophonInterpolationMethod::New
    (
        msh,
        primaryDonorCells_,
        acceptorCells_,
        deltas_,
        donorItemsPtr_,
        dict()
    );

    // Set topology up to date
    topologyUpToDate_ = true;
}

Foam::label Foam::bellerophon::addHole
(
    const Foam::labelList& map,
    const label donorZoneID,
    const label oversetPatchID
) const
{
    if(debug)
    {
        Info<<"Adding implicit hole."<<endl;
    }

    updateTopology();

    const fvMesh& msh = mesh();

    const labelList& own = msh.owner();
    const labelList& nei = msh.neighbour();

    // Count interface faces
    label interfaceFaceI = 0;

    forAll(own, faceI)
    {
        register const label o = own[faceI];
        register const label n = nei[faceI];
        if
        (
            (
                map[o] == REGULAR_CELL
                &&
                map[n] == ACCEPTOR_CELL
            )
            ||
            (
                map[o] == ACCEPTOR_CELL
                &&
                map[n] == REGULAR_CELL
            )
        )
        {
            interfaceFaceI++;
        }
    }

    interfaceFaces_.append(labelList(interfaceFaceI));
    interfaceFlipMap_.append(boolList(interfaceFaceI));

    labelList& newInterfaceFaces =
        interfaceFaces_[interfaceFaces_.size()-1];

    boolList& newInterfaceFlipMap =
        interfaceFlipMap_[interfaceFlipMap_.size()-1];

    interfaceFaceI = 0;

    forAll(own, faceI)
    {
        register const label o = own[faceI];
        register const label n = nei[faceI];
        if
        (
            map[o] == REGULAR_CELL
            &&
            map[n] == ACCEPTOR_CELL
        )
        {
            newInterfaceFaces[interfaceFaceI] = faceI;
            newInterfaceFlipMap[interfaceFaceI] = false;
            interfaceFaceI++;
        }
        else if
        (
            map[o] == ACCEPTOR_CELL
            &&
            map[n] == REGULAR_CELL
        )
        {
            newInterfaceFaces[interfaceFaceI] = faceI;
            newInterfaceFlipMap[interfaceFaceI] = true;
            interfaceFaceI++;
        }
    }

    label acceptorI = 0;
    label holeCellI = 0;

    forAll(map, cellI)
    {
        if(map[cellI] == ACCEPTOR_CELL)
        {
            acceptorI++;
        }
        else if (map[cellI] == HOLE_CELL)
        {
            holeCellI++;
        }
    }

    if(debug && mesh().time().outputTime())
    {
        cellSet acceptorCellSet
        (
            msh,
            "acceptorCells"+name(interfaceFaces_.size()-cifs_.size()),
            acceptorI
        );

        cellSet holeCellSet
        (
            msh,
            "holeCells"+name(interfaceFaces_.size()-cifs_.size()),
            holeCellI
        );

        forAll(map,cellI)
        {
            if(map[cellI] == ACCEPTOR_CELL) acceptorCellSet.insert(cellI);
            else if(map[cellI] == HOLE_CELL) holeCellSet.insert(cellI);
        }

        acceptorCellSet.instance() = msh.time().timeName();
        holeCellSet.instance() = msh.time().timeName();

        acceptorCellSet.write();
        holeCellSet.write();
    }

    acceptorCells_.setSize(nAcceptors_ + acceptorI);
    acceptorI = nAcceptors_;
    nAcceptors_ = acceptorCells_.size();

    labelList holeCells(holeCellI);
    holeCellI = 0;

    forAll(map, cellI)
    {
        if(map[cellI] == ACCEPTOR_CELL)
        {
            acceptorCells_[acceptorI++] = cellI;
        }
        else if (map[cellI] == HOLE_CELL)
        {
            holeCells[holeCellI++] = cellI;
        }
    }

    holeCells_.append(holeCells);
    zoneIDs_.setSize(nAcceptors_, donorZoneID);

    if(newLiveCellsMarkup_.size() == msh.nCells())
    {
        //- Clear marks on cells that are (still) hole or acceptor cells
        forAll(acceptorCells_, i) newLiveCellsMarkup_[acceptorCells_[i]] =false;
        forAll(holeCells_, i) newLiveCellsMarkup_[holeCells_[i]] = false;

        // Count number of new live cells
        label nNewLiveCells = 0;
        forAll(newLiveCellsMarkup_, i) if(newLiveCellsMarkup_[i]) nNewLiveCells++;

        // Resize list for found cells
        newLiveCells_.setSize(nNewLiveCells);
        nNewLiveCells = 0;

        // Populate list
        forAll(newLiveCellsMarkup_, i)
        {
            if(newLiveCellsMarkup_[i])
            {
                newLiveCells_[nNewLiveCells++] = i;
            }
        }
    }
    else
    {
        WarningIn("Foam::bellerophon::addHole(...)")
        <<"markup for new live cells has different size then the mesh. Was "
        <<"Foam::bellerophon::markHole() called?"
        <<endl;
    }

    // Generate donor to acceptor interpolation
    interpolationMethod_.clear();
    interpolationMethod_ = bellerophonInterpolationMethod::New
    (
        msh,
        primaryDonorCells_,
        acceptorCells_,
        deltas_,
        donorItemsPtr_,
        dict()
    );

    // Clear interpolation information
    clearInterpolation();

    const label holeInterface = interfaceFaces_.size()-1;

    // Set index of the hole interface
    refCast<const bellerophonInterface>
        (
            msh.boundary()[oversetPatchID]
        ).setHoleInterface(holeInterface);

    return holeInterface;

}

void Foam::bellerophon::updateInterpolation() const
{
    updateTopology();

    if(interpolationUpToDate_)
    {
        return;
    }

    if(debug)
    {
        Info<<"Updating interpolation of Bellerophon interpolation."<<nl<<endl;
    }

    // Reference to the mesh
    const fvMesh& mesh_ = mesh();

    // Create gradient search object
    gradientSearch gs(mesh_);

    // If we have no valid search items, generate them
    bool donorItemsInvalid = !donorItemsPtr_.valid();

    Pstream::gather(donorItemsInvalid,orOp<bool>());
    Pstream::scatter(donorItemsInvalid);

    // Collect cell centres for acceptor cells
    const vectorField& ccs = mesh_.cellCentres();

    if(donorItemsInvalid)
    {
        donorItemsPtr_.clear();


        // ID of this proc
        const label myProcNo = Pstream::myProcNo();

        if(!Pstream::parRun())
        {
            acceptorMapPtr_.reset(new labelList(mesh_.nCells(),-1));
            labelList& acceptorMap = acceptorMapPtr_();
            forAll(acceptorCells_,cellI)
            {
                acceptorMap[acceptorCells_[cellI]] = cellI;
            }
        }

        donorItemsPtr_.set(new List<searchItem>(nAcceptors_));
        List<searchItem>& donorItems = donorItemsPtr_();
        forAll(donorItems,itemI)
        {
            searchItem& curItem = donorItems[itemI];
            curItem.procID()    = myProcNo;
            curItem.cellLabel() = acceptorCells_[itemI];
            curItem.position()  = ccs[acceptorCells_[itemI]];
            curItem.seed()      = zoneSeeds_[zoneIDs_[itemI]];
            curItem.zoneID()    = zoneIDs_[itemI];
        }
    }

    // Find donor cells
    gs.search(donorItemsPtr_);

    // Local donor items
    const List<searchItem>& donorItems = donorItemsPtr_();

    if(Pstream::parRun())
    {
        // Sorting does not need to coincide, because the acceptors and
        // donors might be split across procs. Values will be sorted when
        // used...

        // Set primary donor cells
        nDonors_ = donorItems.size();
        primaryDonorCells_.setSize(nDonors_);
        forAll(donorItems,donorI)
        {
            const searchItem& curItem = donorItems[donorI];
            primaryDonorCells_[donorI] =
                curItem.seed();
        }
    }
    else
    {
        // Sorting must coincide to allow fast assignment

        // Cell to acceptor label mapping
        const labelList& acceptorMap = acceptorMapPtr_();

        nDonors_ = donorItems.size();
        if(nDonors_ != nAcceptors_)
        {
            boolList cellMask(mesh_.nCells(), false);
            forAll(acceptorCells_, cellI)
            {
                cellMask[acceptorCells_[cellI]] = true;
            }

            forAll(donorItems,donorI)
            {
                const label cellI = donorItems[donorI].cellLabel();
                cellMask[cellI] = false;
            }

            label nOrphanCells = 0;
            forAll(cellMask, cellI)
            {
                if(cellMask[cellI] == true)
                {
                    nOrphanCells++;
                }
            }

            Info<<"Found "<<nOrphanCells<< " orphan cells."<<endl;

            cellSet orphanCellsSet
            (
                mesh_,
                "orphanCells",
                nOrphanCells
            );

            // Acceptor cells
            forAll(cellMask,cellI)
            {
                if(cellMask[cellI] == true)
                {
                    orphanCellsSet.insert(cellI);
                }
            }

            orphanCellsSet.instance() = mesh_.time().timeName();
            orphanCellsSet.write();

            FatalErrorIn("void Foam::bellerophon::updateInterpolation()")
            <<"Number of found donors ("<<nDonors_
            <<")does not match number of acceptors("<<nAcceptors_<<")."<<nl
            <<"Wrote cell set orphanCells."
            <<exit(FatalError);
        }

        primaryDonorCells_.setSize(nDonors_,-1);
        forAll(donorItems,donorI)
        {
            const searchItem& curItem = donorItems[donorI];
            primaryDonorCells_[acceptorMap[curItem.cellLabel()]]
                = curItem.seed();
        }
    }


    deltas_.setSize(nDonors_);

    if(Pstream::parRun())
    {
        forAll(primaryDonorCells_,donorI)
        {
            deltas_[donorI] =
                donorItems[donorI].position()
                -
                ccs[primaryDonorCells_[donorI]];
        }
    }
    else
    {
        forAll(primaryDonorCells_,donorI)
        {
            deltas_[donorI] =
                ccs[acceptorCells_[donorI]]
                -
                ccs[primaryDonorCells_[donorI]];
        }
    }


    interpolationMethod_().update();

    interpolationUpToDate_ = true;

    if(debug>1)
    {
        volScalarField cellState
        (
            IOobject
            (
                "cellState",
                mesh_.thisDb().time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("null", dimensionSet(0,0,0,0,0,0,0), 0.0)
        );
        scalarField& cells = cellState.internalField();

        volVectorField cellCs
        (
            IOobject
            (
                "cellCentres",
                mesh_.thisDb().time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("null", dimensionSet(0,0,0,0,0,0,0), vector::zero)
        );
        const vectorField& ccs = mesh_.cellCentres();
        cellCs.internalField() = ccs;
        cellCs.write();

        cellSet acceptorCellSet
        (
            mesh_,
            "proc"+name(Pstream::myProcNo())+"acceptorCells",
            acceptorCells_.size()
        );

        cellSet donorCellSet
        (
            mesh_,
            "proc"+name(Pstream::myProcNo())+"donorCells",
            primaryDonorCells_.size()
        );

        // Acceptor cells
        forAll(acceptorCells_,acceptorI)
        {
            cells[acceptorCells_[acceptorI]]=1.0;
            acceptorCellSet.insert(acceptorCells_[acceptorI]);
        }

        // Donor cells
        forAll(primaryDonorCells_,donorI)
        {
            cells[primaryDonorCells_[donorI]]=2.0;
            donorCellSet.insert(primaryDonorCells_[donorI]);
        }
        acceptorCellSet.instance() = mesh_.thisDb().time().timeName();
        donorCellSet.instance() = mesh_.thisDb().time().timeName();

        acceptorCellSet.write();
        donorCellSet.write();

        cellState.write();

        const List<interpolationItem>& ownItems = ownInterpolationItems();
        cellSet ownDonorItemSet
        (
            mesh_,
            "proc"+name(Pstream::myProcNo())+"ownInterpolationCells",
            ownItems.size()
        );

        forAll(ownItems,itemI)
        {
            ownDonorItemSet.insert(ownItems[itemI].cellID());
        }
        ownDonorItemSet.write();

        const List< List<interpolationItem> >& neighbourItems =
            neighbourInterpolationItems();
        forAll(neighbourItems, procI)
        {
            const List<interpolationItem>& procItems = neighbourItems[procI];
            cellSet procDonorItemSet
            (
                mesh_,
                "proc"+name(Pstream::myProcNo())+"to"+name(procI)+"InterpolationCells",
                procItems.size()
            );
            forAll(procItems, itemI)
            {
                procDonorItemSet.insert(procItems[itemI].cellID());
            }
            procDonorItemSet.write();
        }

        word filename;
        if(Pstream::parRun())
        {
            filename = "donorFaceLinks_" + name(Pstream::myProcNo()) + ".obj";
            OFstream os(filename);
            label count = 1;
            forAll(donorItems, itemI)
            {
                const searchItem& curItem = donorItems[itemI];
                const point& curItemPosition = curItem.position();
                const label& curItemCell = curItem.seed();
                os  <<"v "<<curItemPosition.x()<<" "
                    <<curItemPosition.y()<<" "
                    <<curItemPosition.z()<<nl
                    <<"v "<<ccs[curItemCell].x()<<" "
                    <<ccs[curItemCell].y()<<" "
                    <<ccs[curItemCell].z()<<nl
                    <<"l "<< count <<  " " << count+1 << endl;
                    count+=2;
            }
        }
        else
        {
            filename = "donorFaceLinks.obj";
            OFstream os(filename);
            label count = 1;
            forAll(primaryDonorCells_, cellI)
            {
                const vector& donor = ccs[primaryDonorCells_[cellI]];
                const vector& acceptor = ccs[acceptorCells_[cellI]];
                os  <<"v "<<donor.x()<<" "
                            <<donor.y()<<" "
                            <<donor.z()<<nl
                    <<"v "<<acceptor.x()<<" "
                            <<acceptor.y()<<" "
                            <<acceptor.z()<<nl
                    <<"l "<< count <<  " " << count+1 << endl;
                    count+=2;
            }
        }
    }
}

