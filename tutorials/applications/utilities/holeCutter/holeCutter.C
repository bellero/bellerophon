/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    parallelTestFoam

Description
    Shows how to send information to specific CPUs instead of reconstructing
    global field.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "holeSource.H"

#include "bellerophonInterface.H"
#include "gradientSearch.H"
#include "cellSet.H"
#include "SLList.H"
#include "argList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

enum cellState
{
    REGULAR_CELL = 0,            // Regular fluid cell, value from equations
    DONOR_CELL = 1,              // Donor cell, value from equations
    EXTENDED_ACCEPTOR_CELL = 2,  // Extended acceptor cell, next to acceptor,
                                 // may not become primary donor cell as it
                                 // is neighbour to acceptor cell
    ACCEPTOR_CELL = 3,           // Acceptor cell, value is interpolated
    POTENTIAL_ACCEPTOR_CELL = 4, // Potential acceptor, needs to be verified
    HOLE_CELL = 5                // Hole cell, not involved in solution
};

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "cellState",
        "write cell state as volScalarField for post-processing."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<<"Reading dictionary."<<nl<<endl;

    IOdictionary dict
    (
        IOobject
            (
                "holeCutterDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ
            )
    );

    autoPtr<holeSource> hsPtr = holeSource::New(mesh, dict);
    holeSource& hs = hsPtr();

    // ID of the overset zone (no hole cells here)
    const label oversetZoneID =
        mesh.cellZones().findZoneID(dict.lookup("oversetZone"));

    if(oversetZoneID < 0)
    {
        FatalErrorIn("holeCutter.C")
        << "Unknown cell zone "<<dict.lookup("oversetZone")
        << exit(FatalError);
    }

    // Cells in the overset zone;
    const labelList& oversetCells = mesh.cellZones()[oversetZoneID];

    // Centres of cells
    const vectorField& centres = mesh.cellCentres();

    // State of cells
    labelList cellState(centres.size(), REGULAR_CELL);

    // Number of hole cells returned by marking
    label nHoleCells = hs.markCells(cellState, HOLE_CELL);

    // Removing overset cells from set
    forAll(oversetCells,cellI)
    {
        if(cellState[oversetCells[cellI]] != REGULAR_CELL)
        {
            cellState[oversetCells[cellI]] = REGULAR_CELL;
            nHoleCells--;
        }
    }

    Info<<"Marked "<<nHoleCells<<" cells inside hole geometry."<<nl<<endl;

    const labelList& acceptorCells =
        bellerophon::Interpolation().acceptorCells();

    forAll(acceptorCells, acceptorI)
    {
        label& curState = cellState[acceptorCells[acceptorI]];
        if(curState == HOLE_CELL)
        {
            nHoleCells--;
        }
        curState = ACCEPTOR_CELL;
    }

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    // Mark neighbours of acceptor cells
    forAll(own, faceI)
    {
        const label owner = own[faceI];
        const label neighbour = nei[faceI];

        if
        (
            (cellState[owner] == ACCEPTOR_CELL)
            &&
            (cellState[neighbour] == REGULAR_CELL)
        )
        {
            cellState[neighbour] = EXTENDED_ACCEPTOR_CELL;
        }
        else if
        (
            (cellState[owner] == REGULAR_CELL)
            &&
            (cellState[neighbour] == ACCEPTOR_CELL)
        )
        {
            cellState[owner] = EXTENDED_ACCEPTOR_CELL;
        }
    }
    const labelList& primaryDonorCells =
        bellerophon::Interpolation().primaryDonorCells();

    // Mark primary donor cells
    forAll(primaryDonorCells, donorI)
    {
        label& curState = cellState[primaryDonorCells[donorI]];
        if(curState == HOLE_CELL)
        {
            nHoleCells--;
        }
        curState = DONOR_CELL;
    }

    const List<interpolationItem>& ownItems =
        bellerophon::Interpolation().ownInterpolationItems();

    forAll(ownItems, itemI)
    {
        label& curState = cellState[ownItems[itemI].cellID()];
        if(curState == HOLE_CELL)
        {
            nHoleCells--;
        }
        curState = DONOR_CELL;
    }

    const List<List <interpolationItem> >& neighbourItems =
        bellerophon::Interpolation().neighbourInterpolationItems();

    forAll(neighbourItems, procI)
    {
        const List<interpolationItem>& procItems = neighbourItems[procI];
        forAll(procItems, itemI)
        {
            label& curState = cellState[procItems[itemI].cellID()];
            if(curState == HOLE_CELL)
            {
                nHoleCells--;
            }
            curState = DONOR_CELL;
        }
    }

    forAll( bellerophon::Interpolation().neighbourInterpolationItems(), procI)
    {
        const List<interpolationItem>& neighbourItems =
            bellerophon::Interpolation().neighbourInterpolationItems()[procI];

        forAll(neighbourItems, itemI)
        {
            label& curState = cellState[neighbourItems[itemI].cellID()];
            if(curState == HOLE_CELL)
            {
                nHoleCells--;
            }
            curState = DONOR_CELL;
        }
    }

    // Create List of hole cells and "shrink" the hole, so all the first row of
    Info<<"Marked "<<nHoleCells
        <<" cells after removing donor and acceptor cells."<<nl<<endl;

    if(args.optionFound("cellState"))
    {
        Info<<"Writing initial cell state field."<<nl<<endl;
        volScalarField cellStateField
        (
            IOobject
            (
                "initialCellState",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("null",dimless,0.0)
        );

        scalarField& cellStateIField = cellStateField.primitiveFieldRef();
        forAll(cellStateIField,cellI)
        {
            cellStateIField[cellI] = cellState[cellI];
        }
        cellStateField.write();
    }

    // List of potential acceptor cells
    SLList<searchItem> potentialAcceptorItems;

    // Find potential interface faces:
    // Check all faces: is owner XOR neighbour hole cell?
    //  -> Search potential donor
    //  -> If found and donor is regular cell -> hole cell remains hole cell
    //     else neighbours of face cannot stay hole cell, becomes regular cell

    const vectorField& cc = mesh.cellCentres();

    const labelList& zoneSeeds = bellerophon::Interpolation().zoneSeeds();


    nHoleCells = 0;
    forAll(cellState,cellI)
    {
        if(cellState[cellI] == HOLE_CELL)
        {
            nHoleCells++;
        }
    }

    bool doSearch = true;
    while(doSearch)
    {
        // Reset potential acceptor items
        potentialAcceptorItems.clear();

        forAll(nei, faceI)
        {
            const label owner = own[faceI];
            const label neighbour = nei[faceI];

            label& ownerState = cellState[owner];
            label& neighbourState = cellState[neighbour];

            label potentialAcceptor = -1;
            if
            (
                (ownerState == HOLE_CELL)
                &&
                (neighbourState < ACCEPTOR_CELL)
            )
            {
                potentialAcceptor = owner;
                ownerState = POTENTIAL_ACCEPTOR_CELL;
            }
            else if
            (
                (neighbourState == HOLE_CELL)
                &&
                (ownerState < ACCEPTOR_CELL)
            )
            {
                potentialAcceptor = neighbour;
                neighbourState = POTENTIAL_ACCEPTOR_CELL;
            }

            if(potentialAcceptor != -1)
            {
                potentialAcceptorItems.append
                (
                    searchItem
                    (
                        Pstream::myProcNo(),
                        potentialAcceptor,
                        cc[potentialAcceptor],
                        zoneSeeds[oversetZoneID],
                        oversetZoneID
                    )
                );

                nHoleCells--;
            }
        }
        // TODO: ADD PARALLEL STUFF HERE

        Info<<"Marked "<<nHoleCells<<" after shrinking hole and "<<potentialAcceptorItems.size()<<" potential acceptors."<<nl<<endl;

        // Possible acceptor cells as autoPtr< List<searchItem> > for gradient
        // search

        autoPtr< List<searchItem> > possibleAcceptorSearchItemsPtr;

        // This is ugly but directly creating a copy causes an error in of-dev - AG 2016/11/30
        possibleAcceptorSearchItemsPtr.set
        (
            new List<searchItem>(potentialAcceptorItems.size())
        );
        List<searchItem>& possibleAcceptorSearchItems = possibleAcceptorSearchItemsPtr();
        label itemI = 0;
        forAllConstIter(SLList<searchItem>, potentialAcceptorItems, iter)
        {
            possibleAcceptorSearchItems[itemI++] = *iter;;
        }
        

        // Search for donor cells of possible acceptor cells
        gradientSearch gs(mesh);
        gs.search(possibleAcceptorSearchItemsPtr);

        Info<<"Finished search"<<endl;

        List<searchItem>& foundAcceptorCells = possibleAcceptorSearchItemsPtr();
        forAll(foundAcceptorCells,itemI)
        {
            searchItem& curItem = foundAcceptorCells[itemI];
            const label potentialAcceptor = curItem.cellLabel();

            label& donorState = cellState[curItem.seed()];
            label& acceptorState = cellState[potentialAcceptor];

            // TODO: communicate success/fail back to originating proc of search
            //       item

            if(donorState == REGULAR_CELL)
            {
                // Found a regular cell as donor cell, mark as acceptor
                donorState = DONOR_CELL;
                acceptorState = ACCEPTOR_CELL;

                Info<<"Found donor for acceptor cell."<<endl;
            }
            else if(donorState == DONOR_CELL)
            {
                // Donor cell already is donor itself
                acceptorState = ACCEPTOR_CELL;
            }
            else if(donorState == ACCEPTOR_CELL)
            {
                // Donor cell for potential acceptor cell is already acceptor
                // cell; potential acceptor will stay regular cell and the
                // hole neighbour of the face will become a regular cell.
                Info<<"Acceptor cell cannot be a donor cell for other cells"
                    <<endl;
                acceptorState = REGULAR_CELL;
            }
            else if(donorState == EXTENDED_ACCEPTOR_CELL)
            {
                // Donor cell for potential acceptor cell is already acceptor
                // cell; potential acceptor will stay regular cell and the
                // hole neighbour of the face will become a regular cell.
                Info<<"Extended acceptor cell cannot be a donor cell for other "
                    <<"cells"<<endl;
                acceptorState = REGULAR_CELL;
            }
            else if(donorState == POTENTIAL_ACCEPTOR_CELL)
            {
                FatalError<<"Potential acceptor at "<<curItem.position()
                    <<" should not happen to be "
                    <<"identified as donor cell."<<nl<<"Check cell zones."
                    <<abort(FatalError);
            }
            else if(donorState == HOLE_CELL)
            {
                FatalError<<"Hole cell should not happen to be "
                    <<"identified as donor cell."<<nl<<"Check cell zones."
                    <<abort(FatalError);
            }
            else
            {
                FatalError<<"Unknown cell state!"<<abort(FatalError);

            }
        }

        doSearch = potentialAcceptorItems.size() > 0;

        reduce(doSearch, orOp<bool>());
    }

    nHoleCells = 0;
    forAll(cellState,cellI)
    {
        if(cellState[cellI] == HOLE_CELL)
        {
            nHoleCells++;
        }
    }

    Info<<"Found "<<nHoleCells<<" hole cells."<<endl;

    labelList holeCells(nHoleCells,-1);
    nHoleCells = 0;
    forAll(cellState,cellI)
    {
        if(cellState[cellI] == HOLE_CELL)
        {
            holeCells[nHoleCells++] = cellI;
        }
    }


    cellSet holeCellSet(runTime,"holeCells",holeCells);
    cellSet liveCellSet(runTime,"liveCells",holeCellSet);
    liveCellSet.invert(cellState.size());

    holeCellSet.write();
    liveCellSet.write();

    if(args.optionFound("cellState"))
    {
        Info<<"Writing cell state field."<<nl<<endl;
        volScalarField cellStateField
        (
            IOobject
            (
                "cellState",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("null",dimless,0.0)
        );

        scalarField& cellStateIField = cellStateField.primitiveFieldRef();
        forAll(cellStateIField,cellI)
        {
            cellStateIField[cellI] = cellState[cellI];
        }
        cellStateField.write();
    }
    Info<<"End.\n"<<endl;

    return 0;
}


// ************************************************************************* //
