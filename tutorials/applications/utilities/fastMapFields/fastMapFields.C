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
    fastMapFields

Description
    mapFields implementation using gradient search for interpolation cells;
    map method is direct mapping (value of the cell containing the cell centre)
    currently only supports mapping to internal cells, patches are not changed

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "gradientSearch.H"
#include "IOobjectList.H"
#include <../db/IOstreams/Pstreams/PstreamBuffers.H>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void mapVolField
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const word& name,
    const List<searchItem>& items,
    const labelList& nItemsPerProc,
    const labelListList& itemToCell
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    IOobjectList objectsSource(meshSource, meshSource.time().timeName());
    IOobjectList objectsTarget(meshTarget, meshTarget.time().timeName());

    Info << "    interpolating field "<<name<<endl;
    const fieldType fieldSource
    (
        IOobject
        (
            name,
            meshSource.time().timeName(),
            meshSource,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        meshSource
    );
    fieldType fieldTarget
    (
        IOobject
        (
            name,
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshTarget
    );


    const Field<Type>& ifSource = fieldSource.internalField();
    Field<Type>& ifTarget = fieldTarget.internalField();

    List<List < Type > > transfer(nItemsPerProc.size());

    labelList itemPerProcI(nItemsPerProc.size(),0);

    forAll(nItemsPerProc, procI)
    {
        transfer[procI].setSize(nItemsPerProc[procI]);
    }

    forAll(items, itemI)
    {
        const searchItem& curItem = items[itemI];
        const label procI = curItem.procID();
        transfer[procI][itemPerProcI[procI]++] = ifSource[curItem.seed()];
    }

    PstreamBuffers pBufs(Pstream::nonBlocking);

    for(label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if(procI != Pstream::myProcNo())
        {
            UOPstream toDomain(procI, pBufs);
            toDomain << transfer[procI];
        }
    }

    pBufs.finishedSends();

    for(label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        const labelList& itc = itemToCell[procI];
        Field<Type> fromDomain;

        if(procI != Pstream::myProcNo())
        {
            UIPstream str(procI, pBufs);
            str >> fromDomain;
        }
        else
        {
            fromDomain = transfer[procI];
        }

        forAll(fromDomain, valueI)
        {
            ifTarget[itc[valueI]] = fromDomain[valueI];
        }
    }

    fieldTarget.write();
}

#define mapTypeField(Type)                                                     \
    {                                                                          \
        IOobjectList typedObjectsSource =                                      \
            objectsSource.lookupClass(GeometricField<Type, fvPatchField, volMesh>::typeName);\
        IOobjectList typedObjectsTarget =                                      \
            objectsTarget.lookupClass(GeometricField<Type, fvPatchField, volMesh>::typeName);\
                                                                               \
        forAllIter(IOobjectList, typedObjectsSource, objectIter)               \
        {                                                                      \
            const word fieldName = objectIter()->name();                       \
            if                                                                 \
            (                                                                  \
                (                                                              \
                    selectedFields.found(fieldName)                            \
                    ||                                                         \
                    selectedFields.empty()                                     \
                )                                                              \
                &&                                                             \
                typedObjectsTarget.found(fieldName)                            \
            )                                                                  \
            {                                                                  \
                mapVolField<Type>                                              \
                (                                                              \
                    meshSource,                                                \
                    meshTarget,                                                \
                    fieldName,                                                 \
                    items,                                                     \
                    nItemsPerProc,                                             \
                    itemsToCell                                                \
                );                                                             \
            }                                                                  \
        }                                                                      \
    }

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "map internal fields from one mesh to another, cell  search using "
        "gradient search, direct mapping without interpolation"
    );

    argList::validArgs.append("sourceCase");

    argList::addOption
    (
        "sourceTime",
        "scalar|'latestTime'",
        "specify the source time"
    );
    argList::addOption
    (
        "sourceRegion",
        "word",
        "specify the source region"
    );
    argList::addOption
    (
        "targetRegion",
        "word",
        "specify the target region"
    );
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be mapped. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );

    argList::addBoolOption
    (
        "useZones",
        "use cell zones (labels have to coincide)"
    );

    argList args(argc, argv);

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    const fileName casePath = args[1];
    const fileName rootDirSource = casePath.path();
    const fileName caseDirSource = casePath.name();

    Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;
    word sourceRegion = fvMesh::defaultRegion;
    if (args.optionFound("sourceRegion"))
    {
        sourceRegion = args["sourceRegion"];
        Info<< "Source region: " << sourceRegion << endl;
    }

    Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;
    word targetRegion = fvMesh::defaultRegion;
    if (args.optionFound("targetRegion"))
    {
        targetRegion = args["targetRegion"];
        Info<< "Target region: " << targetRegion << endl;
    }

    Info<< "\nCreate databases as time" << endl;

    HashTable<string> srcOptions(args.options());
    srcOptions.erase("case");
    srcOptions.insert("case", fileName(rootDirSource/caseDirSource));

    argList argsSrc(args, srcOptions, false, false, false);

    Time runTimeSource(Time::controlDictName, argsSrc);

    Time runTimeTarget(Time::controlDictName, args);

    HashTable<word> patchMap;
    wordList cuttingPatches;

    const bool useZones = args.optionFound("useZones");

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    {
        instantList sourceTimes = runTimeSource.times();
        label sourceTimeIndex = runTimeSource.timeIndex();
        if (args.optionFound("sourceTime"))
        {
            if (args["sourceTime"] == "latestTime")
            {
                sourceTimeIndex = sourceTimes.size() - 1;
            }
            else
            {
                sourceTimeIndex = Time::findClosestTimeIndex
                (
                    sourceTimes,
                    args.optionRead<scalar>("sourceTime")
                );
            }
        }
        else
        {
            sourceTimeIndex = Time::findClosestTimeIndex
            (
                sourceTimes,
                runTimeTarget.time().value()
            );
        }

        runTimeSource.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);

        Info<< "\nSource time: " << runTimeSource.value()
            << "\nTarget time: " << runTimeTarget.value()
            << endl;
    }

    Info<< "\nCreate meshes\n" << endl;

    fvMesh meshSource
    (
        IOobject
        (
            sourceRegion,
            runTimeSource.timeName(),
            runTimeSource
        )
    );

    fvMesh meshTarget
    (
        IOobject
        (
            targetRegion,
            runTimeTarget.timeName(),
            runTimeTarget
        )
    );

    Info<< "Source mesh size: " << meshSource.nCells() << tab
        << "Target mesh size: " << meshTarget.nCells() << nl << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    autoPtr< List< searchItem > > itemsPtr(new List< searchItem >(meshTarget.nCells()));
    List<searchItem>& items = itemsPtr();

    const vectorField& cc = meshTarget.C().internalField();

    forAll(cc, cellI)
    {
        searchItem& item = items[cellI];
        item.cellLabel() = cellI;
        item.procID() = Pstream::myProcNo();
        item.position() = cc[cellI];
        item.seed() = 0;
        item.zoneID() = useZones ? meshTarget.cellZones().whichZone(cellI) : -1;
        item.groupID() = -1;
    }

    Info << "Starting gradient search for "<<cc.size()<<" items."<<nl<<endl;

    gradientSearch gs(meshSource);
    gs.search(itemsPtr);

    autoPtr< List< searchItem> > failedItemsPtr = gs.failItems();

    Info << "Finished search."<<nl<<endl;

    // Number of search items per processor
    labelList nItemsPerProc(Pstream::nProcs(), 0);
    labelListList itemsToCell(Pstream::nProcs(), labelList(0));

    // Count items per proc
    forAll(items, itemI)
    {
        nItemsPerProc[items[itemI].procID()]++;
    }

    Info <<"Items per proc: "<<nItemsPerProc<<endl;

    // Allocate
    forAll(itemsToCell, procI)
    {
        itemsToCell[procI].setSize(nItemsPerProc[procI]);
        nItemsPerProc[procI] = 0;
        Info<<"Size "<<procI<<": "<<itemsToCell.size()<<endl;
    }

    // Sort cell label for items
    forAll(items, itemI)
    {
        const label procI = items[itemI].procID();
        itemsToCell[procI][nItemsPerProc[procI]++] =
            items[itemI].cellLabel();
    }

    // Communicate
    PstreamBuffers pBufs(Pstream::nonBlocking);

    for(label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if(procI != Pstream::myProcNo())
        {
            UOPstream toDomain(procI, pBufs);
            toDomain << itemsToCell[procI];
        }
    }

    pBufs.finishedSends();

    for(label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if(procI != Pstream::myProcNo())
        {
            UIPstream str(procI, pBufs);
            str >> itemsToCell[procI];
        }
    }

    IOobjectList objectsSource(meshSource, meshSource.time().timeName());
    IOobjectList objectsTarget(meshTarget, meshTarget.time().timeName());

    mapTypeField(scalar);
    mapTypeField(vector);
    mapTypeField(sphericalTensor);
    mapTypeField(symmTensor);
    mapTypeField(tensor);

    Info<< "\nEnd.\n" << endl;

    return 0;
}
