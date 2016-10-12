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
#include "mapDistribute.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const label myProcNo = Pstream::myProcNo();
    const label nProcs = Pstream::nProcs();
    const label nLoops = 50;

    // Generate some data to send
    labelList labels(100000);
    forAll(labels,labelI)
    {
        labels[labelI]=myProcNo;
    }

    labelListList receivedLabels(nProcs);
/*
    labelListList receivedLabels2(nProcs);
    Info<<"Starting sending using operator<<()"<<endl;
    std::clock_t a = std::clock();
*/
    for(label i = 0; i < nLoops; i++)
    {
        PstreamBuffers pBufs(Pstream::nonBlocking);

        for (label domain = 0; domain < nProcs; domain++)
        {
            if (domain != myProcNo)
            {
                UOPstream toDomain(domain, pBufs);
                toDomain << labels;
            }
        }

        pBufs.finishedSends();

        receivedLabels[myProcNo] = labels;

        // Consume
        for (label domain = 0; domain < nProcs; domain++)
        {
            if (domain != myProcNo)
            {
                UIPstream str(domain, pBufs);

                str >> receivedLabels[domain];
            }
        }

    }

    Info << "Received: " << receivedLabels << endl;
/*    std::clock_t b = std::clock();
    Info << "Finished sending using operator<<()\nDifference: "
         << (b - a) << endl;

    Info<<"Starting sending as char*"<<endl;
    a = std::clock();


    for(label i = 0; i < nLoops; i++)
    {
        const label nReq = UPstream::nRequests();
        const label tag = UPstream::msgType();

        // Writing
        for (label procI = 0; procI < nProcs; procI++)
        {
            if (procI != myProcNo)
            {
                const label nBytes = labels.byteSize();

                receivedLabels2[procI].setSize((procI+10+nProcs)*10);

                IPstream::read
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<char*>(receivedLabels2[procI].begin()),
                    receivedLabels2[procI].byteSize(),
                    tag,
                    0
                );

                OPstream::write
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<const char*>(labels.begin()),
                    nBytes,
                    tag,
                    0
                );
            }
        }

        receivedLabels2[myProcNo] = labels;
        Pstream::waitRequests(nReq);
    }

    b = std::clock();
    Info << "Finished sending as char*\nDifference: "
         << (b - a) << endl;
*/
    Info<< "\nEnd.\n" << endl;

    return 0;
}


// ************************************************************************* //
