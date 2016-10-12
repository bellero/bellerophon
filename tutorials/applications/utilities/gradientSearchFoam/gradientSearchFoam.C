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
    gradientSearchFoam

Description
    Searches for points from dictionary. Points are moved for every time step to
    check ability of using seed cells.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "gradientSearch.H"
#include "List.H"

// Get the accurate cpu time
#include "time.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    gradientSearch search(mesh);

    // Read control settings

    Info<<"\nReading Dictionary"<<endl;
    IOdictionary gradientSearchDict
    (
        IOobject
        (
            "gradientSearchDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const vector velocity = gradientSearchDict.lookup("velocity");
    List<vector> searchablePoints(0);
    if(Pstream::myProcNo() == 0)
    {
        searchablePoints = List<vector>(gradientSearchDict.lookup("points"));
    }

    List<vector> searchablePoints2 = List<vector>(gradientSearchDict.lookup("points"));

    Info<<"Searching points: "<<searchablePoints<<endl;

    const label nIter =
        ((runTime.endTime()-runTime.startTime())/runTime.deltaT()).value();

    scalar startTime = runTime.elapsedCpuTime();
    for(label iterI = 0; iterI<nIter; iterI++)
    {
        Info<< "Gradient search iteration "<< iterI << endl;

        autoPtr< List< searchItem > > resultPtr =
            search.search(searchablePoints);

        List<searchItem>& result = resultPtr();

        forAll(result, resultI)
        {
            const searchItem& curItem = result[resultI];
            Info<<"Point "<<curItem.cellLabel()<<" from proc "<<curItem.procID()
                <<" in cell "<<curItem.seed()<<endl;
        }
    }
    const scalar timeForGradientSearch = runTime.elapsedCpuTime()-startTime;

    startTime = runTime.elapsedCpuTime();
    for(label iterI = 0; iterI<nIter; iterI++)
    {
        Info<< "Bruteforce search iteration "<< iterI << endl;

        forAll(searchablePoints2,pointI)
        {
            const label cellI = mesh.findCell(searchablePoints2[pointI]);
            Info<<"Point "<<pointI<<" in cell "<<cellI<<endl;
        }
    }

    const scalar timeForLinearSearch = runTime.elapsedCpuTime()-startTime;

    Info<<"\nElapsed CPU Time for "<<nIter<<" searches:"<<nl
        <<"    Gradient search: "<<timeForGradientSearch*1000.0<<nl
        <<"    Linear search:   "<<timeForLinearSearch*1000.0<<endl;

    Info<< "\nEnd.\n" << endl;

    return 0;
}


// ************************************************************************* //
