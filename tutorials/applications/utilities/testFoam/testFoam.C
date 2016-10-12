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
#include "fvPatchField.H"
#include "bellerophonFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "Hello world, time is "<<runTime.timeName()<<"."<<nl<<endl;

    volVectorField cellC
    (
        IOobject
        (
            "cellC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0",dimLength,vector::zero)
    );

    Info<<"1"<<endl;

//     forAll(C.boundaryField(),patchI)
//     {
//         scalarField y1 = C.boundaryField()[patchI].component(1);
//         scalarField y2 = mesh.boundary()[patchI].Cf().component(1);
//         Info<<"Delta: "<< (y1-y2) <<endl;
//     }

    Info<<"2"<<endl;

    cellC = mesh.C();

    Info<<"3"<<endl;

    forAll(cellC.boundaryField(), patchI)
    {
        if(isA<bellerophonFvPatchField<vector> >(cellC.boundaryField()[patchI]))
        {
            cellC.boundaryField()[patchI].initEvaluate();
        }
    }

    Info<<"4"<<endl;

    cellC.write();

    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
