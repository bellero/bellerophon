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

\*---------------------------------------------------------------------------*/

#include "bellerophonFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bellerophonFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, bellerophonFvPatch, polyPatch);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<labelField> bellerophonFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    Info<<"interfaceInternalField"<<endl;

    return patchInternalField(internalData);
}



tmp<labelField> bellerophonFvPatch::internalFieldTransfer
(
    const UPstream::commsTypes commsType,
    const labelUList& iF
) const
{
    /*
     * TODO: interpolate field on "neighbour"?
     */
    Info<<"internalFieldTransfer"<<endl;

    return patchInternalField(iF);
}


void Foam::bellerophonFvPatch::makeWeights(scalarField& w) const
{
    w=0.5;
}

Foam::tmp<Foam::vectorField> Foam::bellerophonFvPatch::delta() const
{
//     return 2.0*(Cf()-Cn());
    return Cf()-Cn();
}


/*
Foam::tmp<Foam::labelField> Foam::bellerophonFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    Info<< "bellerophonFvPatch::interfaceInternalField(...)"
        << " at bellerophonFvPatch.C:118" << nl << endl;eval

    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::bellerophonFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    //return neighbFvPatch().patchInternalField(iF);

    FatalErrorIn("bellerophonFvPatch::internalFieldTransfer(...)")
    <<"internalFieldTransfer(...) not implemented yet for bellerophonFvPatch"
    <<" on Patch " << name()
    <<exit(FatalError);
}
 */

} // End namespace Foam

// ************************************************************************* //
