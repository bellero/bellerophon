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

Class
    Foam::bellerophonPointPatch

Description
    Cyclic bellerophon point patch - place holder only

SourceFiles
    bellerophonPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonPointPatch_H
#define bellerophonPointPatch_H

#include "coupledFacePointPatch.H"
#include "bellerophonPolyPatch.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class bellerophonPointPatch Declaration
\*---------------------------------------------------------------------------*/

class bellerophonPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic AMI patch
        const bellerophonPolyPatch& bellerophonPolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        bellerophonPointPatch(const bellerophonPointPatch&);

        //- Disallow default assignment
        void operator=(const bellerophonPointPatch&);


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName(bellerophonPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        bellerophonPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~bellerophonPointPatch();


    // Member Functions

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return type();
        }

        //- Return the underlying bellerophonPolyPatch
        const bellerophonPolyPatch& bellerophonPatch() const
        {
            return bellerophonPolyPatch_;
        }

/*
 * one nice day...

//- Return face transformation tensor
        const tensorField& forwardT() const
        {
            Info<< "bellerophonPointPatch::forwardT()"
                << " at bellerophonPointPatch.H:164" << nl << endl;

            return bellerophonPolyPatch_.forwardT();
        }

        //- Return neighbour-cell transformation tensor
        const tensorField& reverseT() const
        {
            Info<< "bellerophonPointPatch::reverseT()"
                << " at bellerophonPointPatch.H:173" << nl << endl;

            return bellerophonPolyPatch_.reverseT();
        }
 */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
