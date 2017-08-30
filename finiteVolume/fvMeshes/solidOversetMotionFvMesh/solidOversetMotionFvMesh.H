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
    Foam::solidOversetMotionFvMesh

Description
    Foam::solidOversetMotionFvMesh

SourceFiles
    solidOversetMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef solidOversetMotionFvMesh_H
#define solidOversetMotionFvMesh_H

#include "bellerophonFvMesh.H"

#include "solidBodyMotionFunction.H"
#include "triSurface.H"
#include "gradientSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class solidOversetMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class solidOversetMotionFvMesh
:
    public bellerophonFvMesh
{
private:
        // Private Member Functions

        //- Disallow default bitwise copy construct
        solidOversetMotionFvMesh(const solidOversetMotionFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const solidOversetMotionFvMesh&);


protected:
    // Protected data
        //- Dictionary
        const dictionary dict_;

        //- Motion function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- Transformation from the previous step
        septernion oldTransf_;

        //- Untransformed points
        pointIOField points0_;

        //- ID of overset patch
        const label oversetPatchID_;

        //- ID of donor zone
        label oversetZoneID_;

        //- Name of velocity field to correct BCs
        const word UName_;

        //- Check live cells for being hole cells in previous step
        const Switch checkLiveCells_;

        //- Interpolate live cells that were hole cells in previous step
        const Switch interpolateLiveCells_;

        //- Name of the interpolation scheme to use for new live cells
        const word interpolationSchemeName_;

        //- Correct velocity on surfaces after mesh motion
        const Switch correctUf_;

        //- Name of velocity field on surfaces
        const word UfName_;

        //- Name of face fluxes
        const word phiName_;

        //- Index of the interface of the hole boundary
        label interfaceIndex_;

        //- Hole boundary surface
        autoPtr<triSurface> holeBoundaryPtr_;

        //- Untransformed boundary points
        pointField boundaryPoints0_;

        //- Labels of moving points
        labelList movingPointIDs_;

        //- Labels of moving points
        labelList movingCellIDs_;

    // Protected member functions
        //- Mark moving cells and points
        void generateMotionMap();

        //- Add hole cells to interpolation
        label addHole();

        //- Interpolate values on new live cells
        void interpolateNewLiveCells
        (
            const boolList& map,
            const pointField& transformedPoss
        );

        //- Interpolate fields by given search items and id to proc map
        template<class Type>
        void interpolateFields
        (
            const List<searchItem>& items,
            const labelListList& IDtoProc
        );

public:

    //- Runtime type information
    TypeName("solidOversetMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        solidOversetMotionFvMesh(const IOobject& io);


    //- Destructor
    ~solidOversetMotionFvMesh();


    // Member Functions
        // Generate hole map from triangular surface
        autoPtr<labelList> createHoleMap() const;

        //- Update the mesh
        virtual bool update();

        bool writeLiveCells(const word&) const;

        //- Writing the mesh
        virtual bool writeObjects
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
