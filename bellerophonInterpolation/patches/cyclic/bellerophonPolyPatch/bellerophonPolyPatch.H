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
    Foam::bellerophonPolyPatch

Description
    Cyclic patch for overlapping grid interpolation (bellerophon)

SourceFiles
    bellerophonPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonPolyPatch_H
#define bellerophonPolyPatch_H

#include "coupledPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class bellerophonPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class bellerophonPolyPatch
:
    public coupledPolyPatch
{

private:

    // Private data

        // In which zone shall we search for a donor
        const word donorZoneName_;

        // Name of the overset zone
        mutable word oversetZoneName_;

        // Zone ID of overset zone
        mutable label oversetZoneID_;

        // Transformations

        /*
         * not to be implemented yet, but could be interesting...


            //- Type of transformation - rotational or translational
            transformType transform_;


            // For rotation

                //- Axis of rotation for rotational cyclics
                vector rotationAxis_;

                //- point on axis of rotation for rotational cyclics
                point rotationCentre_;


            // For translation

                //- Translation vector
                vector separationVector_;

         */


    // Private Member Functions

        /*
        void calcTransforms
        (
            const primitivePatch& half0,
            const pointField& half0Ctrs,
            const vectorField& half0Areas,
            const pointField& half1Ctrs,
            const vectorField& half1Areas
        );
        */

protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers& pBufs, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers& pBufs, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("bellerophon");


    // Constructors

        //- Construct from (base couped patch) components
        bellerophonPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word &patchType,
            const transformType transform = UNKNOWN
        );

        //- Construct from dictionary
        bellerophonPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word &patchType
        );

        //- Construct as copy, resetting the boundary mesh
        bellerophonPolyPatch
        (
            const bellerophonPolyPatch& pp,
            const polyBoundaryMesh& bm
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        bellerophonPolyPatch
        (
            const bellerophonPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct given the original patch and a map
        bellerophonPolyPatch
        (
            const bellerophonPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new bellerophonPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new bellerophonPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new bellerophonPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
        virtual ~bellerophonPolyPatch();


    // Member Functions

        // Access
            //- Geometry is not coupled, only the fields are
            virtual bool coupled() const
            {
                return false;
            }

            //- Does this side own the patch? This is always true, since every
            //  patch is stand-alone
            virtual bool owner() const
            {
                return true;
            }

            //- Does the coupled side own the patch? This is always false, since every
            //  patch is stand-alone
            virtual bool neighbour() const
            {
                return false;
            }

            //- Return cell zone name for donor cells
            const word& donorZoneName() const
            {
                return donorZoneName_;
            }

            //- Return ID of overset zone for this patch
            label oversetZoneID() const;

            //- Return boundary for hole
            virtual autoPtr<triSurface> holeBoundary() const;

            /*
             * transformation is not implemented at the moment

            // Transformations

                //- Type of transform
                inline transformType transform() const;

                //- Axis of rotation for rotational cyclic bellerophon
                inline const vector& rotationAxis() const;

                //- point on axis of rotation for rotational cyclic bellerophon
                inline const point& rotationCentre() const;

                //- Translation vector for translational cyclic bellerophon
                inline const vector& separationVector() const;
                */

        //- Transform patch-based positions from nbr side to this side
        virtual void transformPosition(pointField&) const;

        //- Transform a patch-based position from nbr side to this side
        virtual void transformPosition
        (
            point& l,
            const label faceI
        ) const;

        //- Calculate the patch geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const pointField& thisCtrs,
            const vectorField& thisAreas,
            const pointField& thisCc,
            const pointField& nbrCtrs,
            const vectorField& nbrAreas,
            const pointField& nbrCc
        );

        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder
        (
            PstreamBuffers&,
            const primitivePatch&
        ) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Write the coupledPolyPatch data as a dictionary
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bellerophonPolyPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "bellerophonPolyPatchTemplates.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
