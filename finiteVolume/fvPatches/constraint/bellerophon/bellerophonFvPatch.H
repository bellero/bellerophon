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
    Foam::bellerophonFvPatch

Description
    Cyclic coupled patch for Overlapping Grid Interface (bellerophon)

SourceFiles
    bellerophonFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonFvPatch_H
#define bellerophonFvPatch_H
#include "coupledFvPatch.H"
#include "bellerophonInterface.H"
#include "bellerophonPolyPatch.H"
#include "fvBoundaryMesh.H"
#include "bellerophonLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class bellerophonFvPatch Declaration
\*---------------------------------------------------------------------------*/

class bellerophonFvPatch
:
    public coupledFvPatch,
    public bellerophonInterface
{
    // Private data

        // Reference to the poly patch
        const bellerophonPolyPatch& bellerophonPolyPatch_;

        // Zone to search donors in
        mutable label donorZoneID_;

protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

public:

    //- Runtime type information
    TypeName(bellerophonPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        bellerophonFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            bellerophonInterface(),
            bellerophonPolyPatch_(refCast<const bellerophonPolyPatch>(patch)),
            donorZoneID_(-1)
        {}


    // Member functions

        // Access

            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return coupledFvPatch::faceCells();
            }

            //- Return number of faces on the patch (Reimplemented from
            // bellerophonInterface)
            virtual label nFaces() const { return size(); };

            //- Return ID of the zone to search donors in (Reimplemented from
            // bellerophonInterface)
            virtual label donorZoneID() const
            {
                if(donorZoneID_ == -1)
                {
                    donorZoneID_ = this->boundaryMesh().mesh().cellZones()
                        .findZoneID(bellerophonPolyPatch_.donorZoneName());
                }
                return donorZoneID_;
            };

            virtual const word& name() const { return coupledFvPatch::name(); };

            virtual autoPtr<triSurface> holeBoundary() const
            {
                return bellerophonPolyPatch_.holeBoundary();
            };

            virtual const fvMesh& mesh() const
            {
                return coupledFvPatch::boundaryMesh().mesh();
            };

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Initialise "neighbour" field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
             labelUList& iF
            ) const
            {/*do nothing*/}

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
             const labelUList& iF
            ) const;

/*
 * TODO:
 * no need for this, right?

        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
 */

            //- Are the planes parallel
            virtual bool parallel() const
            {
                Info<<"bellerophonFvPatch::parallel(..)"<<nl<<endl;

                return bellerophonPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                Info<<"bellerophonFvPatch::forwardT(..)"<<nl<<endl;

                return bellerophonPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor.
            virtual const tensorField& reverseT() const
            {
                Info<<"bellerophonFvPatch::reverseT(..)"<<nl<<endl;

                return bellerophonPolyPatch_.reverseT();
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
