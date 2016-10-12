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
    Foam::bellerophonInterface

Description
    An abstract base class for cyclic bellerophon interfaces

SourceFiles
    bellerophonInterface.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonInterface_H
#define bellerophonInterface_H

#include "primitiveFieldsFwd.H"
#include "triSurface.H"
#include "bellerophonInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class bellerophonInterface Declaration
\*---------------------------------------------------------------------------*/

class bellerophonInterface
{

    //- Global index of this interface
    label index_;

    //- Index of the hole interface
    mutable label holeInterface_;

    //- Reduce global index of the interface by one
    void decreaseIndex() { index_--; };

    //- Assign value of the hole interface
    void setHoleInterface(const label i) const { holeInterface_ = i; };

public:

    //- Runtime type information
    TypeName("bellerophonInterface");


    // Constructors

        //- Construct null
        bellerophonInterface();


    //- Destructor
    virtual ~bellerophonInterface();


    // Member Functions

        // Access
        virtual label index() const { return index_; };

        virtual label nFaces() const = 0;

        virtual label donorZoneID() const = 0;

        virtual const word& name() const = 0;

        virtual const fvMesh& mesh() const = 0;

        virtual const labelUList& faceCells() const = 0;

        virtual autoPtr<triSurface> holeBoundary() const = 0;

        virtual label holeInterface() const { return holeInterface_; };

/*
            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return face reverse transformation tensor
            virtual const tensorField& reverseT() const = 0;
 */
    friend void bellerophon::removeInterface(bellerophonInterface&);

    friend label bellerophon::addHole
    (
        const labelList& map,
        const label donorZoneID,
        const label oversetPatchID
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
