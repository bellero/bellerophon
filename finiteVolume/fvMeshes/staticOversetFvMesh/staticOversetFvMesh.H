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
    Foam::staticOversetFvMesh

Description
    Foam::staticOversetFvMesh

SourceFiles
    staticOversetFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef staticOversetFvMesh_H
#define staticOversetFvMesh_H


#include "solidOversetMotionFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class staticOversetFvMesh Declaration
\*---------------------------------------------------------------------------*/

class staticOversetFvMesh
:
public solidOversetMotionFvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        staticOversetFvMesh(const staticOversetFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const staticOversetFvMesh&);

public:

    //- Runtime type information
    TypeName("staticOversetFvMesh");


    // Constructors

        //- Construct from IOobject
        staticOversetFvMesh(const IOobject& io);


    //- Destructor
    ~staticOversetFvMesh();


        //- Update the mesh
        virtual bool update();

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
