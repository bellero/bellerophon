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
    Foam::bellerophonGAMGInterface

Description
    GAMG agglomerated cyclic bellerophon interface.

SourceFiles
    bellerophonGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef bellerophonGAMGInterface_H
#define bellerophonGAMGInterface_H

#include "GAMGInterface.H"
#include "bellerophonInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class bellerophonGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class bellerophonGAMGInterface
:
    public GAMGInterface,
    virtual public bellerophonInterface
{
    // Private data

        //- Reference for the cyclicLduInterface from which this is
        //  agglomerated
        const bellerophonInterface& fineBellerophonInterface_;

        //- bellerophon interface
        autoPtr<bellerophonFieldToPatchInterpolation> bellerophonPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        bellerophonGAMGInterface(const bellerophonGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const bellerophonGAMGInterface&);


public:

    //- Runtime type information
    TypeName("bellerophon");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        bellerophonGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& restrictAddressing,
            const labelField& neighbourRestrictAddressing
        );


    //- Destructor
    virtual ~bellerophonGAMGInterface();


    // Member Functions

        // Interface transfer functions

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;


//         //- Cyclic interface functions
// 
//             //- Return neigbour processor number
//             virtual label neighbPatchID() const
//             {
//                 Info<< "bellerophonGAMGInterface::neighbPatchID()"
//                     << " at bellerophonGAMGInterface.H:116" << nl << endl;
// 
//                 return fineBellerophonInterface_.neighbPatchID();
//             }
// 
//             virtual bool owner() const
//             {
//                 Info<< "bellerophonGAMGInterface::owner()"
//                     << " at bellerophonGAMGInterface.H:121" << nl << endl;
// 
//                 return fineBellerophonInterface_.owner();
//             }
// 
//             virtual const bellerophonGAMGInterface& neighbPatch() const
//             {
//                 FatalErrorIn("bellerophonGAMGInterface::neighbPatch()")
//                 << "neighbPatch not implemnted for patch " << name()
//                 << exit(FatalError);
//             }

            virtual const bellerophonFieldToPatchInterpolation& bellerophon() const
            {
                Info<< "bellerophonGAMGInterface::bellerophon()"
                    << " at bellerophonGAMGInterface.H:140" << nl << endl;

                return bellerophonPtr_();
            }

//             //- Return face transformation tensor
//             virtual const tensorField& forwardT() const
//             {
//                 Info<< "bellerophonGAMGInterface::forwardT()"
//                     << " at bellerophonGAMGInterface.H:149" << nl << endl;
// 
//                 return fineBellerophonInterface_.forwardT();
//             }
// 
//             //- Return neighbour-cell transformation tensor
//             virtual const tensorField& reverseT() const
//             {
//                 Info<< "bellerophonGAMGInterface::reverseT()"
//                     << " at bellerophonGAMGInterface.H:158" << nl << endl;
// 
//                 return fineBellerophonInterface_.reverseT();
//             }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
