/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::wakeVelocity

Group
    grpForcesFunctionObjects

Description
    This function object calculates the wakeVelocity and moments by integrating the
    pressure and skin-friction wakeVelocity over a given list of patches.

    Member function wakeVelocity::write() calculates the wakeVelocity/moments and
    writes the wakeVelocity/moments into the file \<timeDir\>/wakeVelocity.dat and bin
    data (if selected) to the file \<timeDir\>/wakeVelocity_bin.dat

    Example of function object specification:
    \verbatim
    wakeVelocity1
    {
        type        wakeVelocity;
        functionObjectLibs ("libwakeVelocity.so");
        ...
        log         yes;
        patches     (walls);

        binData
        {
            nBin        20;
            direction   (1 0 0);
            cumulative  yes;
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: wakeVelocity       | yes         |
        log          | write force data to standard output | no | no
        patches      | patches included in the wakeVelocity calculation | yes |
        pName        | pressure field name     | no          | p
        UName        | velocity field name     | no          | U
        rhoName      | density field name (see below) | no   | rho
        CofR         | centre of rotation (see below) | no   |
        directForceDensity | force density supplied directly (see below)|no|no
        fDName       | name of force density field (see below) | no | fD
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according o
    \table
        nBin         | number of data bins     | yes         |
        direction    | direction along which bins are defined | yes |
        cumulative   | bin data accumulated with incresing distance | yes |
    \endtable

Note
  - For incompressible cases, set \c rhoName to \c rhoInf.  You will then be
    required to provide a \c rhoInf value corresponding to the free-stream
    constant density.
  - If the force density is supplied directly, set the \c directForceDensity
    flag to 'yes', and supply the force density field using the \c
    fDName entry
  - The centre of rotation (CofR) for moment calculations can either be
    specified by an \c CofR entry, or be taken from origin of the local
    coordinate system.  For example,
    \verbatim
        CofR        (0 0 0);
    \endverbatim
    or
    \verbatim
        coordinateSystem
        {
            origin  (0 0 0);
            e3      (0 0 1);
            e1      (1 0 0);
        }
    \endverbatim

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject
    Foam::forceCoeffs

SourceFiles
    wakeVelocity.C
    IOwakeVelocity.H

\*---------------------------------------------------------------------------*/

#ifndef wakeVelocity_H
#define wakeVelocity_H

#include "functionObjectFile.H"
#include "volFieldsFwd.H"
#include "OFstream.H"
#include "Switch.H"
#include "gradientSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class wakeVelocity Declaration
\*---------------------------------------------------------------------------*/

class wakeVelocity
:
    public functionObjectFile
{
protected:

    // Protected data

        //- Name of this set of wakeVelocity,
        //  Also used as the name of the probes directory.
        word name_;

        //- Database this class is registered to
        const objectRegistry& obr_;

        //- Construction dictionary
        dictionary dict_;

        //- Active flag
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- List of field names to operate on
        word fieldName_;

        //- First point of sample line
        vector min_;

        //- Second point of sample line
        vector max_;

        //- Number of points on sample line
        label nPoints_;

        //- Vector of motion amplitude
        vector amplitude_;

        //- Period of motion
        scalar period_;

        //- Phase lag
        scalar phaseLag_;

        //- Search items for interpolation
        autoPtr< List<searchItem> > sampleItemsPtr_;

        //- Disallow default bitwise copy construct
        wakeVelocity(const wakeVelocity&);

        //- Disallow default bitwise assignment
        void operator=(const wakeVelocity&);


public:

    //- Runtime type information
    TypeName("wakeVelocity");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        wakeVelocity
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


    //- Destructor
    virtual ~wakeVelocity();


    // Function object member Functions

        //- Read the wakeVelocity data
        virtual void read(const dictionary&);

        //- Write to screen/file
        virtual void write();

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
