/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "wakeVelocity.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "mathematicalConstants.H"
#include "PstreamCombineReduceOps.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wakeVelocity, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wakeVelocity::wakeVelocity
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    functionObjectFile(obr, typeName, name),
    name_(name),
    obr_(obr),
    dict_(dict),
    active_(true),
    log_(true),
    fieldName_(word::null),
    min_(vector::zero),
    max_(vector::zero),
    nPoints_(0),
    amplitude_(vector::zero),
    period_(0),
    phaseLag_(0)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (isA<fvMesh>(obr_))
    {
        if (readFields)
        {
            read(dict);
            Info<< endl;
        }
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "Foam::wakeVelocity::wakeVelocity"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wakeVelocity::~wakeVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wakeVelocity::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", true);

        if(log_) Info<< type() << " " << name_ << ":" << nl;

        min_      = dict.lookup("min");
        max_      = dict.lookup("max");
        nPoints_  = label(readScalar(dict.lookup("nPoints")));
        amplitude_= dict.lookup("amplitude");
        period_   = readScalar(dict.lookup("period"));
        phaseLag_ = readScalar(dict.lookup("phaseLag"));

        if(nPoints_ < 2)
        {
            Info<<"    deactivating "<<name_<<" because number of points is "
                <<"to small."<<endl;
            active_ = false;
            return;
        }

        fieldName_ = dict.lookupOrDefault<word>("UName", "U");

        if(log_) Info<<"    Sampling field "<<fieldName_<<"."<<nl
            <<"    Using "<<nPoints_<<" points between initial positions "
            <<min_<<" and "<<max_<<"."<<nl
            <<"    Point offset: "<<amplitude_<<" * sin( 2*pi / "<<period_
            <<" * t + "<<phaseLag_<<" / 180 * pi )"<<endl;

        if(Pstream::master())
        {
            const scalar& pi = Foam::constant::mathematical::pi;
            const scalar& t = obr_.time().value();
            const vector offset = amplitude_ *
                sin(2.0 * pi / period_ * t + pi / 180.0 * phaseLag_);

            sampleItemsPtr_.reset(new List<searchItem>(nPoints_));
            List<searchItem>& sampleItems = sampleItemsPtr_();
            forAll(sampleItems, itemI)
            {
                searchItem& curItem = sampleItems[itemI];
                curItem.procID()    = Pstream::myProcNo();
                curItem.cellLabel() = itemI;
                curItem.position()  =
                    min_+(max_-min_)/(nPoints_-1)*itemI + offset;
                curItem.seed()      = 0;
                curItem.zoneID()    = -1;
            }
        }
        else
        {
            sampleItemsPtr_.reset(new List<searchItem>(0));
        }

        if(log_) Info<<"    Inital interpolation cell search"<<endl;

        gradientSearch gs(refCast<const fvMesh>(obr_));

        gs.search(sampleItemsPtr_);
    }
}


void Foam::wakeVelocity::execute()
{
    // Do nothing - only valid on write
}


void Foam::wakeVelocity::end()
{
    // Do nothing - only valid on write
}


void Foam::wakeVelocity::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::wakeVelocity::write()
{
    if (active_)
    {
        const scalar& pi = Foam::constant::mathematical::pi;
        const scalar& t = obr_.time().value();
        const vector offset =
            amplitude_ * sin(2.0 * pi / period_ * t + pi / 180.0 * phaseLag_);

        {
            List<searchItem>& sampleItems = sampleItemsPtr_();
            forAll(sampleItems,itemI)
            {
                searchItem& curItem = sampleItems[itemI];
                const label& sampleI = curItem.cellLabel();
                curItem.position() =
                    min_+(max_-min_)/(nPoints_-1)*sampleI + offset;
            }
        }

        gradientSearch gs(refCast<const fvMesh>(obr_));

        gs.search(sampleItemsPtr_);

        vectorField sampledValues(nPoints_, vector::zero);

        const volVectorField& Ufield =
            obr_.lookupObject<volVectorField>(fieldName_);

        interpolationCellPoint<vector> interp(Ufield);

        {
            List<searchItem>& sampleItems = sampleItemsPtr_();
            forAll(sampleItems,itemI)
            {
                searchItem& curItem = sampleItems[itemI];
                sampledValues[curItem.cellLabel()] =
                    interp.interpolate(curItem.position(), curItem.seed());
            }
        }

        Pstream::listCombineGather<vector>(sampledValues,plusEqOp<vector>());

        if(Pstream::master())
        {
            if(log_) Info<<"    Writing wake velocities."<<endl;

            functionObjectFile::write();

            file()<<t;
            forAll(sampledValues,valueI)
            {
                file()<<"\t"<<sampledValues[valueI];
            }
            file()<<endl;

        }
    }
}


// ************************************************************************* //
