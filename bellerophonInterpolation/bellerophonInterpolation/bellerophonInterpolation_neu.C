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

#include "bellerophonInterpolation.H"
#include "meshTools.H"
#include "mergePoints.H"
#include "mapDistribute.H"
#include "labelList.H"
#include "remainingLabelList.H"
#include "coupledPolyPatch.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    //- Combine operator for interpolateToSource/Target
    template<class Type, class BinaryOp>
    class MyCombineBinaryOp
    {
        const BinaryOp& bop_;

        public:

            MyCombineBinaryOp(const BinaryOp& bop)
            :
                bop_(bop)
            {}

            void operator()
            (
                Type& x,
                const label faceI,
                const Type& y,
                const scalar weight
            ) const
            {
                x = bop_(x, weight*y);
            }
    };

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Patch, class Mesh>
bellerophonInterpolation<Patch, Mesh>::bellerophonInterpolation
(
    const Patch& patch,
    const Mesh& mesh
)
:
mesh_(mesh),
patch_(patch)
//     reverseTarget_(reverseTarget),
//     singlePatchProc_(-999),
//     srcAddress_(),
//     srcWeights_(),
//     tgtAddress_(),
//     tgtWeights_(),
//     treePtr_(NULL),
//     startSeedI_(0),
//     triMode_(triMode),
//     srcMapPtr_(NULL),
//     tgtMapPtr_(NULL)
{
     update(patch);
}

/*
template<class Patch, class Mesh>
Foam::bellerophonInterpolation<Patch, Mesh>::bellerophonInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget
)
:
    reverseTarget_(reverseTarget),
    singlePatchProc_(-999),
    srcAddress_(),
    srcWeights_(),
    tgtAddress_(),
    tgtWeights_(),
    treePtr_(NULL),
    startSeedI_(0),
    triMode_(triMode),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    Info<< "bellerophonInterpolation::bellerophonInterpolation(...)"
        << " at bellerophonInterpolation.C:1529" << nl << endl;

    label srcSize = returnReduce(srcPatch.size(), sumOp<label>());
    label tgtSize = returnReduce(tgtPatch.size(), sumOp<label>());

    Info<< "bellerophon: Creating addressing and weights between "
        << srcSize << " source faces and " << tgtSize << " target faces"
        << endl;

    if (surfPtr.valid())
    {
        // create new patches for source and target
        pointField srcPoints = srcPatch.points();
        primitivePatch srcPatch0
        (
            SubList<face>
            (
                srcPatch,
                srcPatch.size(),
                0
            ),
            srcPoints
        );

        if (debug)
        {
            OFstream os("amiSrcPoints.obj");
            forAll(srcPoints, i)
            {
                meshTools::writeOBJ(os, srcPoints[i]);
            }
        }

        pointField tgtPoints = tgtPatch.points();
        primitivePatch tgtPatch0
        (
            SubList<face>
            (
                tgtPatch,
                tgtPatch.size(),
                0
            ),
            tgtPoints
        );

        if (debug)
        {
            OFstream os("amiTgtPoints.obj");
            forAll(tgtPoints, i)
            {
                meshTools::writeOBJ(os, tgtPoints[i]);
            }
        }

        Info<<"bellerophonInterpolation::bellerophonInterpolation(...) [1]"<<nl<<endl;

        // map source and target patches onto projection surface
        projectPointsToSurface(surfPtr(), srcPoints);

        Info<<"bellerophonInterpolation::bellerophonInterpolation(...) [2]"<<nl<<endl;

        projectPointsToSurface(surfPtr(), tgtPoints);

        Info<<"bellerophonInterpolation::bellerophonInterpolation(...) [3]"<<nl<<endl;

        // calculate bellerophon interpolation
        update(srcPatch0, tgtPatch0);

        Info<<"bellerophonInterpolation::bellerophonInterpolation(...) [4]"<<nl<<endl;
    }
    else
    {
        update(srcPatch, tgtPatch);
    }
}
*/


template<class Patch, class Mesh>
bellerophonInterpolation<Patch, Mesh>::bellerophonInterpolation
(
    const bellerophonInterpolation<Patch, Mesh>& fineBellerophon,
    const labelList& restrictAddressing,
    const labelList& someList
)
//     reverseTarget_(finebellerophon.reverseTarget_),
//     singlePatchProc_(finebellerophon.singlePatchProc_),
//     srcAddress_(),
//     srcWeights_(),
//     tgtAddress_(),
//     tgtWeights_(),
//     treePtr_(NULL),
//     startSeedI_(0),
//     triMode_(finebellerophon.triMode_),
//     srcMapPtr_(NULL),
//     tgtMapPtr_(NULL)
{

/*    label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    label neighbourCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    if (debug & 2)
    {
        Pout<< "bellerophon: Creating addressing and weights as agglomeration of bellerophon :"
            << " source:" << finebellerophon.srcAddress().size()
            << " target:" << finebellerophon.tgtAddress().size()
            << " coarse source size:" << sourceCoarseSize
            << " neighbour source size:" << neighbourCoarseSize
            << endl;
    }

    if
    (
        finebellerophon.srcAddress().size() != sourceRestrictAddressing.size()
     || finebellerophon.tgtAddress().size() != targetRestrictAddressing.size()
    )
    {
        FatalErrorIn
        (
            "bellerophonInterpolation<Patch, Mesh>::bellerophonInterpolation\n"
            "(\n"
            "    const bellerophonInterpolation<Patch, Mesh>&,\n"
            "    const label,\n"
            "    const labelList&\n"
            ")"
        )   << "Size mismatch." << nl
            << "Source patch size:" << finebellerophon.srcAddress().size() << nl
            << "Source agglomeration size:"
            << sourceRestrictAddressing.size() << nl
            << "Target patch size:" << finebellerophon.tgtAddress().size() << nl
            << "Target agglomeration size:"
            << targetRestrictAddressing.size()
            << exit(FatalError);
    }


    // Agglomerate addresses and weights

    agglomerate
    (
        finebellerophon.tgtMapPtr_,
        finebellerophon.srcMagSf(),
        finebellerophon.srcAddress(),
        finebellerophon.srcWeights(),

        sourceRestrictAddressing,
        targetRestrictAddressing,

        srcMagSf_,
        srcAddress_,
        srcWeights_,
        tgtMapPtr_
    );

    //if (tgtMapPtr_.valid())
    //{
    //    Pout<< "tgtMap:" << endl;
    //    string oldPrefix = Pout.prefix();
    //    Pout.prefix() = oldPrefix + "  ";
    //    tgtMapPtr_().printLayout(Pout);
    //    Pout.prefix() = oldPrefix;
    //}

    agglomerate
    (
        finebellerophon.srcMapPtr_,
        finebellerophon.tgtMagSf(),
        finebellerophon.tgtAddress(),
        finebellerophon.tgtWeights(),

        targetRestrictAddressing,
        sourceRestrictAddressing,

        tgtMagSf_,
        tgtAddress_,
        tgtWeights_,
        srcMapPtr_
    );

    //if (srcMapPtr_.valid())
    //{
    //    Pout<< "srcMap:" << endl;
    //    string oldPrefix = Pout.prefix();
    //    Pout.prefix() = oldPrefix + "  ";
    //    srcMapPtr_().printLayout(Pout);
    //    Pout.prefix() = oldPrefix;
    //}
    */
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class Patch, class Mesh>
bellerophonInterpolation<Patch, Mesh>::~bellerophonInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Patch, class Mesh>
void bellerophonInterpolation<Patch, Mesh>::update
(
    const polyPatch& patch
)
{
    Info<<"bellerophonInterpolation: Updating Interpolation Information"
        <<nl<<endl;

#define donorFaceLinks

#ifdef donorFaceLinks
    static label file = 0;
    label count = 0;
    OFstream os("donorFaceLinks_" + name(file++) + ".obj");
#endif

    // TODO:
    // - check wheater it is really necessary to regenerate the interpolation
    //   ab ovo
    // - improve storage allocation
    // - find a neat way to search donor cells

    const label nCells = patch.size();
    donorCells_.setSize(nCells);
    donorWeights_.setSize(nCells);
    nDonors_ = 0;

    // first simple attempt:
    // only find primary donor cell with brute force search

    const cellList& cells_=mesh_.cells();
    const faceList& faces_=mesh_.faces();
    const labelListList& cellCells = mesh_.cellCells();
    const vectorField& cellCentres = mesh_.cellCenters();

    //cell centres instead of face centres
    // const vectorField& centres_=patch.faceCentres();
    const vectorField centres_=patch.faceCellCentres()();

    forAll(patch,faceI)
    {
        label cellI=0;
        bool found = false;
        point cellCentre=centres_[faceI];


        while(cellI<cells_.size() && !found)
        {
            // if the point is inside the cell
            if (mesh_.pointInCell(cellCentre,cellI))
            {
                // check if one of the cells faces is the face we search a donor
                // cell for
                // simple check: no face of the found cell should be the face
                // of the bellerophon patch

                bool ownFace=false;
                forAll(cells_[cellI],cellFaceI)
                {
                    if (faces_[cells_[cellI][cellFaceI]]==patch[faceI])
                    {
                        ownFace=true;
                    }
                }
                if(!ownFace)
                {
                    found=true;
//                     label nLocalDonors = 1 + cellCells[cellI].size();
//                     const labelList& localDonors = cellCells[cellI];
// 
//                     donorCells_[faceI].setSize( nLocalDonors );
//                     donorWeights_[faceI].setSize( nLocalDonors );
//                     labelList& localDonorCells = donorCells_[faceI];
//                     labelList& localDonorWeights = donorWeights_[faceI];
//                     localDonorCells[0] = cellI;
// 
//                     scalar denom = SMALL + magSqr(cellCentre-cellCentres[cellI]);
//                     scalar norm = denom;
//                     localDonorWeights[0] =
//                         1.0 / denom;
// 
//                     forAll(localDonors, donorI)
//                     {
//                         localDonorCells[donorI+1] = localDonors[donorI];
//                         denom = SMALL + magSqr(cellCentre - cellCentres[localDonors[donorI]]);
//                         localDonorWeights[donorI+1] = 1.0 / denom;
//                         norm += denom;
//                     }
// 
//                     norm = 1.0/norm;
// 
//                     forAll(localDonorWeights, donorI)
//                     {
//                         localDonorWeights[donorI] *= norm;
//                     }
// 
//                     nDonors_ += nLocalDonors;
                }
            }

            // if we don't have a donor, continue with the next cell
            if (!found)
            {
                cellI++;
            }
        }

        // gone through all cells for this face
        // if still not found, it is an orphan face
        if(!found)
        {
            FatalErrorIn
            (
                "bellerophonInterpolation<Patch, Mesh>::update"
                "("
                "    const primitivePatch& patch"
                ")"
            )
            << "No donor cell found for face " << faceI
            << " on patch with size " << patch.size() << nl << nl
            << "Please check mesh."
            << exit(FatalError);
        }

#ifdef donorFaceLinks
        os  << "v " << cellCentre.x() << " " << cellCentre.y() << " " << cellCentre.z() << nl
        << "v " << mesh_.cellCentres()[cellI].x() << " "
        << mesh_.cellCentres()[cellI].y() << " "
        << mesh_.cellCentres()[cellI].z() << nl
        << "l " << ++count-1 <<  " " << ++count << endl;
#endif

    }

    Info<<"Finished updating."<<nl<<endl;
}


template<class Patch, class Mesh>
label bellerophonInterpolation<Patch, Mesh>::calcDistribution
(
   const primitivePatch& patch
)
{
    label procI = 0;

    if (Pstream::parRun())
    {
        List<label> facesPresentOnProc(Pstream::nProcs(), 0);
        if (patch.size() > 0)
        {
            facesPresentOnProc[Pstream::myProcNo()] = 1;
        }
        else
        {
            facesPresentOnProc[Pstream::myProcNo()] = 0;
        }

        Pstream::gatherList(facesPresentOnProc);
        Pstream::scatterList(facesPresentOnProc);

        label nHaveFaces = sum(facesPresentOnProc);

        if (nHaveFaces > 1)
        {
            procI = -1;
            if (debug)
            {
                Info<< "bellerophonInterpolation::calcDistribution: "
                << "bellerophon split across multiple processors" << endl;
            }
        }
        else if (nHaveFaces == 1)
        {
            procI = findIndex(facesPresentOnProc, 1);
            if (debug)
            {
                Info<< "bellerophonInterpolation::calcDistribution: "
                << "bellerophon local to processor" << procI << endl;
            }
        }
    }

    // Either not parallel or no faces on any processor
    return procI;
}


// ************************************************************************* //
