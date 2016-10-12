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

#include "bellerophonFvsPatchField.H"

#include "fvcDiv.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::bellerophonFvsPatchField<Type>::bellerophonFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    bellerophonPatch_(refCast<const bellerophonFvPatch>(p)),
    massFlux_(0.0),
    fieldName_("")
{
    if(notNull(this->dimensionedInternalField()))
    {
        fieldName_ = this->dimensionedInternalField().name();
    }
}


template<class Type>
Foam::bellerophonFvsPatchField<Type>::bellerophonFvsPatchField
(
    const bellerophonFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    bellerophonPatch_(refCast<const bellerophonFvPatch>(p)),
    massFlux_(ptf.massFlux_),
    fieldName_("")
{
    if (!isA<bellerophonFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "bellerophonFvsPatchField<Type>::bellerophonFvsPatchField\n"
            "("
                "const bellerophonFvsPatchField<Type>&, "
                "const fvPatch&, "
                "const DimensionedField<Type, surfaceMesh>&"
            ")"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }

    if(notNull(this->dimensionedInternalField()))
    {
        fieldName_ = this->dimensionedInternalField().name();
    }
}


template<class Type>
Foam::bellerophonFvsPatchField<Type>::bellerophonFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    bellerophonPatch_(refCast<const bellerophonFvPatch>(p)),
    massFlux_(dict.lookupOrDefault("massFlux",0.0)),
    fieldName_("")
{
    if (!isA<bellerophonFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "bellerophonFvsPatchField<Type>::bellerophonFvsPatchField"
            "("
                "const fvPatch&, "
                "const Field<Type>&, "
                "const dictionary&"
            ")",
            dict
        )   << "patch " << this->patch().index() << " not bellerophon type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    if(notNull(this->dimensionedInternalField()))
    {
        fieldName_ = this->dimensionedInternalField().name();
    }
}


template<class Type>
Foam::bellerophonFvsPatchField<Type>::bellerophonFvsPatchField
(
    const bellerophonFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    bellerophonPatch_(ptf.bellerophonPatch_),
    massFlux_(ptf.massFlux_),
    fieldName_("")
{
    if(notNull(this->dimensionedInternalField()))
    {
        fieldName_ = this->dimensionedInternalField().name();
    }
}

template<class Type>
Foam::bellerophonFvsPatchField<Type>::bellerophonFvsPatchField
(
    const bellerophonFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    bellerophonPatch_(ptf.bellerophonPatch_),
    massFlux_(ptf.massFlux_),
    fieldName_("")
{
    if(notNull(this->dimensionedInternalField()))
    {
        fieldName_ = this->dimensionedInternalField().name();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
namespace Foam
{

template<>
void bellerophonFvsPatchField<scalar>::enforceContinuity
(
    const surfaceScalarField& bField
)
{
    if(debug)
    {
        Info<<"Enforcing continuity..."<<endl;
    }

    // Methodology
    // 1. Correct fluxes to ensure fluxes across cell next to patch are
    //    divergence free
    //    Try: enforce continuity across faces between acceptor and regular
    //         cells - AG20150416

    // 2. Correct fluxes across patch faces to ensure global continuity

    // 1. Integrate over faces and sum up face are

    // Reference to the mesh
    const fvMesh& mesh = bField.mesh();

    // Faces areas
    const scalarField& internalMagSf = mesh.magSf().internalField();

    // Faces of the interface
    const labelList& interfaceFaces =
        bellerophon::Interpolation().interfaceFaces
        (
            bellerophonPatch_.bellerophonInterface::index()
        );

    // Flip map of the interface faces
    const boolList& interfaceFlipMap =
        bellerophon::Interpolation().interfaceFlipMap
        (
            bellerophonPatch_.bellerophonInterface::index()
        );

    // Reference to the field this patch belongs to
    // De-consting the field. DON'T TRY THIS AT HOME
    scalarField& iField =
        const_cast<scalarField&>(fvsPatchField<scalar>::internalField());

    // Sum of fluxes, positive into acceptor cells
    scalar interfaceSum = 0.0;
    scalar interfaceArea = 0.0;
    scalar interfaceAbsSum = 0.0;

    forAll(interfaceFaces,faceI)
    {
        const scalar& flux = iField[interfaceFaces[faceI]];
        if(interfaceFlipMap[faceI])
        {
            interfaceSum -= flux;
        }
        else
        {
            interfaceSum += flux;
        }
        interfaceAbsSum += fabs(flux);
        interfaceArea+=internalMagSf[faceI];
    }

    if(Pstream::parRun())
    {
        Pstream::gather(interfaceSum, sumOp<scalar>());
        interfaceSum-=massFlux_;
        Pstream::scatter(interfaceSum);
        reduce(interfaceAbsSum, sumOp<scalar>());
        reduce(interfaceArea, sumOp<scalar>());
    }
    else
    {
        interfaceSum-=massFlux_;
    }

    const scalar corr = interfaceSum / interfaceArea;

    Info<<"Interface error : "<<interfaceSum<<" / "<<interfaceAbsSum<<endl;

    interfaceSum = -massFlux_;

    forAll(interfaceFaces,faceI)
    {
        if(interfaceFlipMap[faceI])
        {
            iField[interfaceFaces[faceI]] += corr * internalMagSf[faceI];
            interfaceSum -= iField[interfaceFaces[faceI]];
        }
        else
        {
            iField[interfaceFaces[faceI]] -= corr * internalMagSf[faceI];
            interfaceSum += iField[interfaceFaces[faceI]];
        }
    }

    // List of cells adjoint to the patch faces
    const labelList& fcs = bellerophonPatch_.faceCells();

    // List of Cells of the mesh
    const cellList& cs = mesh.cells();

    // Number of internal faces
    const label nIntFaces = mesh.nInternalFaces();

    // Reference to owner labels of the faces
    const labelList& own = mesh.owner();

    // Sum of continuity error per adjoint cell
    scalar sum;

    forAll(fcs, localFaceI)
    {
        // Current internal cell
        const label curCell = fcs[localFaceI];

        // Sum of fluxes across faces of current internal cell
        sum = 0.0;

        // Labels of the faces of the current internal cell
        const labelList& curCellFaces = cs[curCell];
        forAll(curCellFaces,faceI)
        {
            // Label of the face of the current cell
            const label curFace = curCellFaces[faceI];

            if( curFace < nIntFaces )
            {
                if( own[curFace] == curCell )
                {
                    // On owner side
                    sum += iField[curFace];
                }
                else
                {
                    // On neighbour side
                    sum -= iField[curFace];
                }
            }
            else
            {
                // boundary faces are always owned by adjoint cell

                // Get ID patch of the boundary face
                const label patchID = mesh.boundaryMesh().whichPatch(curFace);

                const label localLabel = curFace - mesh.boundaryMesh()[patchID].start();
                const scalarList& bValues = bField.boundaryField()[patchID];
                if(bValues.size() > localLabel)
                {
                    sum += bValues[localLabel];
                }
            }
        }

        // Correct flux on the patch
        (*this)[localFaceI] -= sum;
    }

    if(bellerophonPatch_.holeInterface() != -1)
    {
        // Reference to the boundary Mesh (access to patches)
        const polyBoundaryMesh& bm = mesh.boundaryMesh();

        // Reference to owner labels of the faces
        const labelList& nei = mesh.neighbour();

        const label holeInterface = bellerophonPatch_.holeInterface();

        // Faces of the interface
        const labelList& holeFaces =
            bellerophon::Interpolation().interfaceFaces(holeInterface);

        // Flip map of the interface faces
        const boolList& holeFlipMap =
            bellerophon::Interpolation().interfaceFlipMap(holeInterface);

        forAll(holeFaces, faceI)
        {
            const label f = holeFaces[faceI];
            const label cellI = holeFlipMap[faceI] ? nei[f] : own[f];
            scalar fluxSum = 0.0;
            const labelList& faces = cs[cellI];
            forAll(faces, i)
            {
                const label curFace = faces[i];
                if( curFace < nIntFaces )
                {
                    if( own[curFace] == cellI )
                    {
                        // On owner side
                        fluxSum += iField[curFace];
                    }
                    else
                    {
                        // On neighbour side
                        fluxSum -= iField[curFace];
                    }
                }
                else
                {
                    // boundary faces are always owned by attached cell

                    // Get ID patch of the boundary face
                    const label patchID = bm.whichPatch(curFace);

                    const label localLabel = curFace - bm[patchID].start();
                    const scalarList& bValues = bField.boundaryField()[patchID];
                    if(bValues.size() > localLabel)
                    {
                        fluxSum += bValues[localLabel];
                    }
                }
            }
            iField[f] -= holeFlipMap[faceI] ? -fluxSum: fluxSum;
        }
    }

    // 2. Enforce global continuity
    const scalarField& areas = bellerophonPatch_.magSf();
    scalar area = 0.0;
    scalar massDefect = 0.0;

    forAll(areas, areaI)
    {
        area += areas[areaI];
        massDefect += (*this)[areaI];
    }

    Pstream::gather(massDefect, sumOp<scalar>());
    massDefect-=massFlux_;
    Pstream::scatter(massDefect);

    reduce(area, sumOp<scalar>());

    // Calculate specific continuity error
    massDefect/=area;

    // Correct values
    forAll(areas, areaI)
    {
        (*this)[areaI] -= massDefect*areas[areaI];
    }

}

template<>
void bellerophonFvsPatchField<vector>::enforceContinuity
(
    const surfaceVectorField& bField
)
{
    if(debug)
    {
        Info<<"Enforcing continuity..."<<endl;
    }

    // TODO is this necessary????

    // Reference to the mesh
    const fvMesh& mesh = bField.mesh();

    // Faces areas
    const scalarField& internalMagSf = mesh.magSf().internalField();

    // Faces vectors
    const vectorField& internalSf = mesh.Sf().internalField();

    // Faces of the interface
    const labelList& interfaceFaces =
        bellerophon::Interpolation().interfaceFaces
        (
            bellerophonPatch_.bellerophonInterface::index()
        );

    // Flip map of the interface faces
    const boolList& interfaceFlipMap =
        bellerophon::Interpolation().interfaceFlipMap
        (
            bellerophonPatch_.bellerophonInterface::index()
        );

    // Reference to the field this patch belongs to
    // De-consting the field. DON'T TRY THIS AT HOME
    vectorField& iField =
        const_cast<vectorField&>(fvsPatchField<vector>::internalField());

    // List of Cells of the mesh
    const cellList& cs = mesh.cells();

    // Number of internal faces
    const label nIntFaces = mesh.nInternalFaces();

    // Reference to owner labels of the faces
    const labelList& own = mesh.owner();

    // Reference to neighbour labels of the faces
    const labelList& nei= mesh.neighbour();

    // Sum of continuity error per adjoint cell
    scalar sum;

    forAll(interfaceFaces, interfaceFaceI)
    {
        // face of the interface
        const label f = interfaceFaces[interfaceFaceI];

        // Current internal cell
        const label curCell =
            interfaceFlipMap[interfaceFaceI] ? nei[f] : own[f];

        // Sum of fluxes across faces of current internal cell
        sum = 0.0;

        // Labels of the faces of the current internal cell
        const labelList& curCellFaces = cs[curCell];
        forAll(curCellFaces,faceI)
        {
            // Label of the face of the current cell
            const label curFace = curCellFaces[faceI];

            if( curFace < nIntFaces )
            {
                if( own[curFace] == curCell )
                {
                    // On owner side
                    sum += iField[curFace] & internalSf[curFace];
                }
                else
                {
                    // On neighbour side
                    sum -= iField[curFace] & internalSf[curFace];
                }
            }
            else
            {
                // boundary faces are always owned by adjoint cell

                // Get ID patch of the boundary face
                const label patchID = mesh.boundaryMesh().whichPatch(curFace);

                const label localLabel =
                    curFace - mesh.boundaryMesh()[patchID].start();
                const vectorField& bValues = bField.boundaryField()[patchID];
                if(bValues.size() > localLabel)
                {
                    sum += bValues[localLabel]
                             & mesh.Sf().boundaryField()[patchID][localLabel];
                }
            }
        }

        if(f < nIntFaces)
        {
            const vector corr = internalSf[f]*sum/sqr(internalMagSf[f]);
            if(interfaceFlipMap[interfaceFaceI])
            {
                // Cell in focus, i.e. cell behind interface, is neighbour cell.
                // Normal points inwards, positive sum of fluxes, so more flux
                // has to go into the cell
                iField[f] += corr;
            }
            else
            {
                // Cell in focus, i.e. cell behind interface, is owner cell
                // The other way round
                iField[f] -= corr;
            }
        }
        else
        {
            // Get ID patch of the boundary face
            const label patchID = mesh.boundaryMesh().whichPatch(f);

            const label localLabel =
                f - mesh.boundaryMesh()[patchID].start();

            const vector corr =
                         sum
                         *
                         mesh.Sf().boundaryField()[patchID][localLabel]
                         /
                         sqr(mesh.magSf().boundaryField()[patchID][localLabel]);

            vectorField& bValues =
                const_cast<vectorField&>(refCast<const vectorField>(bField.boundaryField()[patchID]));

            if(bValues.size() > localLabel)
            {
                bValues[localLabel] -= corr;
            }
        }
    }

    if(bellerophonPatch_.holeInterface() != -1)
    {
        // Reference to the boundary Mesh (access to patches)
        const polyBoundaryMesh& bm = mesh.boundaryMesh();

        const label holeInterface = bellerophonPatch_.holeInterface();

        // Faces of the interface
        const labelList& holeFaces =
            bellerophon::Interpolation().interfaceFaces(holeInterface);

        // Flip map of the interface faces
        const boolList& holeFlipMap =
            bellerophon::Interpolation().interfaceFlipMap(holeInterface);

        forAll(holeFaces, interfaceFaceI)
        {
            // face of the interface
            const label f = holeFaces[interfaceFaceI];

            // Current internal cell
            const label curCell =
                holeFlipMap[interfaceFaceI] ? nei[f] : own[f];

            // Sum of fluxes across faces of current internal cell
            sum = 0.0;

            // Labels of the faces of the current internal cell
            const labelList& curCellFaces = cs[curCell];
            forAll(curCellFaces,faceI)
            {
                // Label of the face of the current cell
                const label curFace = curCellFaces[faceI];

                if( curFace < nIntFaces )
                {
                    if( own[curFace] == curCell )
                    {
                        // On owner side
                        sum += iField[curFace] & internalSf[curFace];
                    }
                    else
                    {
                        // On neighbour side
                        sum -= iField[curFace] & internalSf[curFace];
                    }
                }
                else
                {
                    // boundary faces are always owned by adjoint cell

                    // Get ID patch of the boundary face
                    const label patchID = mesh.boundaryMesh().whichPatch(curFace);

                    const label localLabel =
                        curFace - mesh.boundaryMesh()[patchID].start();
                    const vectorField& bValues = bField.boundaryField()[patchID];
                    if(bValues.size() > localLabel)
                    {
                        sum += bValues[localLabel]
                                 & mesh.Sf().boundaryField()[patchID][localLabel];
                    }
                }
            }

            if(f < nIntFaces)
            {
                const vector corr = internalSf[f]*sum/sqr(internalMagSf[f]);
                if(holeFlipMap[interfaceFaceI])
                {
                    // Cell in focus, i.e. cell behind interface, is neighbour cell.
                    // Normal points inwards, positive sum of fluxes, so more flux
                    // has to go into the cell
                    iField[f] += corr;
                }
                else
                {
                    // Cell in focus, i.e. cell behind interface, is owner cell
                    // The other way round
                    iField[f] -= corr;
                }
            }
            else
            {
                // Get ID patch of the boundary face
                const label patchID = bm.whichPatch(f);

                const label localLabel =
                    f - bm[patchID].start();

                const vector corr =
                             sum
                             *
                             mesh.Sf().boundaryField()[patchID][localLabel]
                             /
                             sqr
                             (
                                 mesh.magSf()
                                     .boundaryField()[patchID][localLabel]
                             );

                vectorField& bValues =
                    const_cast<vectorField&>
                    (
                        refCast<const vectorField>
                            (bField.boundaryField()[patchID])
                    );

                if(bValues.size() > localLabel)
                {
                    bValues[localLabel] -= corr;
                }
            }
        }
    }
}

template<class Type>
void bellerophonFvsPatchField<Type>::enforceContinuity
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& bField
)
{}

template<>
void bellerophonFvsPatchField<scalar>::enforceContinuity()
{
    const objectRegistry& reg =
        bellerophonPatch_.boundaryMesh().mesh().thisDb();

    if(reg.foundObject<surfaceScalarField>(fieldName_))
    {
        const surfaceScalarField& bField = reg.lookupObject<surfaceScalarField>
            (fieldName_);

        enforceContinuity(bField);
    }
}

template<>
void bellerophonFvsPatchField<vector>::enforceContinuity()
{
    const objectRegistry& reg =
        bellerophonPatch_.boundaryMesh().mesh().thisDb();

    if(reg.foundObject<surfaceVectorField>(fieldName_))
    {
        const surfaceVectorField& bField = reg.lookupObject<surfaceVectorField>
            (fieldName_);

        enforceContinuity(bField);
    }
}

template<class Type>
void bellerophonFvsPatchField<Type>::enforceContinuity()
{}

template<class Type>
void bellerophonFvsPatchField<Type>::assignHoleInterface(const Type& rhs)
{
    if(bellerophonPatch_.holeInterface() != -1)
    {
        if(debug)
        {
            Info<<"Assigning "<<rhs<<" to hole in field "
                <<fieldName_<<endl;
        }

        const label holeInterfaceID = bellerophonPatch_.holeInterface();

        const labelList& interfaceFaces =
            bellerophon::Interpolation().interfaceFaces(holeInterfaceID);

        const boolList& interfaceFlipMap =
            bellerophon::Interpolation().interfaceFlipMap(holeInterfaceID);

        Field<Type>& iField =
            const_cast<Field<Type>&>(fvsPatchField<Type>::internalField());

        forAll(interfaceFaces, faceI)
        {
            if(interfaceFlipMap[faceI])
            {
                iField[interfaceFaces[faceI]] = -rhs;
            }
            else
            {
                iField[interfaceFaces[faceI]] = rhs;
            }
        }
    }
    else if (debug)
    {
        Info<<"Assigning "<<rhs<<" to hole in field "<<fieldName_<<endl;
    }
}

template<class Type>
void bellerophonFvsPatchField<Type>::correctSum()
{}

template<>
void bellerophonFvsPatchField<scalar>::correctSum()
{
    if
    (
        fieldName_ != ""
        &&
        bellerophon::Interpolation().enforceContinuity
        (
            fieldName_
        )
    )
    {
        enforceContinuity();
    }
}

template<>
void bellerophonFvsPatchField<vector>::correctSum()
{
    if
    (
        fieldName_ != ""
        &&
        bellerophon::Interpolation().enforceContinuity
        (
            fieldName_
        )
    )
    {
        enforceContinuity();
    }
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator=
(
    const UList<Type>& rhs
)
{
    fvsPatchField<Type>::operator=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator=
(
    const fvsPatchField<Type>& rhs
)
{
    fvsPatchField<Type>::operator=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator+=
(
    const fvsPatchField<Type>& rhs
)
{
    fvsPatchField<Type>::operator+=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator-=
(
    const fvsPatchField<Type>& rhs
)
{
    fvsPatchField<Type>::operator-=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator*=
(
    const fvsPatchField<scalar>& rhs
)
{
    fvsPatchField<Type>::operator*=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator/=
(
    const fvsPatchField<scalar>& rhs
)
{
    fvsPatchField<Type>::operator/=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator+=
(
    const Field<Type>& rhs
)
{
    fvsPatchField<Type>::operator+=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator-=
(
    const Field<Type>& rhs
)
{
    fvsPatchField<Type>::operator-=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator*=
(
    const Field<scalar>& rhs
)
{
    fvsPatchField<Type>::operator*=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator/=
(
    const Field<scalar>& rhs
)
{
    fvsPatchField<Type>::operator/=(rhs);

    correctSum();
}

template<>
void bellerophonFvsPatchField<scalar>::operator=
(
    const scalar& rhs
)
{
    fvsPatchField<scalar>::operator=(rhs);

    // This is a hack for the transient flux correction
    // ddtScheme::fvcDdtPhiCoeff shall be zero on this patch and on the hole
    // as well faces

    assignHoleInterface(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator=
(
    const Type& rhs
)
{
    fvsPatchField<Type>::operator=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator+=
(
    const Type& rhs
)
{
    fvsPatchField<Type>::operator+=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator-=
(
    const Type& rhs
)
{
    fvsPatchField<Type>::operator-=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator*=
(
    const scalar rhs
)
{
    fvsPatchField<Type>::operator*=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::operator/=
(
    const scalar rhs
)
{
    fvsPatchField<Type>::operator/=(rhs);

    correctSum();
}

template<class Type>
void bellerophonFvsPatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("massFlux")<<massFlux_<< token::END_STATEMENT << nl;

    coupledFvsPatchField<Type>::write(os);
}

}

// ************************************************************************* //
