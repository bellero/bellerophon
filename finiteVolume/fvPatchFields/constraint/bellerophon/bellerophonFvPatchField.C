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

#include "bellerophonFvPatchField.H"
#include "transformField.H"
#include "fvMatrix.H"
#include "fvcGrad.H"
#include "bellerophonLduMatrix.H"

#include "cellSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
bellerophonFvPatchField<Type>::bellerophonFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    bellerophonInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    interface_(refCast<const bellerophonInterface>(p)),
    forceInterpolation_(false),
    bound_(false)/*,
    fixesValue_(false)*/
{}


template<class Type>
bellerophonFvPatchField<Type>::bellerophonFvPatchField
(
    const bellerophonFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    bellerophonInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    interface_(refCast<const bellerophonInterface>(p)),
    forceInterpolation_(ptf.forceInterpolation_),
    bound_(ptf.bound_)/*,
    fixesValue_(ptf.fixesValue_)*/
{
    if (!isA<bellerophonFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "bellerophonFvPatchField<Type>::bellerophonFvPatchField"
            "("
                "const bellerophonFvPatchField<Type>& ,"
                "const fvPatch&, "
                "const DimensionedField<Type, volMesh>&"
                "const fvPatchFieldMapper& mapper"
                ")"
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
bellerophonFvPatchField<Type>::bellerophonFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    bellerophonInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict),
    interface_(refCast<const bellerophonInterface>(p)),
    forceInterpolation_(dict.lookupOrDefault("forceInterpolation",true)),
    bound_(dict.lookupOrDefault("bound",false))/*,
    fixesValue_(dict.lookupOrDefault("fixesValue",false))*/
{
    if (!isA<bellerophonFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "bellerophonFvPatchField<Type>::bellerophonFvPatchField"
            "("
                "const fvPatch&, "
                "const DimensionedField<Type, volMesh>&, "
                "const dictionary&"
            ")",
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    if(debug && forceInterpolation_)
    {
        Info<<"Forcing interpolation for field "
            <<this->dimensionedInternalField().name()<<endl;
    }

    if(debug && bound_)
    {
        Info<<"Bounding field "
            <<this->dimensionedInternalField().name()
            <<" during interpolation."<<endl;
    }

//     if(debug && fixesValue_)
//     {
//         Info<<"Will pretend to fix value on  "
//             <<this->dimensionedInternalField().name()
//             <<" if asked."<<endl;
//     }

    // Don't evaluate here, as other patch fields might not be constructed yet
    // This causes problems when calculating the gradient on evaluation
//     if (!dict.found("value") && this->coupled())
//     {
//         this->evaluate(Pstream::blocking);
//     }
}


template<class Type>
bellerophonFvPatchField<Type>::bellerophonFvPatchField
(
    const bellerophonFvPatchField<Type>& ptf
)
:
    bellerophonInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    interface_(ptf.interface_),
    forceInterpolation_(ptf.forceInterpolation_),
    bound_(ptf.bound_)/*,
    fixesValue_(ptf.fixesValue_)*/
{}


template<class Type>
bellerophonFvPatchField<Type>::bellerophonFvPatchField
(
    const bellerophonFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    bellerophonInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    interface_(ptf.interface_),
    forceInterpolation_(ptf.forceInterpolation_),
    bound_(ptf.bound_)/*,
    fixesValue_(ptf.fixesValue_)*/
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::bellerophonFvPatchField<Type>::coupled() const
{
    if (Pstream::parRun() || this->interface_.nFaces())
    {
        return true;
    }
    else
    {
        return false;
    }
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::bellerophonFvPatchField<Type>::patchNeighbourField() const
{
    return this->patchInternalField();
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::bellerophonFvPatchField<Type>::snGrad() const
{
    return Field<Type>(this->interface_.nFaces(),pTraits<Type>::zero);
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::bellerophonFvPatchField<Type>::snGrad(const scalarField&) const
{
    return this->snGrad();
}

template<>
Foam::tmp<Foam::Field<scalar> >
Foam::bellerophonFvPatchField<scalar>::snGrad() const
{
    // TODO gradient might be calculated several times

    const word iFieldName = fvPatchField<scalar>::dimensionedInternalField().name();

    const GeometricField<scalar, fvPatchField, volMesh>& iField =
        this->db().lookupObject<GeometricField<scalar, fvPatchField, volMesh> >
        (
            iFieldName
        );

    tmp<volVectorField> tIGrad = fvc::grad(iField);

    volVectorField& iGrad = tIGrad();

    forAll(iGrad.boundaryField(), patchI)
    {
        if(isA<bellerophonFvPatchField<vector> >(iGrad.boundaryField()[patchI]))
        {
            iGrad.boundaryField()[patchI].initEvaluate(Pstream::nonBlocking);
        }
    }

    const label myIndex = this->patch().index();

    return Field<scalar>(this->patch().Sf() & iGrad.boundaryField()[myIndex]);
}

template<>
Foam::tmp<Foam::Field<vector> >
Foam::bellerophonFvPatchField<vector>::snGrad() const
{
    // TODO gradient might be calculated several times

    const word iFieldName = fvPatchField<vector>::dimensionedInternalField().name();

    const GeometricField<vector, fvPatchField, volMesh>& iField =
        this->db().lookupObject<GeometricField<vector, fvPatchField, volMesh> >
        (
            iFieldName
        );

    tmp<volTensorField> tIGrad = fvc::grad(iField);

    volTensorField& iGrad = tIGrad();

    forAll(iGrad.boundaryField(), patchI)
    {
        if(isA<bellerophonFvPatchField<tensor> >(iGrad.boundaryField()[patchI]))
        {
            iGrad.boundaryField()[patchI].initEvaluate(Pstream::nonBlocking);
        }
    }

    const label myIndex = this->patch().index();

    return Field<vector>(this->patch().Sf() & iGrad.boundaryField()[myIndex]);
}

template<class Type>
void bellerophonFvPatchField<Type>::initEvaluate(const Pstream::commsTypes commsType)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Interpolate internal field
    if(interface_.index() == 0)
    {
        bellerophon::Interpolation().interpolate<Type>
        (
            iFieldName(),
            commsType,
            bound_
        );
    }

    if
    (
        bellerophon::Interpolation().forceZero(iFieldName())
    )
    {
        if(debug || true)
        {
            Info<<"Forcing zero to cells in field "<<iFieldName()<<endl;
        }

        // Faces between live cells and acceptor cells
        const labelList& interfaceFaces =
             bellerophon::Interpolation().interfaceFaces(interface_.index());

        // Flip map for interface faces
        const boolList& interfaceFlipMap =
             bellerophon::Interpolation().interfaceFlipMap(interface_.index());

        // Access to the values in the field
        Field<Type>& iField = const_cast<Field<Type>&>(this->internalField());

        // Zero
        const Type& zero = pTraits<Type>::zero;

        // The mesh
        const fvMesh& mesh = this->patch().boundaryMesh().mesh();

        const labelList& own = mesh.owner();
        const labelList& nei = mesh.neighbour();

        forAll(interfaceFaces, faceI)
        {
            const label f = interfaceFaces[faceI];
            const label c = interfaceFlipMap[faceI] ? nei[f] : own[f];
            iField[c] = zero;
        }

        if(interface().holeInterface() >= 0)
        {
            if(debug || true)
            {
                Info<<"Forcing zero to cells around hole in field "
                    <<iFieldName()<<endl;
            }

            // Index of the hole interface
            const label holeIndex = interface().holeInterface();

            // Faces between live cells and acceptor cells of the interface
            const labelList& holeFaces =
                bellerophon::Interpolation().interfaceFaces(holeIndex);

            // Flip map of the interface faces of the hole
            const boolList& holeFlipMap =
                bellerophon::Interpolation().interfaceFlipMap(holeIndex);

            forAll(holeFaces, faceI)
            {
                const label f = holeFaces[faceI];
                const label c = holeFlipMap[faceI] ? nei[f] : own[f];
                iField[c] = zero;
            }

//             // Faces between live cells and acceptor cells of the interface
//             const labelList& newLiveCells =
//                 bellerophon::Interpolation().newLiveCells();
//
// //             forAll(newLiveCells, c) iField[newLiveCells[c]] = zero;
//
//             cellSet newLiveSet(mesh,"newLiveCells", newLiveCells.size());
//             forAll(newLiveCells, c)
//             {
//                 iField[newLiveCells[c]] = zero;
//                 newLiveSet.insert(newLiveCells[c]);
//             }
//             newLiveSet.instance() = mesh.time().timeName();
//             newLiveSet.write();
        }
    }
    else if(debug)
    {
        Info<<"Not forcing zero for field "<<iFieldName()<<endl;
    }

    fvPatchField<Type>::operator==(this->patchInternalField());

    fvPatchField<Type>::initEvaluate(commsType);
}

template<class Type>
const Foam::word& bellerophonFvPatchField<Type>::iFieldName() const
{
    return fvPatchField<Type>::dimensionedInternalField().name();
}

template<class Type>
void bellerophonFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    // Ordinary coupled patches do stuff here, as there is a coupling in both
    // directions e.g. processor, ami, cyclic...
    // Overlapping grids do not have a direct coupling from the patch to the
    // overlapped grid, only the other way round. So therefore this mechanism
    // will not work. Coupling of the matrix system is done by the solver and
    // the matrix, i.e. bellerophonLduMatrix and bellerophonPBiCG
}

template<class Type>
void bellerophonFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    // Ordinary coupled patches do stuff here, as there is a coupling in both
    // directions e.g. processor, ami, cyclic...
    // Overlapping grids do not have a direct coupling from the patch to the
    // overlapped grid, only the other way round. So therefore this mechanism
    // will not work. Coupling of the matrix system is done by the solver and
    // the matrix, i.e. bellerophonLduMatrix and bellerophonPBiCG
}

template<class Type>
void bellerophonFvPatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("forceInterpolation")<<forceInterpolation_
            << token::END_STATEMENT << nl;
    os.writeKeyword("bound")<<bound_
            << token::END_STATEMENT << nl;
    coupledFvPatchField<Type>::write(os);
}

template<class Type>
tmp<Field<Type> > bellerophonFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one) * w;
}


template<class Type>
tmp<Field<Type> > bellerophonFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*(1.0 - w);
}

template<class Type>
tmp<Field<Type> > bellerophonFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return -Type(pTraits<Type>::one)*deltaCoeffs;
}

template<class Type>
tmp<Field<Type> > bellerophonFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}

template<class Type>
tmp<Field<Type> > bellerophonFvPatchField<Type>::gradientBoundaryCoeffs(
    const scalarField& deltaCoeffs
) const
{
    return Type(pTraits<Type>::one)*deltaCoeffs;
}

template<class Type>
tmp<Field<Type> > bellerophonFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}

template<class Type>
void bellerophonFvPatchField<Type>::operator=(const UList<Type>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator=(v);
    }
}

template<class Type>
void bellerophonFvPatchField<Type>::operator=(const fvPatchField<Type>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator+=(const fvPatchField<Type>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator+=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator-=(const fvPatchField<Type>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator-=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator*=(const fvPatchField<scalar>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator*=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator/=(const fvPatchField<scalar>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator/=(v);
    }
}

template<class Type>
void bellerophonFvPatchField<Type>::operator+=(const Field<Type>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator+=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator-=(const Field<Type>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator-=(v);
    }
}

template<class Type>
void bellerophonFvPatchField<Type>::operator*=(const Field<scalar>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator*=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator/=(const Field<scalar>& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator/=(v);
    }
}

template<class Type>
void bellerophonFvPatchField<Type>::operator=(const Type& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator+=(const Type& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator+=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator-=(const Type& v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator-=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator*=(const scalar v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator*=(v);
    }
}
template<class Type>
void bellerophonFvPatchField<Type>::operator/=(const scalar v)
{
    if(forceInterpolation_)
    {
        this->initEvaluate();
    }
    else
    {
        fvPatchField<Type>::operator/=(v);
    }
}


} // End namespace Foam

// ************************************************************************* //
