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

Description
    Multiply a given vector (second argument) by the matrix or its transpose
    and return the result in the first argument.

\*---------------------------------------------------------------------------*/

#include "bellerophonLduMatrix.H"
#include "bellerophonInterpolation.H"
#include "interpolationItem.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::bellerophonLduMatrix::Amul
(
    scalarField& Apsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    scalar __restrict__ * ApsiPtr = Apsi.begin();

    const scalarField& psi = tpsi();

    const scalar* const __restrict__ diagPtr = diag().begin();

    const scalar* const __restrict__ psiPtr = psi.begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = modUpper().begin();
    const scalar* const __restrict__ lowerPtr = modLower().begin();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    register const label nCells = diag().size();
    for (register label cell=0; cell<nCells; cell++)
    {
        ApsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }

    register const label nFaces = modUpper().size();

    for (register label face=0; face<nFaces; face++)
    {
        ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
        ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
    }

    // Cache values of acceptor cells, because they might be changed by other
    // interfaces
    autoPtr<scalarField> ApsiCache;
    const labelList& acceptors = bellerophon::Interpolation().acceptorCells();
    if(acceptors.size() > 0)
    {
        const label nAcceptors = acceptors.size();
        ApsiCache.set(new scalarField(nAcceptors));
        scalar* __restrict__ cachePtr = ApsiCache().begin();
        const label* const __restrict__ acceptorRowPtr = acceptors.begin();
        for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
        {
            cachePtr[acceptorI] = ApsiPtr[acceptorRowPtr[acceptorI]];
        }
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    // Restore cached values
    if(acceptors.size() > 0)
    {
        const label nAcceptors = acceptors.size();
        scalar* __restrict__ cachePtr = ApsiCache().begin();
        const label* const __restrict__ acceptorRowPtr = acceptors.begin();
        for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
        {
            ApsiPtr[acceptorRowPtr[acceptorI]] = cachePtr[acceptorI];
        }
        ApsiCache.clear();
    }

    // Interpolate Psi to donor items
    scalar __restrict__ *interpolatedPsiPtr =
        interpolatedPsi_.begin();

    // Collect primary donor influence
    {
        const label* const __restrict__ pDC =
            primaryDonorCells_.begin();

        const scalar* const __restrict__ pDW =
            primaryDonorWeights_.begin();

        const label nDonors = primaryDonorCells_.size();

        for(label donorI = 0; donorI<nDonors; donorI++)
        {
            interpolatedPsiPtr[donorI] = pDW[donorI] * psiPtr[pDC[donorI]];
        }
    }

    // Collect influence from cells in own domain
    {
        const label nItems = ownInterpolationItems_.size();

        const interpolationItem* const __restrict__ oII =
            ownInterpolationItems_.begin();

        for(label itemI = 0; itemI<nItems; itemI++)
        {
            const interpolationItem& curItem = oII[itemI];
            interpolatedPsiPtr[curItem.donorID()] +=
                curItem.weight() * psiPtr[curItem.cellID()];
        }
    }

    // Send to neighbour procs
    // Number of this proc
    const label myProcNo = Pstream::myProcNo();

    if(Pstream::parRun())
    {
        const label nReq = UPstream::nRequests();
        const label tag = UPstream::msgType();

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                const label nItems = neighbourInterpolationItems_[procI].size();

                const interpolationItem* const __restrict__ nII =
                neighbourInterpolationItems_[procI].begin();

                scalarList& sendBuf = neighbourInterpolationBuf_[procI];
                {
                    scalar __restrict__ *pIP = sendBuf.begin();

                    for(label itemI = 0; itemI<nItems; itemI++)
                    {
                        const interpolationItem& curItem = nII[itemI];
                        pIP[itemI] =
                        curItem.weight() * psiPtr[curItem.cellID()];
                    }
                }

                scalarList& receiveBuf = neighbourValueBuf_[procI];

                IPstream::read
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<char*>(receiveBuf.begin()),
                    receiveBuf.byteSize(),
                    tag,
                    0
                );

                OPstream::write
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<const char*>(sendBuf.begin()),
                    sendBuf.byteSize(),
                    tag,
                    0
                );
            }
        }

        Pstream::waitRequests(nReq);

        forAll(neighbourValueToFieldMap_,procI)
        {
            if(procI != myProcNo)
            {
                const scalar* const __restrict__ pFP =
                    neighbourValueBuf_[procI].begin();

                const scalar nValues = neighbourValueToFieldMap_[procI].size();

                const label* const __restrict__ nVM =
                    neighbourValueToFieldMap_[procI].begin();
                for(label valueI = 0; valueI < nValues; valueI++)
                {
                    interpolatedPsiPtr[nVM[valueI]] += pFP[valueI];
                }
            }
        }
    }

    if(Pstream::parRun())
    {
        const label nReq = UPstream::nRequests();
        const label tag = UPstream::msgType();

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                const labelList& curDonorCols = donorCols_[procI];
                const label* const __restrict__ curDonorColsPtr =
                curDonorCols.begin();
                register const label nDonors = curDonorCols.size();

                scalarList& sendBuf = donorColBuf_[procI];
                {
                    scalar __restrict__ *buf = sendBuf.begin();
                    for(register label donorI = 0; donorI < nDonors; donorI++)
                    {
                        buf[donorI] = interpolatedPsiPtr[curDonorColsPtr[donorI]];
                    }
                }

                scalarList& receiveBuf = acceptorRowBuf_[procI];

                IPstream::read
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<char*>(receiveBuf.begin()),
                    receiveBuf.byteSize(),
                    tag,
                    0
                );

                OPstream::write
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<const char*>(sendBuf.begin()),
                    sendBuf.byteSize(),
                    tag,
                    0
                );
            }
        }

        Pstream::waitRequests(nReq);

        forAll(acceptorRows_,procI)
        {
            if(procI != myProcNo)
            {
                const labelList& curAcceptorRows = acceptorRows_[procI];
                const label* const __restrict__ curAcceptorRowPtr =
                    curAcceptorRows.begin();
                register const label nAcceptors = curAcceptorRows.size();
                const scalar* const __restrict__ receiveBufPtr =
                    acceptorRowBuf_[procI].begin();
                for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
                {
                    register label row = curAcceptorRowPtr[acceptorI];
                    ApsiPtr[row] -=
                        receiveBufPtr[acceptorI] * diagPtr[row];
                }
            }
        }

        const labelList& curDonorCols = donorCols_[myProcNo];
        const labelList& curAcceptorRows = acceptorRows_[myProcNo];

        const label* const __restrict__ curDonorColsPtr =
            curDonorCols.begin();
        const label* const __restrict__ curAcceptorRowPtr =
            curAcceptorRows.begin();

        register const label nAcceptors = curAcceptorRows.size();

        for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
        {
            register const label row = curAcceptorRowPtr[acceptorI];
            ApsiPtr[row] -=
                interpolatedPsiPtr[curDonorColsPtr[acceptorI]] * diagPtr[row];

        }
    }
    else
    {
        // Utilizing that donors are presorted during serial runs

        const labelList& curAcceptorRows = acceptorRows_[myProcNo];

        const label* const __restrict__ curAcceptorRowPtr =
            curAcceptorRows.begin();

        register const label nAcceptors = curAcceptorRows.size();

        for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
        {
            register const label row = curAcceptorRowPtr[acceptorI];
            ApsiPtr[row] -=
                interpolatedPsiPtr[acceptorI] * diagPtr[row];
        }
    }

    tpsi.clear();
}


void Foam::bellerophonLduMatrix::Tmul
(
    scalarField& Tpsi,
    const tmp<scalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    scalar __restrict__ * TpsiPtr = Tpsi.begin();

    const scalarField& psi = tpsi();
    const scalar* const __restrict__ diagPtr = diag().begin();

    scalar* __restrict__ psiPtr = const_cast<scalarField&>(psi).begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ lowerPtr = modLower().begin();
    const scalar* const __restrict__ upperPtr = modUpper().begin();

    // Doing nasty stuff with the internal field to remove influence of
    // acceptor cells on neighbour on other procs...
    autoPtr<scalarField> psiCache;
    const labelList& acceptors = bellerophon::Interpolation().acceptorCells();
    const label nAcceptors = acceptors.size();
    if(nAcceptors > 0)
    {
        //TODO: make this list persistent for the live time of the matrix
        psiCache.set(new scalarField(nAcceptors));
        scalar* __restrict__ cachePtr = psiCache().begin();
        const label* const __restrict__ acceptorRowPtr = acceptors.begin();
        for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
        {
            cachePtr[acceptorI] = psiPtr[acceptorRowPtr[acceptorI]];
            psiPtr[acceptorRowPtr[acceptorI]] = 0.0;
        }
    }

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    // Restore cached values
    if(nAcceptors > 0)
    {
        scalar* __restrict__ cachePtr = psiCache().begin();
        const label* const __restrict__ acceptorRowPtr = acceptors.begin();
        for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
        {
            psiPtr[acceptorRowPtr[acceptorI]] = cachePtr[acceptorI];
        }
    }

    register const label nCells = diag().size();
    for (register label cell = 0; cell<nCells; cell++)
    {
        TpsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }

    register const label nFaces = modUpper().size();
    for (register label face = 0; face<nFaces; face++)
    {
        TpsiPtr[uPtr[face]] += upperPtr[face]*psiPtr[lPtr[face]];
        TpsiPtr[lPtr[face]] += lowerPtr[face]*psiPtr[uPtr[face]];
    }

    // Label of the current proc
    const label myProcNo = Pstream::myProcNo();

    scalar __restrict__ * interpolatedPsiPtr = interpolatedPsi_.begin();

    // Adding local items to interpolated Psi
    {
        register const label nLocalAcceptors = acceptorRows_[myProcNo].size();

        const label* const __restrict__ donorCols =
            donorCols_[myProcNo].begin();

        const label* const __restrict__ acceptorRows =
            acceptorRows_[myProcNo].begin();

        for(register label acceptorI = 0; acceptorI < nLocalAcceptors; acceptorI++)
        {
            register const label row = acceptorRows[acceptorI];
            interpolatedPsiPtr[donorCols[acceptorI]] =
                - psiPtr[row] * diagPtr[row];
        }
    }

    // Send Psi to other primary donor proc
    if(Pstream::parRun())
    {
        const label nReq = UPstream::nRequests();
        const label tag = UPstream::msgType();

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                const labelList& curAcceptorRows = acceptorRows_[procI];
                const label* const __restrict__ curAcceptorRowPtr =
                curAcceptorRows.begin();
                register const label nAcceptors = curAcceptorRows.size();

                scalarList& sendBuf = acceptorRowBuf_[procI];
                {
                    scalar __restrict__ *buf = sendBuf.begin();
                    for
                    (
                        register label acceptorI = 0;
                        acceptorI < nAcceptors;
                        acceptorI++
                    )
                    {
                        register const label row = curAcceptorRowPtr[acceptorI];
                        buf[acceptorI] = -diagPtr[row] * psiPtr[row];
                    }
                }

                scalarList& receiveBuf = donorColBuf_[procI];

                IPstream::read
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<char*>(receiveBuf.begin()),
                    receiveBuf.byteSize(),
                    tag,
                    0
                );

                OPstream::write
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<const char*>(sendBuf.begin()),
                    sendBuf.byteSize(),
                    tag,
                    0
                );
            }
        }

        Pstream::waitRequests(nReq);

        forAll(donorCols_,procI)
        {
            if(procI != myProcNo)
            {
                const labelList& curDonorCols = donorCols_[procI];

                const label* const __restrict__ curDonorColPtr =
                    curDonorCols.begin();

                register const label nDonors = curDonorCols.size();

                const scalar* const __restrict__ buf =
                    donorColBuf_[procI].begin();

                for
                (
                    register label donorI = 0;
                    donorI < nDonors;
                    donorI++
                )
                {
                    interpolatedPsiPtr[curDonorColPtr[donorI]] = buf[donorI];
                }
            }
        }
    }

    // Add primary donor component to Tpsi
    {
        const label* const __restrict__ pDC =
            primaryDonorCells_.begin();

        const scalar* const __restrict__ pDW =
            primaryDonorWeights_.begin();

        const label nDonors = primaryDonorCells_.size();

        for(label donorI = 0; donorI<nDonors; donorI++)
        {
            TpsiPtr[pDC[donorI]] +=
                pDW[donorI] * interpolatedPsiPtr[donorI];
        }
    }

    // Add component of own interpolation items to Tpsi
    {
        const scalar nItems = ownInterpolationItems_.size();
        const interpolationItem* const __restrict__ oII =
            ownInterpolationItems_.begin();

        for(label itemI = 0; itemI < nItems; itemI++)
        {
            const interpolationItem& curItem = oII[itemI];
            TpsiPtr[curItem.cellID()] +=
                curItem.weight() * interpolatedPsiPtr[curItem.donorID()];
        }
    }

    // Send and receive neighbour interpolation items to Tpsi
    {
        const label nReq = UPstream::nRequests();
        const label tag = UPstream::msgType();

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                const label* const __restrict__ nPM =
                    neighbourValueToFieldMap_[procI].begin();
                register const label nValues =
                    neighbourValueToFieldMap_[procI].size();

                scalarList& sendBuf = neighbourValueBuf_[procI];

                {
                    scalar __restrict__ *buf = sendBuf.begin();
                    for
                        (
                            register label valueI = 0;
                            valueI < nValues;
                            valueI++
                        )
                        {
                            buf[valueI] = interpolatedPsiPtr[nPM[valueI]];
                        }
                }

                scalarList& receiveBuf = neighbourInterpolationBuf_[procI];

                IPstream::read
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<char*>(receiveBuf.begin()),
                    receiveBuf.byteSize(),
                    tag,
                    0
                );

                OPstream::write
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<const char*>(sendBuf.begin()),
                    sendBuf.byteSize(),
                    tag,
                    0
                );
            }
        }

        Pstream::waitRequests(nReq);
        forAll(donorCols_,procI)
        {
            if(procI != myProcNo)
            {
                const interpolationItem* const __restrict__ nII =
                    neighbourInterpolationItems_[procI].begin();

                register const label nItems =
                    neighbourInterpolationItems_[procI].size();

                const scalar* const __restrict__ buf =
                    neighbourInterpolationBuf_[procI].begin();

                for
                (
                    register label itemI = 0;
                    itemI < nItems;
                    itemI++
                )
                {
                    const interpolationItem& curItem = nII[itemI];
                    TpsiPtr[curItem.cellID()] += curItem.weight() * buf[itemI];
                }
            }
        }
    }

    // Repeat caching if there are cyclic patches TODO if(cyclics)...
    if(nAcceptors > 0)
    {
        const label* const __restrict__ acceptorRowPtr = acceptors.begin();
        for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
        {
            psiPtr[acceptorRowPtr[acceptorI]] = 0.0;
        }
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    // Restore cached values
    if(nAcceptors > 0)
    {
        scalar* __restrict__ cachePtr = psiCache().begin();
        const label* const __restrict__ acceptorRowPtr = acceptors.begin();
        for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
        {
            psiPtr[acceptorRowPtr[acceptorI]] = cachePtr[acceptorI];
        }
        psiCache.clear();
    }

    tpsi.clear();
}

void Foam::bellerophonLduMatrix::sumA
(
    scalarField& sumA,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
) const
{
    lduMatrix::sumA(sumA, interfaceBouCoeffs, interfaces);
    const labelList& acceptorRows = acceptorRows_[Pstream::myProcNo()];

    scalar* __restrict__ sumAPtr = sumA.begin();
    const scalar* const __restrict__ diagPtr = diag().begin();
    const label* const __restrict__ rowPtr = acceptorRows.begin();

    const label nRows = acceptorRows.size();
    for(register label rowI = 0; rowI < nRows; rowI++)
    {
        register const label row = rowPtr[rowI];
        sumAPtr[row] -= diagPtr[row];
    }
}


void Foam::bellerophonLduMatrix::residual
(
    scalarField& rA,
    const scalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    // TODO: Check if this is implemented correct

    // Correct source
    tmp<scalarField> tCorrSource =
        bellerophon::Interpolation().correctSource(source, diag(), cmpt);

    scalarField& source2 = tCorrSource();

    scalar* __restrict__ rAPtr = rA.begin();

    const scalar* const __restrict__ psiPtr = psi.begin();
    const scalar* const __restrict__ diagPtr = diag().begin();
    const scalar* const __restrict__ sourcePtr = source2.begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = modUpper().begin();
    const scalar* const __restrict__ lowerPtr = modLower().begin();

    // Parallel boundary initialisation.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    FieldField<Field, scalar> mBouCoeffs(interfaceBouCoeffs.size());

    forAll(mBouCoeffs, patchi)
    {
        if (interfaces.set(patchi))
        {
            mBouCoeffs.set(patchi, -interfaceBouCoeffs[patchi]);
        }
    }

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        mBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt
    );

    register const label nCells = diag().size();
    for (register label cell=0; cell<nCells; cell++)
    {
        rAPtr[cell] = sourcePtr[cell] - diagPtr[cell]*psiPtr[cell];
    }


    register const label nFaces = modUpper().size();

    for (register label face=0; face<nFaces; face++)
    {
        rAPtr[uPtr[face]] -= lowerPtr[face]*psiPtr[lPtr[face]];
        rAPtr[lPtr[face]] -= upperPtr[face]*psiPtr[uPtr[face]];
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        mBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt
    );

    const label myProcNo = Pstream::myProcNo();

    if(Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::nonBlocking);

        forAll(donorCols_,procI)
        {
            if(procI != myProcNo)
            {
                UOPstream toDomain(procI, pBufs);
                const labelList& curDonorCols = donorCols_[procI];
                const label* const __restrict__ curDonorColsPtr =
                    curDonorCols.begin();
                register const label nDonors = curDonorCols.size();
                for(register label donorI = 0; donorI < nDonors; donorI++)
                {
                    toDomain << psiPtr[curDonorColsPtr[donorI]];
                }
            }
        }

        pBufs.finishedSends();

        forAll(acceptorRows_,procI)
        {
            if(procI != myProcNo)
            {
                UIPstream str(procI, pBufs);
                const labelList& curAcceptorRows = acceptorRows_[procI];
                const label* const __restrict__ curAcceptorRowPtr =
                    curAcceptorRows.begin();
                register const label nAcceptors = curAcceptorRows.size();
                for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
                {
                    register label row = curAcceptorRowPtr[acceptorI];
                    scalar received;
                    str >> received;
                    rAPtr[row] -=
                        received * diagPtr[row];
                }
            }
        }
    }

    const labelList& curDonorCols = donorCols_[myProcNo];
    const labelList& curAcceptorRows = acceptorRows_[myProcNo];

    const label* const __restrict__ curDonorColsPtr =
        curDonorCols.begin();
    const label* const __restrict__ curAcceptorRowPtr =
        curAcceptorRows.begin();

    register const label nAcceptors = curAcceptorRows.size();

    for(register label acceptorI = 0; acceptorI < nAcceptors; acceptorI++)
    {
        register const label row = curAcceptorRowPtr[acceptorI];
        rAPtr[row] -=
            psiPtr[curDonorColsPtr[acceptorI]] * diagPtr[row];
    }

}


Foam::tmp<Foam::scalarField> Foam::bellerophonLduMatrix::residual
(
    const scalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    tmp<scalarField> trA(new scalarField(psi.size()));
    residual(trA(), psi, source, interfaceBouCoeffs, interfaces, cmpt);
    return trA;
}


// ************************************************************************* //
