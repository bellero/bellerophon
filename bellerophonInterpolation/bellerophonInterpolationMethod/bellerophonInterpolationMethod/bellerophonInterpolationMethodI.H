/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef bellerophonInterpolationMethodI_H
#define bellerophonInterpolationMethodI_H

#include "bellerophonInterpolationMethod.H"

//- Access to primary donor weights
const Foam::scalarField&
Foam::bellerophonInterpolationMethod::primaryDonorWeights() const
{
    return primaryDonorWeightsPtr_();
};

//- Interpolation items in own domain
const Foam::List<Foam::interpolationItem>&
Foam::bellerophonInterpolationMethod::ownInterpolationItems() const
{
    return ownInterpolationItemsPtr_();
};

//- Interpolation items for neighbour domains
const Foam::List< Foam::List<Foam::interpolationItem> >& 
Foam::bellerophonInterpolationMethod::neighbourInterpolationItems() const
{
    return neighbourInterpolationItemsPtr_();
};

//- Maps interpolation contribution of interpolation items on
//  neighbour proc to interpolated field
const Foam::labelListList&
Foam::bellerophonInterpolationMethod::neighbourValueToFieldMap() const
{
    return neighbourValueToFieldMapPtr_();
}

//- Index of the donor on the local domain stored separatly for
//  each domain the acceptor cell is in
const Foam::labelListList&
Foam::bellerophonInterpolationMethod::donorCols() const
{
    return donorColPtr_();
}

//- Rows of the acceptor cells on the local domain stored separatly for
//  each domain the donor cell is in
const Foam::labelListList&
Foam::bellerophonInterpolationMethod::acceptorRows() const
{
    return acceptorRowPtr_();
}

#endif // bellerophonInterpolationMethodI_H