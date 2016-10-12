/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef bellerophonInterpolation_H
#define bellerophonInterpolation_H

#include "typeInfo.H"
#include "Pstream.H"
#include "labelList.H"
#include "scalarList.H"


#include "primitiveFieldsFwd.H"
#include "volFields.H"
#include "searchItem.H"
#include "interpolationItem.H"
#include "bellerophonInterpolationMethod.H"

namespace Foam
{

//- Forward definition of some classes
class bellerophonInterface;
class bellerophonInterfaceField;
class bellerophonLduMatrix;

class bellerophon
{
private:
    //- Pointers to relevant objects

        //- Pointer to dictionary containing settings for interpolation
        mutable autoPtr<IOdictionary> dictPtr_;

        //- Pointer to fvMesh
        mutable fvMesh const* meshPtr_;

        //- Pointer to patches
        mutable List<bellerophonInterface*> patchPtrs_;

        //- Constant pointers to patches for returning to callers
        mutable List<const bellerophonInterface*> cifs_;

        //- Total number of Acceptor cells
        mutable label nAcceptors_;

        //- List of Acceptor cells
        mutable labelList acceptorCells_;

        //- Total number of Donor cells
        mutable label nDonors_;

        //- List of primary Donor cells
        //  Used as seed cells for donor search if present
        mutable labelList primaryDonorCells_;

        //- Implicit hole cells
        mutable labelList holeCells_;

        //- Markup for new live cells
        mutable boolList newLiveCellsMarkup_;

        //- Live cells that haven't been live cells in the previous time step
        mutable labelList newLiveCells_;

        // Hold a mapping from acceptors to donors if their order was changed
        // during gradient search (only needed for serial runs)
        mutable autoPtr<labelList> acceptorMapPtr_;

        //- Vector from primary donor cell to acceptor cell
        mutable vectorField deltas_;

        //- Interpolation method holding interpolation data
        mutable autoPtr<bellerophonInterpolationMethod> interpolationMethod_;

        //- Seed cells per zone id
        mutable labelList zoneSeeds_;

        //- Zone IDs for each acceptor
        mutable labelList zoneIDs_;

        //- Faces between acceptor and regular cells
        mutable labelListList interfaceFaces_;

        //- Flip map for interface faces (true if acceptor cell is owner)
        mutable boolListList interfaceFlipMap_;

        //- State of the interpolation information
        mutable bool topologyUpToDate_;
        mutable bool interpolationUpToDate_;

        //- Donor Items on this proc
        mutable autoPtr< List<searchItem> > donorItemsPtr_;

        //- Name of field to enforce continuity on
        mutable List<word> continuityFields_;

        //- Name of field to force zero in first live cells
        mutable List<word> forceZeroFields_;

    //- Hide common construction and copy functions
        //- Disallow Constructor
        bellerophon();

        //- Disallow Copy-Constructor
        bellerophon(const bellerophon&);

        //- Disallow another copy
        void operator=(const bellerophon&);

        //- Private destructor
        ~bellerophon();

public:

    TypeName("bellerophon");

    enum cellState
    {
        REGULAR_CELL = 0,            // Regular fluid cell, value from equations
        DONOR_CELL = 1,              // Donor cell, value from equations
        EXTENDED_ACCEPTOR_CELL = 2,  // Extended acceptor cell, next to acceptor,
                                     // may not become primary donor cell as it
                                     // is neighbour to acceptor cell
        ACCEPTOR_CELL = 3,           // Acceptor cell, value is interpolated
        HOLE_CELL = 4                // Hole cell, not involved in solution
    };


    //- Return the Interpolation
    static bellerophon& Interpolation();

    //- Update functions

        //- Update the topological properties
        //  Updates number and labels of acceptor cells
        void updateTopology() const;

        //- Mark old acceptor and hole cells before clearing topology
        //  Must be called separately
        void markHole() const;

        //- Add hole from overset domain from map
        label addHole
        (
            const labelList& holeMap,
            const label donorZoneID,
            const label oversetPatchID
        ) const;

        //- Update the Interpolation
        //  Updates the primary donor cells
        void updateInterpolation() const;


        //- Clear topological information
        void clearTopology() const;


        //- Clear interpolation information
        void clearInterpolation() const;


    //- Data access
        //- Return fvMesh
        const fvMesh& mesh() const;

        //- Return dictionary
        const dictionary& dict() const;

        //- Return pointers to interfaces
        const List<const bellerophonInterface*> interfaces() const
        {
            if(cifs_.size() != patchPtrs_.size())
            {
                cifs_.setSize(patchPtrs_.size());
                forAll(patchPtrs_,patchI)
                {
                    cifs_[patchI] = patchPtrs_[patchI];
                }
            }
            return cifs_;
        }

        //- Acceptor cell indices
        const labelList& acceptorCells() const
        {
            this->updateTopology();
            return acceptorCells_;
        };

        //- Number of Donor cells
        label nDonors() const
        {
            updateInterpolation();
            return nDonors_;
        };

        //- Primary donor cell indices
        const labelList& primaryDonorCells() const
        {
            updateInterpolation();
            return primaryDonorCells_;
        };

        //- Hole cell indices
        const labelList& holeCells() const
        {
            return holeCells_;
        };

        //- New live cell indices
        const labelList& newLiveCells() const
        {
            return newLiveCells_;
        };

        //- Primary donor weights
        const labelList& zoneSeeds() const
        {
            updateInterpolation();
            return zoneSeeds_;
        };

        //- Interface faces for given patch
        const labelList& interfaceFaces(const label patchI)
        {
            updateTopology();
            return interfaceFaces_[patchI];
        }

        //- Flipmap for interface faces, true if acceptor cell is owner
        const boolList& interfaceFlipMap(const label patchI)
        {
            updateTopology();
            return interfaceFlipMap_[patchI];
        }

        //- Primary donor weights
        const scalarField& primaryDonorWeights() const
        {
            updateInterpolation();
            return interpolationMethod_().primaryDonorWeights();
        };

        //- Interpolation items in own domain
        const List<interpolationItem>& ownInterpolationItems() const
        {
            updateInterpolation();
            return interpolationMethod_().ownInterpolationItems();
        };

        //- Interpolation items for neighbour domains
        const List< List<interpolationItem> >& neighbourInterpolationItems() const
        {
            updateInterpolation();
            return interpolationMethod_().neighbourInterpolationItems();
        };

        //- Maps interpolation contribution of interpolation items on
        //  neighbour proc to interpolated field
        const labelListList& neighbourValueToFieldMap() const
        {
            updateInterpolation();
            return interpolationMethod_().neighbourValueToFieldMap();
        }

        //- Index of the donor on the local domain stored separately for
        //  each domain the acceptor cell is in
        const labelListList& donorCols() const
        {
            updateInterpolation();
            return interpolationMethod_().donorCols();
        }

        //- Rows of the acceptor cells on the local domain stored separately for
        //  each domain the donor cell is in
        const labelListList& acceptorRows() const
        {
            updateInterpolation();
            return interpolationMethod_().acceptorRows();
        }

        //- Enforce continuity on given field?
        bool enforceContinuity(const word&) const;

        //- Enforce continuity on given field?
        bool forceZero(const word&) const;

    //- Information Collection

        //- Add a patch interface
        label addInterface(bellerophonInterface&);

        //- Remove a patch interface from list
        void removeInterface(bellerophonInterface&);

    //- Overset interpolation
        //- Interpolate field to acceptor cells
        template<class Type>
        void interpolate
        (
            const word& fieldName,
            Pstream::commsTypes commsType,
            bool bound = false
        ) const;

    //- Matrix manipulation
        //- Add additional coefficients
        void updateMatrix
        (
            bellerophonLduMatrix&,
            const bellerophonInterfaceField*
        ) const;

        //- Modifiy RHS of matrix for initial residuals
        tmp<scalarField> correctSource
        (
            const scalarField&,
            const scalarField&,
            const direction
        ) const;
};

} // End namespace Foam

#endif // bellerophonInterpolation_H