/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef gradientSearch_H
#define gradientSearch_H

#include "searchItem.H"
#include "fvMesh.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

namespace Foam
{


class gradientSearch
{
private:
    //- Mesh to search on
    const fvMesh& mesh_;

    const labelList& owner_;

    const labelList& neighbour_;

    //- Local labels of the cells behind the processor patches for each proc
    labelListList neighbourFaceCells_;

    // List of ten last cells
    mutable labelList visited_;

    //- Cell trees for cell zones and overall 
    mutable List< autoPtr< indexedOctree<treeDataCell> > > cellTreePtrs_;

    //- Boundboxes for zones
    mutable List< autoPtr<boundBox> > zoneBoundPtrs_;

    mutable autoPtr< List<searchItem> > failItemsPtr_;

    mutable autoPtr< List<boundBox> > bBoxesPtr_;

    //- Disable empty contructor
    gradientSearch();

    //- Disable copy contructor
    gradientSearch(const gradientSearch&);

    //- Disable copy-assignment operator
    gradientSearch& operator=(const gradientSearch&);

    // Search for a searchItem in local mesh and with seed cell and reference
    // to the next proc to search on
    label searchLocal(searchItem&, label) const;

    // Get possible seed cells for a given zone ID
    labelList zoneSeeds(const label) const;

    bool visited(const label) const;

    void setVisited(const label) const;

    void clearVisited() const;

    // Collect global succes for items
    void collectGlobalSuccess
    (
        labelListList& finallyFoundPerProc
    ) const;

public:
    //- Construct from mesh
    gradientSearch(const fvMesh&);

    ~gradientSearch();

    // making these non-virtual...
//     TypeName("gradientSearch");
    ClassName("gradientSearch");
    const word& type() const { return typeName; };

    //- Search for cell containing point
    autoPtr< searchItem > search(const point&, const label, const label) const;

    //- Search for cells containing points without seed cells
    autoPtr< List< searchItem > > search(const List<point>&) const;

    //- Search for cells containing points without seed cells
    autoPtr< List< searchItem> > search
    (
        const List<point>&,
        const List<label>&,
        const List<label>&
    ) const;

    //- Search for autoPtr to List of searchItems
    void search ( autoPtr< List< searchItem > >& ) const;

    autoPtr< List<searchItem> > failItems() { return failItemsPtr_; };

    const indexedOctree<treeDataCell>& cellTree(const label zoneID) const;

    const boundBox& zoneBounds(const label zoneID) const;
};

} // End namespace Foam

#endif // gradientSearch_H