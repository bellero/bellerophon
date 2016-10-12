/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef searchItem_H
#define searchItem_H

#include "Ostream.H"
#include "Istream.H"
#include "point.H"

namespace Foam
{

class searchItem;
inline Ostream& operator<<(Ostream&, const searchItem&);
inline Istream& operator>>(Istream&, searchItem&);
inline bool operator==(const searchItem&, const searchItem&);
inline bool operator!=(const searchItem&, const searchItem&);

class searchItem
{
private:
    //- ID of the processor the cell containing the point is on
    label procID_;

    //- Label of the cell the point is in
    label cellLabel_;

    //- Position of the point in question
    point position_;

    //- Cell in focus, after search: cell the point is in
    label seed_;

    //- Zone the point should be in
    label zoneID_;

    //- Group the item belongs to, i.e. processor the search was started from
    label groupID_;

    //- Running number in the group
    label groupIndex_;

public:
    //- Empty Constructor
    searchItem()
    :
    procID_(-1),
    cellLabel_(-1),
    position_(vector::zero),
    seed_(0),
    zoneID_(-1),
    groupID_(-1)
    {};

    //- Construct from proc ID and cell label
    searchItem
    (
        const label proc,
        const label cell,
        const point& pos,
        const label seed = 0,
        const label zoneID = -1,
        const label groupID = -1,
        const label groupIndex = -1
    )
    :
    procID_(proc),
    cellLabel_(cell),
    position_(pos),
    seed_(seed),
    zoneID_(zoneID),
    groupID_(groupID),
    groupIndex_(groupIndex)
    {};

    //- Destructor
    ~searchItem() {};

    //- Constant access to procID
    label procID() const { return procID_; };

    //- Non-constant access to procID
    label& procID() { return procID_; };

    //- Constant access to cell label
    label cellLabel() const { return cellLabel_; };

    //- Non-constant access to cell label
    label& cellLabel() { return cellLabel_; };

    //- Constant access to point position
    const point& position() const { return position_; };

    //- Non-constant access to point position
    point& position() { return position_; };

    //- Constant access to seed label
    label seed() const { return seed_; };

    //- Non-constant access to seed label
    label& seed() { return seed_; };

    //- Constant access to zone ID
    label zoneID() const { return zoneID_; };

    //- Non-constant access to zone ID
    label& zoneID() { return zoneID_; };

    //- Constant access to zone ID
    label groupID() const { return groupID_; };

    //- Non-constant access to zone ID
    label& groupID() { return groupID_; };

    //- Constant access to zone ID
    label groupIndex() const { return groupIndex_; };

    //- Non-constant access to zone ID
    label& groupIndex() { return groupIndex_; };

    friend Ostream& operator<<(Ostream&, const searchItem&);

    friend Istream& operator>>(Istream&, searchItem&);

    friend bool operator==(const searchItem&, const searchItem&);

    friend bool operator!=(const searchItem&, const searchItem&);
};

inline bool operator==(const searchItem& a, const searchItem& b)
{
    return ( a.procID_    == b.procID_    )
           &&
           ( a.cellLabel_ == b.cellLabel_ )
           &&
           ( a.position_  == b.position_  )
           &&
           ( a.seed_      == b.seed_      )
           &&
           ( a.zoneID_    == b.zoneID_    )
           &&
           ( a.groupID_   == b.groupID_   )
           &&
           ( a.groupIndex_== b.groupIndex_);
}

inline bool operator!=(const searchItem& a, const searchItem& b)
{
    return ( a.procID_    != b.procID_    )
           ||
           ( a.cellLabel_ != b.cellLabel_ )
           ||
           ( a.position_  != b.position_  )
           ||
           ( a.seed_      != b.seed_      )
           ||
           ( a.zoneID_    != b.zoneID_    )
           ||
           ( a.groupIndex_!= b.groupIndex_);
}

inline Ostream& operator<<
(
    Ostream& os,
    const searchItem& result
)
{
    os  <<token::BEGIN_LIST
        <<result.procID_ << token::SPACE
        <<result.cellLabel_ << token::SPACE
        <<result.position_ << token::SPACE
        <<result.seed_ << token::SPACE
        <<result.zoneID_ << token::SPACE
        <<result.groupID_ << token::SPACE
        <<result.groupIndex_
        << token::END_LIST;
    return os;
}

inline Istream& operator>>(Istream& is, searchItem& result)
{
    is.readBegin("searchItem");
    is>>result.procID_
      >>result.cellLabel_
      >>result.position_
      >>result.seed_
      >>result.zoneID_
      >>result.groupID_
      >>result.groupIndex_;
    is.readEnd("searchItem");

    is.check("operator>>(Istream&, searchItem&)");

    return is;
}

} // End namespace Foam

#endif // searchItem_H