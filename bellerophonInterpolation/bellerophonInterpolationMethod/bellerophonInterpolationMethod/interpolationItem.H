/*
 * TODO: add funky header and license here...
 */

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef interpolationItem_H
#define interpolationItem_H

#include "scalar.H"
#include "OFstream.H"
#include "IFstream.H"

namespace Foam
{

//- Interpolation data for matrix multiplications
class interpolationItem
{
private:
    // Label of donor item, i.e. position in interpolated field
    label donorID_;

    // Label of cell containing value, i.e. donor cell
    label cellID_;

    // Weight for interpolation
    scalar weight_;

public:
    // Empty default contructor
    interpolationItem()
    :
    donorID_(-1),
    cellID_(-1),
    weight_(0.0)
    {};

    // Construct from donor, cell and weight
    interpolationItem
    (
        const label d,
        const label c,
        const scalar w
    )
    :
    donorID_(d),
    cellID_(c),
    weight_(w)
    {};

    // Non-constant access to donorID
    inline label& donorID() { return donorID_; };

    // Constant access to donorID
    inline label donorID() const { return donorID_; };

    // Non-constant access to cellID
    inline label& cellID() { return cellID_; };

    // Constant access to cellID
    inline label cellID() const { return cellID_; };

    // Non-constant access to weight
    inline scalar& weight() { return weight_; };

    // Constant access to weight
    inline scalar weight() const { return weight_; };

    // Set all three attributes at once
    inline void set(const label d, const label c, const scalar& w)
    {
        donorID_ = d;
        cellID_  = c;
        weight_  = w;
    }

    friend Ostream& operator<<(Ostream&, const interpolationItem&);

    friend Istream& operator>>(Istream&, interpolationItem&);

    friend bool operator==(const interpolationItem&, const interpolationItem&);

    friend bool operator!=(const interpolationItem&, const interpolationItem&);
};

inline bool operator==(const interpolationItem& a, const interpolationItem& b)
{
    return ( a.donorID_ == b.donorID_ )
           &&
           ( a.cellID_  == b.cellID_  )
           &&
           ( a.weight_  == b.weight_  );
}

inline bool operator!=(const interpolationItem& a, const interpolationItem& b)
{
    return ( a.donorID_ != b.donorID_ )
           ||
           ( a.cellID_  != b.cellID_  )
           ||
           ( a.weight_  != b.weight_  );
}

inline Ostream& operator<<
(
    Ostream& os,
    const interpolationItem& result
)
{
    os  <<token::BEGIN_LIST
        <<result.donorID_ << token::SPACE
        <<result.cellID_ << token::SPACE
        <<result.weight_
        << token::END_LIST;
    return os;
}

inline Istream& operator>>(Istream& is, interpolationItem& result)
{
    is.readBegin("interpolationItem");
    is>>result.donorID_
      >>result.cellID_
      >>result.weight_;
    is.readEnd("interpolationItem");

    is.check("operator>>(Istream&, interpolationItem&)");

    return is;
}

} // End namespace Foam

#endif // interpolationItem_H