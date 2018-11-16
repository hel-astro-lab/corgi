#pragma once
  
// Tags
//
// Note: 
//    Naming should be Name{} name;
//
// Usage: 
//    Define functions that depend on tags as
//    func(xx, yy, ..., Mode::TagName)
//    and call the function as
//    yy = func(aa, bb,..., Mode::tagName)
//

namespace corgi { namespace tags {

/// I/O Write modes 
class WriteMode
{
  public:

    /// Standard/default writing mode
    static struct Standard{} standard;



};



}} // ns corgi::tags
