#pragma once
  
// Tags
//
// Note: 
//    Naming should be Name{} name;
//
// Usage: 
//    Define functions that depend on tags as
//    func(xx, yy, ..., Mode::Tag_name)
//    and call the function as
//    yy = func(aa, bb,..., Mode::tag_name)
//

namespace corgi { namespace tags {

/// I/O Write modes 
class Write_mode
{
  public:

    /// Standard/default writing mode
    static struct Standard{} standard;



};



}} // ns corgi::tags
